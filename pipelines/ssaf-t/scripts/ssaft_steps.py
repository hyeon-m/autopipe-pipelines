"""The SSAF-T stages after load: zstack, tile, panel_vf, downsample, cluster, map, filter, segment
(annotate is annotate.py). Each takes (run, unit), works on units/<unit>/, and returns the short
summary run_one records in state.json.

  zstack      the image a 3D map's downsample samples: the z-stack (only the planes chosen in
              z_planes, each at its own number), or for a 3D map on a 2D image that 2D image
  tile        the section cut into tiles (Xenium, MERSCOPE) with SSAF-T's slice_spatialdata.py
  panel_vf    count platforms: gene panel (ontology markers + HVGs) and the blurred vector field
  downsample  anchors sampled inside the Otsu tissue mask of the image chosen (image_key)
  cluster     PCA, kNN graph and Leiden on the anchors (GPU when attached)
  map         the annotated anchors mapped back over the whole tissue
  filter      the density cut, the cell-type map written into the SpatialData
  segment     cells cut from the map, seeded by nuclei
"""
from __future__ import annotations

import contextlib
import glob
import re
import subprocess
import sys
import time
from pathlib import Path

from ssaft_core import analysis, device, image_key_2d, log, persist, read_sdata, scale0, transcripts_key
from ssaft_load import ZKEY

ZSRC = "ssaft_zstack_src"      # the in-memory element a 3D downsample samples


# ------------------------------------------------------------------------------------------- z-stack
def chosen_planes(run):
    p = run.cfg.get("zstack", {}).get("planes")
    return [int(v) for v in p] if p else None


def flat_image(run) -> bool:
    """A 3D map downsampled on a 2D image (no z-stack)."""
    return run.cfg["dimension"] == 3 and (run.cfg.get("downsample_image") or {}).get("kind") == "2d"


def zstack_files_element(run, sdata):
    """The z-stack file / glob (zstack.source) as a lazy 3D image element: nothing is copied, the
    planes are read from the original ome.tif / per-plane tifs when downsampling needs them. The
    planes chosen in config are passed to ssaft's zstack_image as z=, so it holds those planes only,
    each at its own number; with none chosen, every plane."""
    import ssaft
    Zc = run.cfg["zstack"]
    src = Zc["source"]
    source = sorted(glob.glob(src)) if "*" in str(src) else src
    kw = {}
    if Zc.get("z_step"):
        kw["z_step"] = float(Zc["z_step"])
    planes = chosen_planes(run)
    if planes:
        kw["z"] = planes
    log(f"z-stack: {'planes ' + str(planes) if planes else 'every plane'} of {src}")
    return ssaft.Analysis.zstack_image(source, like=sdata.images[image_key_2d(run, sdata)], **kw), src


def zstack_element(run, sdata):
    """The image a 3D map's downsample samples: the z-stack, or for a 3D map on a 2D image that 2D
    image (SSAF-T places it at the middle of the transcripts' z range)."""
    if flat_image(run):
        key = image_key_2d(run, sdata)
        log(f"3D map on the 2D image images[{key!r}]: one tissue mask, placed at the middle of the "
            f"transcripts' z range (no z-stack)")
        return sdata.images[key], f"images[{key!r}] (2D)"
    if ZKEY in sdata.images:
        return sdata.images[ZKEY], f"images[{ZKEY!r}]"
    return zstack_files_element(run, sdata)


def st_zstack(run, unit):
    """Check the z-stack and its placement, with the chosen planes. A tiled run writes the z-stack
    into the zarr so every tile carries it; otherwise it is attached lazily at downsample time."""
    sdata = read_sdata(run.udir(unit) / "spatialdata.zarr")
    copied = False
    if (run.cfg.get("zarr") or {}).get("zstack_element") and ZKEY in sdata.images:
        el, src = sdata.images[ZKEY], f"images[{ZKEY!r}] (from the input zarr)"
    else:
        if ZKEY in sdata.images:
            # written by an earlier attempt with other planes or tiling: made again below if needed
            with contextlib.suppress(Exception):
                sdata.delete_element_from_disk(ZKEY)
            del sdata.images[ZKEY]
        el, src = zstack_files_element(run, sdata)
        if "tile" in run.cfg["stages"]:
            sdata.images[ZKEY] = el
            persist(sdata, [ZKEY])
            copied = True
    shape = tuple(scale0(el).shape)
    return {"source": str(src), "shape_czyx": list(shape), "copied": copied,
            "planes": chosen_planes(run) or "all"}


# ---------------------------------------------------------------------------------------------- tile
def st_tile(run, unit):
    import ssaft
    cfg = run.cfg
    zarr = run.udir(unit) / "spatialdata.zarr"
    sdata = read_sdata(zarr)
    keys = [image_key_2d(run, sdata)] + ([ZKEY] if cfg["dimension"] == 3 else [])
    crops = run.udir(unit) / "crops"
    if crops.exists() and any(crops.glob("*.zarr")):
        # this stage never finished for this run (state says so), so the tiles are its own partial
        # output from an interrupted attempt: start the tiling over
        import shutil
        log(f"removing partial tiles from an interrupted attempt: {crops}")
        shutil.rmtree(crops)
    script = Path(ssaft.__file__).parent / "scripts" / "slice_spatialdata.py"
    cmd = [sys.executable, str(script), "--sdata-path", str(zarr), "--tile-size", str(cfg["tile"]["tile_size"]),
           "--output-dir", str(crops), "--image-key", *keys, "--keep-points", transcripts_key(sdata, cfg["platform"]),
           "--drop-labels", "--drop-shapes", "--drop-tables"]
    log(" ".join(cmd))
    r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(r.stdout[-20000:])
    if r.returncode != 0:
        raise RuntimeError(f"slice_spatialdata.py exited {r.returncode}")
    return {"tiles": len(list(crops.glob("*.zarr"))), "dir": str(crops)}


# ------------------------------------------------------------------------------------------ panel_vf
def st_panel_vf(run, unit):
    import ssaft
    cfg = run.cfg
    sdata = read_sdata(run.udir(unit) / "spatialdata.zarr")
    grid = ssaft.GridDataset.from_sdata(sdata)
    P, V = cfg["panel"], cfg["vf"]
    hvg_input = grid.raw
    if P.get("hvg_bin_um"):
        hvg_input = ssaft.aggregate_to_grid(grid.raw, bin_um=float(P["hvg_bin_um"]), origin=grid.origin)
    panel, onto = ssaft.select_gene_panel(hvg_input, cfg["species"], cfg["organ"], hvg_n=int(P["hvg_n"]),
                                          server_url=cfg.get("cellman_server"))
    (run.udir(unit) / "panel.txt").write_text("\n".join(panel))
    del hvg_input
    an = analysis(run, unit)
    an.compute_vf(sdata, bin_um=float(V["bin_um"]), panel=panel, sampling_distance=float(V["sampling_distance"]),
                  bandwidth_um=float(V["bandwidth_um"]), kernel=V["kernel"], device=device(run) or "cpu")
    return {"panel_genes": len(panel), "ontology_markers": len(onto), "bins": int(grid.raw.n_obs)}


# ---------------------------------------------------------------------------------------- downsample
def st_downsample(run, unit):
    """downsample_vector with the one image it samples as image_key= (SSAF-T 23dd889): the 2D image
    for a 2D map; for a 3D map the z-stack (a tiled run: the ssaft_zstack element every tile
    carries), or the 2D image when the 3D map is on a 2D image. otsu_mult on every platform."""
    cfg = run.cfg
    D = dict(cfg["downsample"])
    an = analysis(run, unit)
    zarr = run.udir(unit) / "spatialdata.zarr"
    gpu = device(run)
    common = dict(gpu_device=gpu or "cuda", auto_gpu=gpu is not None)
    if D.get("otsu_mult") is not None:
        common["otsu_mult"] = float(D["otsu_mult"])
    if cfg["regime"] == "count":
        sdata = read_sdata(zarr)
        adata = an.downsample_vector(sdata, z=2, n_samples_per_block=int(D["n_samples_per_block"]),
                                     block_um=float(D["block_um"]), vf_norm_percentile=D.get("vf_norm_percentile"),
                                     **common)
    else:
        img = cfg["image"]
        kw = dict(n_samples=int(D["n_samples"]), block_um=float(D["block_um"]), bandwidth=float(D["bandwidth"]),
                  sampling_distance=float(D["sampling_distance"]), kernel=D["kernel"],
                  revert_image=bool(img["revert_image"]), intensity=img.get("intensity", "luminance"), **common)
        chunked = "tile" in cfg["stages"]
        sdata = read_sdata(zarr)
        # name the molecules explicitly: CosMx units also carry the vendor cell centroids as points
        kw["element_name"] = transcripts_key(sdata, cfg["platform"])
        if cfg["dimension"] == 3:
            el, _ = zstack_element(run, sdata)
            if chunked and not flat_image(run):
                kw["image_key"] = ZKEY                  # the element every tile carries
            else:
                # the lazy element over the source file: its planes are decoded straight from the
                # file, in parallel
                sdata.images[ZSRC] = el                 # in memory only
                kw["image_key"] = ZSRC
        else:
            kw["image_key"] = image_key_2d(run, sdata)
        if chunked:
            adata = an.downsample_vector(z=cfg["dimension"], chunked=True,
                                         sdata_dir=str(run.udir(unit) / "crops"), **kw)
        else:
            adata = an.downsample_vector(sdata, z=cfg["dimension"], **kw)
    if adata.n_obs == 0:
        raise RuntimeError("no anchors were sampled: check revert_image (stain polarity), otsu_mult and "
                           "whether the image overlays the molecules")
    adata.write_h5ad(run.udir(unit) / "anchors.h5ad")
    return {"anchors": int(adata.n_obs), "genes": int(adata.n_vars)}


# ------------------------------------------------------------------------------------------- cluster
def _gpu_oom(e) -> bool:
    return isinstance(e, MemoryError) or "memory" in str(e).lower()


def chunked_pca(adata, n_comps, chunk=500_000, gpu=True):
    """Exact PCA (covariance eigendecomposition) that never holds more than one chunk on the device.

    rsc's chunked=True is cuML IncrementalPCA, whose fit() copies the whole dense matrix to the GPU
    first, so it fails exactly when it is needed. With a few thousand genes the covariance is small:
    accumulate X^T X and the column sums chunk by chunk, eigendecompose, then project each chunk.
    On the GPU an out-of-memory halves the chunk; below 50k rows the same computation runs on the CPU.

    The retry runs after the except block has ended: inside it, the exception's traceback keeps the
    failed attempt's frame, and with it its GPU blocks, alive, so free_all_blocks() could not return
    them and every smaller chunk failed too (a 500k-row block, ~20 GB, held while 250k, 125k and 62.5k
    were tried on a 24 GB card)."""
    import numpy as np
    if gpu:
        import cupy as cp
        while True:
            try:
                return _covariance_pca(adata, n_comps, chunk, cp)
            except (MemoryError, RuntimeError) as e:
                if not _gpu_oom(e):
                    raise
            cp.get_default_memory_pool().free_all_blocks()
            if chunk // 2 < 50_000:
                log("GPU still short: covariance PCA on the CPU")
                break
            log(f"GPU out of memory in chunked PCA at {chunk:,} rows; {chunk // 2:,}")
            chunk //= 2
    return _covariance_pca(adata, n_comps, chunk, np)


def _covariance_pca(adata, n_comps, chunk, xp):
    import numpy as np
    import scipy.sparse as sp
    X = adata.X
    n, g = X.shape
    on_gpu = xp is not np

    def block(i):
        b = X[i:i + chunk]
        b = b.toarray() if sp.issparse(b) else np.asarray(b)
        return xp.asarray(b, dtype=xp.float64)

    s = xp.zeros(g, xp.float64)
    xtx = xp.zeros((g, g), xp.float64)
    for i in range(0, n, chunk):
        b = block(i)
        s += b.sum(0)
        xtx += b.T @ b
        del b
    mu = s / n
    cov = (xtx - n * xp.outer(mu, mu)) / (n - 1)
    host = (lambda x: x.get()) if on_gpu else (lambda x: x)
    w, v = np.linalg.eigh(host(cov))                     # genes x genes: trivial on the CPU
    order = np.argsort(w)[::-1][:n_comps]
    w, v = w[order], v[:, order]
    # deterministic sign: the largest-magnitude loading of each component is positive
    v = v * np.sign(v[np.abs(v).argmax(0), np.arange(v.shape[1])])
    total = float(np.trace(host(cov)))
    v_h, v = v, xp.asarray(v)
    out = np.empty((n, n_comps), np.float32)
    for i in range(0, n, chunk):
        b = block(i)
        r = (b - mu) @ v
        out[i:i + chunk] = (r.get() if on_gpu else r).astype(np.float32)
        del b, r
    adata.obsm["X_pca"] = out
    adata.varm["PCs"] = v_h.astype(np.float32)
    adata.uns["pca"] = {"variance": w.astype(np.float32),
                        "variance_ratio": (w / total).astype(np.float32),
                        "params": {"zero_center": True, "method": "chunked_covariance"}}
    if on_gpu:
        xp.get_default_memory_pool().free_all_blocks()
    return "gpu-chunked" if on_gpu else "cpu-covariance"


def st_cluster(run, unit):
    import anndata as ad
    import scanpy as sc
    C = run.cfg["cluster"]
    adata = ad.read_h5ad(run.udir(unit) / "anchors.h5ad")
    rsc = None
    if device(run) is not None:
        try:
            import cupy as cp
            import rapids_singlecell as rsc  # noqa: F811
            dev = int(str(device(run)).split(":")[-1]) if ":" in str(device(run)) else 0
            cp.cuda.Device(dev).use()
        except Exception as e:  # noqa: BLE001
            log(f"rapids_singlecell unavailable ({e}); clustering on the CPU")
            rsc = None
    n_comps = int(min(C["n_comps"], adata.n_vars - 1, adata.n_obs - 1))
    pca_mode = "cpu"
    if rsc is not None:
        import cupy as cp
        free = cp.cuda.runtime.memGetInfo()[0]
        need = adata.n_obs * adata.n_vars * 4 * 3          # X, its centred copy and workspace (float32)
        short = None
        try:
            if need > 0.6 * free:
                raise MemoryError(f"PCA needs ~{need / 2**30:.0f} GiB, {free / 2**30:.0f} GiB free on the GPU")
            rsc.pp.pca(adata, n_comps=n_comps)
            pca_mode = "gpu"
        except Exception as e:  # noqa: BLE001 -- GPU memory taken by another job
            if not _gpu_oom(e):
                raise
            short = str(e)[:120]
        if short is not None:
            # out of the except block, so the failed attempt's GPU blocks can be freed (chunked_pca)
            cp.get_default_memory_pool().free_all_blocks()
            log(f"full GPU PCA does not fit ({short}); chunked covariance PCA")
            pca_mode = chunked_pca(adata, n_comps)
        # the GPU may be shared with another job: on out-of-memory free the pool and wait for the
        # other job to release memory, then retry; the CPU is the last resort
        graph_done = False
        for attempt in range(10):
            try:
                rsc.pp.neighbors(adata, n_neighbors=int(C["n_neighbors"]))
                rsc.tl.leiden(adata, resolution=float(C["resolution"]))
                graph_done = True
                break
            except (MemoryError, RuntimeError) as e:
                if not _gpu_oom(e):
                    raise
            cp.get_default_memory_pool().free_all_blocks()
            log(f"GPU out of memory in neighbors/leiden (attempt {attempt + 1}/10); waiting 120 s")
            time.sleep(120)
        if not graph_done:
            log("GPU still short after 10 attempts: neighbors/leiden on the CPU")
            sc.pp.neighbors(adata, n_neighbors=int(C["n_neighbors"]), random_state=0)
            sc.tl.leiden(adata, resolution=float(C["resolution"]), random_state=0, flavor="igraph", n_iterations=2)
            rsc = None
    else:
        sc.pp.pca(adata, n_comps=n_comps, svd_solver="randomized", random_state=0)
        sc.pp.neighbors(adata, n_neighbors=int(C["n_neighbors"]), random_state=0)
        sc.tl.leiden(adata, resolution=float(C["resolution"]), random_state=0, flavor="igraph", n_iterations=2)
    out = run.cls_path(unit)
    adata.write_h5ad(out)
    return {"clusters": int(adata.obs["leiden"].nunique()), "anchors": int(adata.n_obs), "h5ad": str(out),
            "pca": pca_mode, "graph": "gpu" if rsc is not None else "cpu"}


# ----------------------------------------------------------------------------------- map and filter
def annotation_outputs(run, unit):
    out = run.udir(unit) / "annotation_out"
    anno = sorted(out.glob("*_annotated.h5ad"))
    markers = out / "celltype_markers.yaml"
    return (anno[0] if anno else None), (markers if markers.exists() else None)


def st_map(run, unit):
    cfg = run.cfg
    anno, markers = annotation_outputs(run, unit)
    if anno is None:
        raise FileNotFoundError(f"no annotated h5ad in {run.udir(unit) / 'annotation_out'}; run the annotation first")
    sdata = read_sdata(run.udir(unit) / "spatialdata.zarr")
    an = analysis(run, unit)
    kw = dict(method=cfg["map"]["method"], h5ad=str(anno), gpu_device=device(run) or "cuda",
              force_cpu=device(run) is None)
    if markers is not None:
        kw["marker_genes_path"] = str(markers)
    if cfg["regime"] == "transcript":
        D = cfg["downsample"]
        kw.update(kernel=D["kernel"], bandwidth=float(D["bandwidth"]), sampling_distance=float(D["sampling_distance"]),
                  element_name=transcripts_key(sdata, cfg["platform"]))
    sdata = an.map_celltypes(sdata, **kw)
    persist(sdata, ["ssaft_downsampled_vectors"])
    names, _ = an.dataset.celltype_names_colors()
    # the map splits the Unknown anchors into one class per cluster (unknown1, unknown2, ...): not cell types
    unknown = sum(1 for n in names or [] if re.fullmatch(r"unknown\d+", n))
    return {"celltypes": len(names or []) - unknown, "unknown_groups": unknown, "annotated_h5ad": str(anno),
            "markers": str(markers) if markers else None}


def st_filter(run, unit):
    F = run.cfg["filter"]
    sdata = read_sdata(run.udir(unit) / "spatialdata.zarr")
    an = analysis(run, unit)
    kw = dict(norm_f=F.get("norm_f"), show_unknown=bool(F.get("show_unknown", True)))
    if F.get("min_p") is not None:
        kw["min_p"] = float(F["min_p"])
    if F.get("hsm_f") is not None:
        kw["hsm_f"] = float(F["hsm_f"])
    sdata = an.filter_celltypemaps(sdata, **kw)
    persist(sdata, ["ssaft_celltype_map", "ssaft_celltype_map_rgb"])
    import numpy as np
    fm = np.asarray(an.dataset.zarr_group["filtered_celltype_maps"][:])
    return {"kept_voxel_fraction": round(float((fm >= 0).mean()), 4)}


# ------------------------------------------------------------------------------------------- segment
def nuclei_element(run, unit, sdata):
    S, plat = run.cfg["segment"], run.cfg["platform"]
    key = S.get("nuclei", "auto")
    if key not in ("auto", "file"):
        return key
    if key == "auto":
        default = {"xenium": "nucleus_boundaries", "stereoseq": "cells_polygons", "cosmx": "cell_centroids"}.get(plat)
        if default and (default in sdata.shapes or default in sdata.points):
            return default
    f = S.get("nuclei_file")
    if not f:
        raise ValueError(f"no nuclei for {plat}: set nuclei_file (polygons as parquet / geojson / "
                         f"csv with a WKT 'geometry' column)")
    import geopandas as gpd
    import pandas as pd
    from shapely import from_wkt
    from spatialdata.models import ShapesModel
    from spatialdata.transformations import get_transformation
    f = Path(f)
    if f.suffix == ".parquet":
        gdf = gpd.read_parquet(f)
    elif f.suffix in (".geojson", ".json", ".gpkg", ".shp"):
        gdf = gpd.read_file(f)
    else:
        df = pd.read_csv(f)
        gdf = gpd.GeoDataFrame(geometry=gpd.GeoSeries(from_wkt(df["geometry"].to_numpy())).to_numpy())
    gdf = gpd.GeoDataFrame(geometry=gdf.geometry.to_numpy())
    off = S.get("nuclei_offset_px")
    if off:
        gdf = gpd.GeoDataFrame(geometry=gdf.geometry.translate(float(off[0]), float(off[1])).to_numpy())
    if S.get("nuclei_frame", "image") == "image":
        ref = sdata.images[image_key_2d(run, sdata)] if run.cfg["regime"] == "transcript" else \
            sdata.images[[k for k in sdata.images if k.endswith("_full_image") or k.endswith("_regist")][0]]
    else:                                       # 'um': the transcripts' own frame
        ref = sdata.points[transcripts_key(sdata, plat)]
    sdata.shapes["ssaft_nuclei"] = ShapesModel.parse(gdf, transformations=dict(get_transformation(ref, get_all=True)))
    persist(sdata, ["ssaft_nuclei"])
    return "ssaft_nuclei"


def nuclei_z_from_transcripts(sdata, nuc, tx_key):
    """Each nucleus's z: the median z of the transcripts inside it (Xenium: overlaps_nucleus, matched
    by cell_id), else of all its cell's transcripts; NaN when neither exists. None when the
    transcripts carry no per-cell assignment (the segmentation then seeds on the tissue midline)."""
    import numpy as np
    if tx_key is None or tx_key not in sdata.points:
        return None
    pts = sdata.points[tx_key]
    if "cell_id" not in pts.columns or "z" not in pts.columns:
        return None
    cols = ["z", "cell_id"] + (["overlaps_nucleus"] if "overlaps_nucleus" in pts.columns else [])
    df = pts[cols].compute()
    df["cell_id"] = df["cell_id"].astype(str)
    ids = sdata.shapes[nuc].index.astype(str)
    z = df.groupby("cell_id", observed=True)["z"].median().reindex(ids).to_numpy(dtype=np.float64)
    if "overlaps_nucleus" in cols:
        inn = df[df["overlaps_nucleus"].astype(int) == 1]
        zn = inn.groupby("cell_id", observed=True)["z"].median().reindex(ids).to_numpy(dtype=np.float64)
        z = np.where(np.isfinite(zn), zn, z)
    log(f"nucleus z from transcripts: {int(np.isfinite(z).sum()):,} of {len(ids):,} nuclei")
    return z if np.isfinite(z).any() else None


def st_segment(run, unit):
    sdata = read_sdata(run.udir(unit) / "spatialdata.zarr")
    an = analysis(run, unit)
    nuc = nuclei_element(run, unit, sdata)
    kw = {}
    if run.cfg["regime"] == "transcript":
        kw["transcripts"] = transcripts_key(sdata, run.cfg["platform"])
    seed_rule = (run.cfg.get("segment") or {}).get("seed_rule")
    if run.cfg["dimension"] == 3 and seed_rule:
        kw["seed_rule"] = seed_rule
    if run.cfg["dimension"] == 3 and nuc in sdata.shapes and seed_rule != "all":   # 'all' does not read it
        z = nuclei_z_from_transcripts(sdata, nuc, kw.get("transcripts"))
        if z is not None:
            kw["nuclei_z"] = z
    seg = run.cfg.get("segment") or {}
    for k in ("min_area_um2", "max_area_um2", "simplify"):   # run.json may override the defaults
        if seg.get(k) is not None:
            kw[k] = seg[k]
    sdata = an.segment_cells(sdata, nuclei=nuc, **kw)
    keys = ["ssaft_cell_boundaries_2d", "ssaft_cells"]
    if run.cfg["dimension"] == 3:
        keys += ["ssaft_cell_boundaries_3d", "ssaft_cell_labels_3d"]
    persist(sdata, keys)
    cells = sdata.tables["ssaft_cells"]
    cells.write_h5ad(run.udir(unit) / "ssaft_cells.h5ad")
    return {"cells": int(cells.n_obs), "nuclei": nuc,
            "top_celltypes": cells.obs["cellman_annotation"].astype(str).value_counts().head(8).to_dict()}


def stages():
    from ssaft_load import st_load
    return dict(load=st_load, zstack=st_zstack, tile=st_tile, panel_vf=st_panel_vf, downsample=st_downsample,
                cluster=st_cluster, map=st_map, filter=st_filter, segment=st_segment)
