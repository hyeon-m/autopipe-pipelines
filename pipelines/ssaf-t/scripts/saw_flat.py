"""Stereo-seq SAW output with every file in one folder (<chip>.tissue.gem.gz, <chip>_regist.tif,
<chip>.cellbin.gef, <chip>.tissue.gef, ...). spatialdata_io.stereoseq reads only the numbered SAW
layout (03.register/, 04.tissuecut/, ...); this folder is read here into
the same SpatialData:

  * tables['table']: the DNBs of <chip>.tissue.gem.gz summed into bins of load.bin_um µm: per axis
    the bin index (x - x.min) // step, step = bin_um / 0.5 µm DNBs; each bin placed at its centre in
    the chip's µm frame, (index * step + x.min + (step - 1) / 2) * 0.5. GridDataset.to_spatialdata
    stores the GridDataset metadata in uns['cellman'];
  * images['staining']: <chip>_regist.tif, registered to the DNB grid (1 px = 1 DNB = 0.5 µm), with
    the pixel -> µm affine diag(0.5, 0.5);
  * shapes['cells_polygons'] the cells of
    <chip>.cellbin.gef, each polygon its cellBorder offsets added to the cell's x, y as
    spatialdata_io.stereoseq builds cells_polygons, here in µm (x 0.5) so that they lie in the frame
    of the bins and the image. The segmentation seeds on them.

The stain polarity (he_revert) is true for an H&E image (HE in the file name) and false for ssDNA
(fluorescence: bright tissue on dark); with several registered images the
ssDNA one is used, as ssaft's own Stereo-seq reading picks it (ssaft._io._spatialdata_io).

The bin size (load.bin_um) is vf_bin_um, by default the GEM's own (BinSize x 0.5 µm: 0.5 µm for the
bin1 GEM SAW writes). A folder holding a GEM table but no SAW
.gef file is not a SAW output and is not read (reader_files.near_miss says so).

plan.py calls inspect() (headers and file names only); ssaft_load.py calls load() at the load stage.
"""
from __future__ import annotations

import gzip
from pathlib import Path

DNB_UM = 0.5                       # Stereo-seq DNB pitch
PADDING = 32767                    # cellBorder padding value (SAW / spatialdata_io)


def _first(items):
    return items[0] if items else None


def is_flat(p: Path) -> bool:
    """SAW output in one folder: its .gef files directly in it, no numbered layout."""
    p = Path(p)
    return not (p / "03.register").is_dir() and bool(list(p.glob("*.tissue.gef")) + list(p.glob("*.cellbin.gef")))


def gem_file(p: Path):
    """The tissue-cut GEM, or else any GEM of the chip (the prep script's order)."""
    p = Path(p)
    return _first(sorted(p.glob("*.tissue.gem.gz"))) or _first(sorted(p.glob("*.gem.gz")))


def images(p: Path):
    """The registered images, ssDNA first (an H&E one has HE in its name)."""
    regist = sorted(Path(p).glob("*regist.tif"))
    return sorted(regist, key=lambda f: is_he(f))


def is_he(f) -> bool:
    return "HE" in Path(f).stem.split("_")


def cellbin_file(p: Path):
    """The cell segmentation: <chip>.cellbin.gef, else the expanded <chip>.adjusted.cellbin.gef."""
    gefs = sorted(Path(p).glob("*.cellbin.gef"))
    return _first([g for g in gefs if ".adjusted." not in g.name]) or _first(gefs)


def gem_header(gem: Path) -> tuple[dict, list[str], int]:
    """The '#key=value' header of a GEM, its column names, and how many '#' lines it has."""
    head, cols, skip = {}, [], 0
    with gzip.open(gem, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                skip += 1
                k, _, v = line[1:].strip().partition("=")
                head[k] = v
                continue
            cols = line.rstrip("\n").split("\t")
            break
    return head, cols, skip


def data_bin_um(head: dict) -> float:
    """The GEM's own bin size in µm: BinSize DNBs x 0.5 µm."""
    try:
        return float(head.get("BinSize") or 1) * DNB_UM
    except ValueError:
        return DNB_UM


def inspect(p: Path) -> dict:
    """The inspect_raw-style report of a flat SAW folder (names and the GEM header only)."""
    p = Path(p)
    gem = gem_file(p)
    head, cols, _ = gem_header(gem) if gem else ({}, [], 0)
    imgs = images(p)
    cb = cellbin_file(p)
    chip = head.get("STOmicsChip") or (gem.name.split(".")[0] if gem else None) or \
        (_first(sorted(p.glob("*.tissue.gef"))).name.split(".")[0] if list(p.glob("*.tissue.gef")) else None)
    return dict(platform="stereoseq", supported=True, dimensions=[2], zstack=False, layout="saw_flat",
                files=dict(gem=str(gem) if gem else None, regist=[f.name for f in imgs],
                           cellbin=str(cb) if cb else None,
                           gef=sorted(x.name for x in p.glob("*.gef"))),
                hints=dict(dataset_id=chip, gem_header=head, gem_columns=cols,
                           ssdna=[f.name for f in imgs if not is_he(f)], he=[f.name for f in imgs if is_he(f)]))


# ------------------------------------------------------------------------------- load (run time)
def cell_polygons(path, log=print):
    """shapes['cells_polygons'] from a cellbin.gef, in µm (spatialdata_io.stereoseq's polygons x 0.5)."""
    import geopandas as gpd
    import h5py
    import numpy as np
    import pandas as pd
    from shapely import Polygon
    from spatialdata.models import ShapesModel
    with h5py.File(str(path), "r") as f:
        cell = f["cellBin"]["cell"][:]
        border = f["cellBin"]["cellBorder"][:]
    polys, ids, skipped = [], [], 0
    for i in range(len(cell)):
        b = border[i]
        b = b[(b[:, 0] != PADDING) & (b[:, 1] != PADDING)].astype(np.float64)
        if len(b) < 3:
            skipped += 1
            continue
        b[:, 0] += float(cell["x"][i])
        b[:, 1] += float(cell["y"][i])
        polys.append(Polygon(b * DNB_UM))
        ids.append(str(int(cell["id"][i])))
    gdf = gpd.GeoDataFrame(geometry=polys, index=pd.Index(ids, name="cell_id"))
    log(f"cells: {len(polys):,} polygons from {Path(path).name}"
        + (f" ({skipped} with fewer than 3 border points left out)" if skipped else ""))
    return ShapesModel.parse(gdf), len(polys)


def load(cfg: dict, log=print):
    """The SpatialData of a flat SAW folder (the bins, the registered image and the cells)."""
    import time

    import anndata as ad
    import numpy as np
    import pandas as pd
    import pyarrow.csv as pacsv
    import tifffile
    from scipy.sparse import coo_matrix
    from ssaft._griddataset import GridDataset

    L = cfg["load"]
    bin_um = float(L["bin_um"])
    gem, img_file = Path(L["gem"]), Path(L["image"])
    t0 = time.time()
    # the '#' header lines are counted here: pyarrow's CSV reader has no comment character
    _, _, skip = gem_header(gem)
    tbl = pacsv.read_csv(gem, read_options=pacsv.ReadOptions(use_threads=True, skip_rows=skip,
                                                             block_size=1 << 26),
                         parse_options=pacsv.ParseOptions(delimiter="\t"))
    names = tbl.column_names
    gene_col = "geneName" if "geneName" in names else "geneID"      # GEMv0.2 carries both
    cnt_col = next(c for c in ("MIDCount", "MIDCounts", "UMICount", "count") if c in names)
    genes = tbl.column(gene_col).dictionary_encode().combine_chunks()
    gi = np.asarray(genes.indices if hasattr(genes, "indices") else genes.chunk(0).indices)
    gnames = [str(v) for v in (genes.dictionary if hasattr(genes, "dictionary")
                               else genes.chunk(0).dictionary)]
    x = tbl.column("x").to_numpy().astype(np.int64)
    y = tbl.column("y").to_numpy().astype(np.int64)
    c = tbl.column(cnt_col).to_numpy().astype(np.float32)
    del tbl
    log(f"{gem.name}: {len(c):,} GEM rows, {len(gnames):,} genes ({gene_col}, {cnt_col}) ({time.time() - t0:.0f}s)")

    step = int(round(bin_um / DNB_UM))                              # DNBs per bin side
    bx, by = (x - x.min()) // step, (y - y.min()) // step
    ny = int(by.max()) + 1
    flat = bx.astype(np.int64) * ny + by
    uniq, inv = np.unique(flat, return_inverse=True)
    X = coo_matrix((c, (inv.astype(np.int64), gi)), shape=(len(uniq), len(gnames))).tocsr()
    X.sum_duplicates()
    ux, uy = uniq // ny, uniq % ny
    raw = ad.AnnData(X=X, obs=pd.DataFrame(index=[f"bin_{a}_{b}" for a, b in zip(ux, uy)]),
                     var=pd.DataFrame(index=gnames))
    raw.var_names_make_unique()
    # bin centre in µm, in the chip frame the registered image is in (origin = DNB 0)
    raw.obsm["spatial"] = np.stack([(ux * step + x.min() + (step - 1) / 2) * DNB_UM,
                                    (uy * step + y.min() + (step - 1) / 2) * DNB_UM], 1)
    n_gem = len(c)
    del x, y, c, gi, flat, uniq, inv

    img = tifffile.imread(img_file)
    if img.ndim == 3 and img.shape[0] in (3, 4):
        img = np.moveaxis(img, 0, -1)
    he = bool(L.get("he_revert"))
    log(f"{raw.n_obs:,} bins of {bin_um:g} µm x {raw.n_vars:,} genes | image {img_file.name} {img.shape} "
        f"at {DNB_UM} µm/px, he_revert {he} ({time.time() - t0:.0f}s)")
    gd = GridDataset(raw=raw, image=img, um_per_px=DNB_UM, species=cfg["species"], organ=cfg["organ"],
                     platform="stereoseq",
                     affine_px_to_um=np.array([[DNB_UM, 0.0, 0.0], [0.0, DNB_UM, 0.0]]),
                     he_revert=he, meta={"bin_um": bin_um})
    sdata = gd.to_spatialdata()
    shape = list(img.shape)
    n_bins, n_genes = int(raw.n_obs), int(raw.n_vars)
    del gd, img, raw
    n_cells = None
    if L.get("cellbin"):
        sdata.shapes["cells_polygons"], n_cells = cell_polygons(L["cellbin"], log)
    return sdata, {"gem_rows": n_gem, "bins": n_bins, "genes": n_genes, "bin_um": bin_um,
                   "image": img_file.name, "image_shape": shape, "he_revert": he,
                   "cells": n_cells, "cellbin": Path(L["cellbin"]).name if L.get("cellbin") else None,
                   "loader": "saw_flat.py (the GEM bins, the registered image, the cellbin.gef cells)"}
