"""The load stage: the input as the run's SpatialData zarr (units/<unit>/spatialdata.zarr).

  raw Xenium        spatialdata_io.xenium, only what the pipeline reads (transcripts,
                    morphology_focus, nucleus boundaries, the cell table)
  raw MERSCOPE      spatialdata_io.merscope without its mosaics; every stain of the mask layer read
                    here, the transcripts with pyarrow (merscope_load); elements named after the data
                    folder (ssaft_images.rename_merscope)
  raw CosMx         the whole slide as one unit from the flat files (load_cosmx_slide): FOV images
                    pasted at their positions, each flipped when the export needs it, transcripts at
                    their global position
  raw Visium HD     spatialdata_io.visium_hd at the bin size chosen
  raw Stereo-seq    the numbered SAW layout: spatialdata_io.stereoseq; a SAW output in one folder:
                    saw_flat.py
  SpatialData zarr  copied (load_zarr), SSAF-T outputs of an earlier run dropped; a per-FOV CosMx zarr
                    (spatialdata_io.cosmx) one unit per FOV (load_cosmx_fovs)
"""
from __future__ import annotations

import re
import shutil
import time
from pathlib import Path

from ssaft_core import log, persist, read_sdata, scale0, transcripts_key

ZKEY = "ssaft_zstack"


def fresh(out: Path):
    """The load is not recorded as done (or it would not run): an existing zarr is the partial
    output of an interrupted attempt of this same run, so it is rebuilt."""
    if out.exists():
        log(f"removing partial output of an interrupted load: {out}")
        shutil.rmtree(out)


def elements(sdata) -> dict:
    return {k: sorted(getattr(sdata, k)) for k in ("images", "labels", "points", "shapes", "tables")}


def st_load(run, unit):
    cfg = run.cfg
    if cfg.get("input_kind") == "zarr":
        return load_zarr(run, unit)
    if (cfg.get("load") or {}).get("layout") == "saw_flat":
        return load_saw_flat(run, unit)
    import spatialdata_io
    p = Path(cfg["dataset_path"])
    plat, L = cfg["platform"], cfg["load"]
    if plat == "cosmx":
        return load_cosmx_slide(run)
    out = run.udir(unit) / "spatialdata.zarr"
    fresh(out)
    if plat == "xenium":
        # the label rasters, the MIP, cell outlines and aligned images would roughly double the zarr copy
        sdata = spatialdata_io.xenium(str(p), cells_boundaries=False, cells_labels=False, nucleus_labels=False,
                                      cells_as_circles=False, morphology_mip=False, aligned_images=False)
    elif plat == "merscope":
        sdata = merscope_load(p, L, z_um=float(L.get("z_plane_um", 1.5)) if cfg["dimension"] == 3 else None)
        if cfg.get("merscope_images"):
            import ssaft_images as IM
            sdata = IM.rename_merscope(sdata, cfg)
    elif plat == "visiumhd":
        kw = {}
        if L.get("fullres_image_file"):
            kw["fullres_image_file"] = L["fullres_image_file"]
        sdata = spatialdata_io.visium_hd(str(p), bin_size=int(L.get("bin_size", 2)), load_all_images=True,
                                         load_segmentations_only=False, **kw)
    elif plat == "stereoseq":
        sdata = spatialdata_io.stereoseq(str(p), dataset_id=L.get("dataset_id"), read_square_bin=True)
    else:
        raise ValueError(plat)
    log("elements", elements(sdata))
    sdata.write(str(out))
    sdata = read_sdata(out)
    add_mask_image(run, sdata)
    return {"zarr": str(out), "images": sorted(sdata.images), "points": sorted(sdata.points),
            "shapes": sorted(sdata.shapes), "tables": sorted(sdata.tables)}


# ----------------------------------------------------------------------------------------- MERSCOPE
def merscope_load(p: Path, L: dict, z_um: float | None = None):
    """spatialdata_io.merscope without its mosaic images, plus the stain(s) in load.stains at
    load.z_layers. The reader otherwise stacks every stain of the plane (each ~21 GB on a full
    section) into one element, which the zarr copy would duplicate.

    The transcripts are read here, not by the reader: its dask read_csv of detected_transcripts.csv
    (tens of GB) with every column, and the full compute that the z conversion needs, took ~35 min and
    >200 GB. A multithreaded pyarrow read of the columns used (x, y, z, gene) takes minutes.
    z_um: µm per z plane (3D); global_z is a plane index. None keeps the plane index."""
    import spatialdata_io
    from spatialdata_io.readers.merscope import _get_channel_names, _get_reader
    z = L.get("z_layers", 3)
    z = z[0] if isinstance(z, list) else int(z)
    sdata = spatialdata_io.merscope(str(p), z_layers=z, mosaic_images=False, cells_boundaries=False,
                                    transcripts=False)
    have = _get_channel_names(p / "images")
    want = [s for s in (L.get("stains") or ["DAPI"]) if s in have]
    if not want:
        raise ValueError(f"none of the stains {L.get('stains') or ['DAPI']} in {p / 'images'}; available: {sorted(have)}")
    dsid = f"{p.parent.name}_{p.name}"                    # the reader's dataset_id
    sdata.points[f"{dsid}_transcripts"] = merscope_points(p, z_um)
    sdata.images[f"{dsid}_z{z}"] = _get_reader(None)(p / "images", want, z,
                                                     {"chunks": (1, 4096, 4096), "scale_factors": [2, 2, 2, 2]})
    log(f"mosaic image: stains {want} at z{z} (of {sorted(have)})")
    return sdata


def merscope_points(p: Path, z_um: float | None):
    import dask.dataframe as dd
    import numpy as np
    import pandas as pd
    import pyarrow as pa
    import pyarrow.csv as pacsv
    from spatialdata.models import PointsModel
    from spatialdata.transformations import Affine
    t0 = time.time()
    tbl = pacsv.read_csv(p / "detected_transcripts.csv",
                         read_options=pacsv.ReadOptions(use_threads=True, block_size=1 << 26),
                         convert_options=pacsv.ConvertOptions(
                             include_columns=["global_x", "global_y", "global_z", "gene"],
                             column_types={"global_x": pa.float64(), "global_y": pa.float64(),
                                           "global_z": pa.float64(), "gene": pa.dictionary(pa.int32(), pa.string())}))
    df = pd.DataFrame({"x": tbl.column("global_x").to_numpy(), "y": tbl.column("global_y").to_numpy(),
                       "z": tbl.column("global_z").to_numpy(),
                       "gene": tbl.column("gene").combine_chunks().to_pandas()})
    del tbl
    if z_um is not None:
        df["z"] *= z_um
    df["gene"] = df["gene"].cat.remove_unused_categories()
    m2p = Affine(np.genfromtxt(p / "images" / "micron_to_mosaic_pixel_transform.csv"),
                 input_axes=("x", "y"), output_axes=("x", "y"))
    ddf = dd.from_pandas(df, npartitions=max(1, len(df) // 5_000_000))
    pts = PointsModel.parse(ddf, coordinates={"x": "x", "y": "y", "z": "z"}, feature_key="gene",
                            transformations={"global": m2p})
    log(f"transcripts: {len(df):,} read with pyarrow in {time.time() - t0:.0f}s"
        + (f"; z plane index x {z_um} -> µm" if z_um is not None else ""))
    return pts


def merscope_z_to_um(sdata, z_plane_um):
    from spatialdata.models import PointsModel
    from spatialdata.transformations import get_transformation
    key = transcripts_key(sdata, "merscope")
    pts = sdata.points[key]
    df = pts.compute().reset_index(drop=True)
    df["z"] = df["z"].astype("float64") * z_plane_um
    fk = pts.attrs.get("spatialdata_attrs", {}).get("feature_key", "gene")
    sdata.points[key] = PointsModel.parse(df, coordinates={"x": "x", "y": "y", "z": "z"}, feature_key=fk,
                                          transformations=dict(get_transformation(pts, get_all=True)))
    log(f"points[{key!r}].z: plane index x {z_plane_um} -> µm")


def add_mask_image(run, sdata):
    """image.channel set: a single-channel copy of that channel becomes the tissue-mask image."""
    ch = run.cfg["image"].get("channel")
    if ch is None:
        return
    from spatialdata.models import Image2DModel
    from spatialdata.transformations import get_transformation
    src_key = run.cfg["image"].get("key")
    if not src_key or src_key == "auto":
        src_key = [k for k in sdata.images if not k.startswith("ssaft_")][0]
    el = sdata.images[src_key]
    node = scale0(el)
    names = [str(c) for c in node.coords["c"].values] if "c" in node.coords else []
    idx = names.index(ch) if isinstance(ch, str) and ch in names else int(ch)
    sdata.images["ssaft_mask_image"] = Image2DModel.parse(
        node.data[idx:idx + 1], dims=("c", "y", "x"), transformations=dict(get_transformation(el, get_all=True)),
        scale_factors=[2, 2, 2])
    persist(sdata, ["ssaft_mask_image"])


# -------------------------------------------------------------------------------------------- CosMx
def cosmx_px_um(path: Path, L: dict) -> float:
    """µm per CosMx image pixel, from the FOV image size unless load.px_um_override is set: the
    5472 x 3648 camera (the public NSCLC, breast, colon ... flat files) is 0.18 µm/px, the newer
    4256 x 4256 one 0.12028 µm/px."""
    if L.get("px_um_override"):
        return float(L["px_um_override"])
    from PIL import Image
    f = next(iter(sorted((path / "CellComposite").glob("*"))), None)
    if f is None:
        return float(L.get("px_um", 0.18))
    w, h = Image.open(f).size
    return 0.18 if (w, h) == (5472, 3648) else 0.12028 if (w, h) == (4256, 4256) else float(L.get("px_um", 0.18))


def y_global_sign(path):
    """+1 if a cell's global y rises with its local y in this export, -1 if it falls. From the
    metadata file's CenterY_local_px / CenterY_global_px; +1 when there is nothing to measure."""
    import numpy as np
    import pandas as pd

    meta = next(path.glob("*metadata_file*"), None)
    if meta is None:
        return 1
    try:
        t = pd.read_csv(meta, usecols=lambda c: c.strip('"') in ("fov", "CenterY_local_px",
                                                                 "CenterY_global_px"))
        t.columns = [c.strip('"') for c in t.columns]
        q = t[t["fov"] == t["fov"].iloc[len(t) // 2]]
        if len(q) < 10:
            return 1
        return 1 if np.corrcoef(q["CenterY_local_px"], q["CenterY_global_px"])[0, 1] > 0 else -1
    except Exception as e:  # noqa: BLE001
        log(f"cosmx: could not read the global-y sign ({e}); assuming it rises with local y")
        return 1


def local_y_is_image_row(path):
    """+1 if a FOV's y_local_px IS the image row, -1 if it is measured from the bottom (row h-1-y).

    CosMx exports disagree, and getting it wrong pastes every FOV image upside down. The vendor's own
    CellLabels mask settles it: a label's row in the mask vs its CenterY_local_px in the metadata
    correlate +1 in one convention and -1 in the other. With no mask to ask, assume the
    bottom-origin convention, which is what the pixel-column exports use.
    """
    import numpy as np
    import pandas as pd
    from scipy import ndimage

    labs = sorted((path / "CellLabels").glob("*.tif*"))
    meta = next(path.glob("*metadata_file*"), None)
    if not labs or meta is None:
        return -1
    try:
        import tifffile
        f = labs[len(labs) // 2]
        fov = int("".join(c for c in f.stem.split("_F")[-1] if c.isdigit()))
        m = tifffile.imread(f)
        t = pd.read_csv(meta, usecols=lambda c: c.strip('"') in ("fov", "cell_ID", "CenterY_local_px"))
        t.columns = [c.strip('"') for c in t.columns]
        t = t[t["fov"] == fov]
        ids = np.unique(m)
        ids = ids[ids > 0]
        rows = np.asarray(ndimage.center_of_mass(m > 0, m, ids))[:, 0]
        q = pd.DataFrame({"cell_ID": ids, "row": rows}).merge(t, on="cell_ID")
        if len(q) < 10:
            return -1
        return 1 if np.corrcoef(q["row"], q["CenterY_local_px"])[0, 1] > 0 else -1
    except Exception as e:  # noqa: BLE001
        log(f"cosmx: could not read the CellLabels convention ({e}); assuming a bottom-origin y_local")
        return -1


def fov_corners_from_transcripts(tx, h, path):
    """Global-pixel top corner of every FOV plus whether its image must be flipped, read off the
    transcripts (fov, x/y_local_px, x/y_global_px).

    Used when the vendor's fov_positions file gives stage millimetres rather than pixels. Two signs
    decide the geometry: s = +1 if y_global rises with y_local (median-absolute-deviation picks it),
    and k = +1 if y_local is the image row (local_y_is_image_row). The image is pasted flipped when
    the two disagree, and the tile's top is median(y_global - s*y_local), less h-1 when s = -1.
    """
    import pandas as pd
    import pyarrow.csv as pacsv

    t = pacsv.read_csv(tx, convert_options=pacsv.ConvertOptions(
        include_columns=["fov", "x_local_px", "y_local_px", "x_global_px", "y_global_px"]))
    d = pd.DataFrame({c: t.column(c).to_numpy().astype("float64") for c in
                      ("fov", "x_local_px", "y_local_px", "x_global_px", "y_global_px")})
    del t
    d["dx"] = d["x_global_px"] - d["x_local_px"]
    d["dy_same"] = d["y_global_px"] - d["y_local_px"]
    d["dy_inv"] = d["y_global_px"] + d["y_local_px"]
    g = d.groupby("fov")
    mad = lambda col: float((g[col].transform("median") - d[col]).abs().median())  # noqa: E731
    inverted = mad("dy_inv") < mad("dy_same")
    out = pd.DataFrame({"fov": g["dx"].median().index.astype(int),
                        "x_global_px": g["dx"].median().to_numpy(),
                        "y_global_px": (g["dy_inv"].median().to_numpy() - (h - 1)) if inverted
                        else g["dy_same"].median().to_numpy()})
    k = local_y_is_image_row(Path(tx).parent)
    flip = k != (-1 if inverted else 1)
    log(f"cosmx slide: y_global {'falls' if inverted else 'rises'} with y_local, y_local is "
        f"{'the image row' if k > 0 else 'bottom-origin'} -> image {'flipped' if flip else 'as-is'}")
    return out.reset_index(drop=True), flip


def load_cosmx_slide(run):
    """The whole slide as ONE unit: every FOV's transcripts at their slide (global) position, the FOV
    composite images pasted into one mosaic at their fov_positions, the gaps between FOVs empty (no
    transcripts, a black image -- the KDE is exactly zero there and the map skips it). One annotation
    per dataset instead of one per FOV, and a map that is already the whole slide.

    Reads the flat files directly (tx_file, fov_positions_file, metadata_file, CellComposite/), not
    spatialdata_io.cosmx: that reader also parses exprMat_file, which is malformed in some datasets
    (an object column, a duplicated cell_ID) and is not needed here.
    Global pixel = FOV-local pixel + the FOV's fov_positions offset; whether each FOV image is
    flipped vertically before it is pasted is measured per slide."""
    import numpy as np
    import pandas as pd
    import pyarrow.csv as pacsv
    import spatialdata as sd
    from PIL import Image
    from spatialdata.models import Image2DModel, PointsModel, ShapesModel
    from spatialdata.transformations import Affine, Scale

    cfg, L = run.cfg, run.cfg["load"]
    path = Path(cfg["dataset_path"])
    px_um = cosmx_px_um(path, L)
    # full resolution: downsample_vector sizes its sampling blocks in the transcripts' pixel units and
    # counts mask pixels of the image, so the image pixel must be the global pixel (as on Xenium)
    f = int(L.get("mosaic_downsample", 1))
    tx = next(path.glob("*tx_file*"))
    t0 = time.time()
    tbl = pacsv.read_csv(tx, convert_options=pacsv.ConvertOptions(
        include_columns=["fov", "x_global_px", "y_global_px", "z", "target"]))
    df = pd.DataFrame({"x": tbl.column("x_global_px").to_numpy() * px_um,
                       "y": tbl.column("y_global_px").to_numpy() * px_um,
                       "target": tbl.column("target").to_pandas().astype("category")})
    if cfg["dimension"] == 3:
        df["z"] = tbl.column("z").to_numpy().astype("float64") * float(L.get("z_plane_um", 0.8))
    del tbl
    to_px = {"global": Scale([1 / px_um, 1 / px_um], axes=("x", "y"))}
    coords = {"x": "x", "y": "y"} | ({"z": "z"} if "z" in df else {})
    log(f"cosmx slide: {len(df):,} transcripts read in {time.time() - t0:.0f}s, {px_um} µm/px")

    imgs = {int("".join(ch for ch in p.stem.split("_F")[-1] if ch.isdigit())): p
            for p in (path / "CellComposite").glob("*")}
    w, h = Image.open(next(iter(imgs.values()))).size
    pos = pd.read_csv(next(path.glob("*fov_positions*")))
    pos.columns = [c.strip().strip('"') for c in pos.columns]
    pos = pos.rename(columns={"FOV": "fov"})
    # whether the FOV image is upside down relative to the global frame is a property of the export,
    # not of CosMx: measure it (COS047-054 flip, COS066/072 do not)
    flip = local_y_is_image_row(path) != y_global_sign(path)
    log(f"cosmx slide: FOV image {'flipped' if flip else 'as-is'} (from CellLabels / metadata)")
    if "x_global_px" not in pos.columns:
        # older exports give the stage position in mm (X_mm / Y_mm) instead of the FOV's global pixel
        # corner. The transcripts carry BOTH coordinate systems, so take the corner from them.
        pos, flip = fov_corners_from_transcripts(tx, h, path)
        log("cosmx slide: fov_positions has no pixel columns; corners taken from the transcripts")
    x0, y0 = float(pos["x_global_px"].min()), float(pos["y_global_px"].min())
    W = int(np.ceil((pos["x_global_px"].max() - x0 + w) / f)) + 1
    H = int(np.ceil((pos["y_global_px"].max() - y0 + h) / f)) + 1
    mosaic = np.zeros((3, H, W), dtype=np.uint8)
    for _, r in pos.iterrows():
        p = imgs.get(int(r["fov"]))
        if p is None:
            continue
        # most exports store the FOV image upside down relative to the transcripts' local y (the
        # transcript density correlates with the image only after a vertical flip: r 0.18-0.49 vs ~0
        # as-is); the mm-position exports do not, so `flip` is decided per slide above
        a = np.asarray(Image.open(p).convert("RGB"))
        a = (a[::-1] if flip else a)[::f, ::f]
        oy, ox = int(round((r["y_global_px"] - y0) / f)), int(round((r["x_global_px"] - x0) / f))
        hh, ww = min(a.shape[0], H - oy), min(a.shape[1], W - ox)
        mosaic[:, oy:oy + hh, ox:ox + ww] = np.moveaxis(a[:hh, :ww], -1, 0)
    img_tf = {"global": Affine(np.array([[f, 0, x0], [0, f, y0], [0, 0, 1]], dtype=float),
                               input_axes=("x", "y"), output_axes=("x", "y"))}
    # chunk the mosaic before it is parsed: a 10-gigapixel slide written as one chunk copies the whole
    # array for every pyramid level and peaks at 376 GiB. Tiles keep the peak near one copy.
    import dask.array as da
    mosaic_da = da.from_array(mosaic, chunks=(3, 4096, 4096))
    log(f"cosmx slide: {len(pos)} FOVs pasted into a {W} x {H} mosaic (1/{f} resolution)")

    els = dict(points={"transcripts": PointsModel.parse(df, coordinates=coords, feature_key="target",
                                                         transformations=to_px)},
               images={"composite": Image2DModel.parse(mosaic_da, dims=("c", "y", "x"),
                                                       transformations=img_tf, scale_factors=[2, 2, 2])})
    meta = next(path.glob("*metadata_file*"), None)
    if meta is not None:
        m = pd.read_csv(meta, usecols=lambda c: c.strip('"') in ("fov", "cell_ID", "CenterX_global_px",
                                                                  "CenterY_global_px"))
        m.columns = [c.strip('"') for c in m.columns]
        m = m.loc[:, ~m.columns.duplicated()]
        xy = np.c_[m["CenterX_global_px"].to_numpy() * px_um, m["CenterY_global_px"].to_numpy() * px_um]
        idx = pd.Index([f"{a}_{b}" for a, b in zip(m["fov"], m["cell_ID"])], name="cell_label")
        els["shapes"] = {"cell_centroids": ShapesModel.parse(xy, geometry=0, radius=np.ones(len(xy)),
                                                             index=idx, transformations=to_px)}
    out = run.udir("sample") / "spatialdata.zarr"
    fresh(out)
    sd.SpatialData(**els).write(str(out))
    run.set_units(["sample"])
    add_mask_image(run, read_sdata(out))
    return {"zarr": str(out), "transcripts": int(len(df)), "fovs": int(len(pos)), "px_um": px_um,
            "mosaic": [W, H], "cells": int(len(els.get("shapes", {}).get("cell_centroids", [])))}


def load_cosmx_fovs(run, raw):
    """A per-FOV CosMx SpatialData (as spatialdata_io.cosmx writes it: <fov>_points, <fov>_image,
    <fov>_labels): one unit per FOV, transcripts in µm (Scale back to pixels in 'global'), the
    composite image in pixels, and the vendor cell-label centroids as segmentation seeds."""
    import numpy as np
    import pandas as pd
    import spatialdata as sd
    from scipy import ndimage
    from spatialdata.models import Image2DModel, PointsModel, ShapesModel
    from spatialdata.transformations import Identity, Scale

    cfg, L = run.cfg, run.cfg["load"]
    px_um = float(L["px_um"])
    fovs = sorted({k.split("_")[0] for k in raw.points})
    if L.get("fovs") not in (None, "all"):
        want = {str(int(f)) for f in L["fovs"]}
        fovs = [f for f in fovs if f in want]
    units, done = [], []
    to_px = {"global": Scale([1 / px_um, 1 / px_um], axes=("x", "y"))}
    for fov in fovs:
        unit = f"fov{fov}"
        out = run.udir(unit) / "spatialdata.zarr"
        units.append(unit)
        fresh(out)
        df = raw.points[f"{fov}_points"].compute().reset_index(drop=True)
        df["x"] = df["x"].astype("float64") * px_um
        df["y"] = df["y"].astype("float64") * px_um
        coords = {"x": "x", "y": "y"}
        if cfg["dimension"] == 3:
            df["z"] = df["z_raw"].astype("float64") * float(L["z_plane_um"])
            coords["z"] = "z"
        fk = "target" if "target" in df.columns else raw.points[f"{fov}_points"].attrs["spatialdata_attrs"]["feature_key"]
        img = raw.images[f"{fov}_image"]
        node = scale0(img)
        els = dict(
            points={"transcripts": PointsModel.parse(df, coordinates=coords, feature_key=fk, transformations=to_px)},
            images={"composite": Image2DModel.parse(node.data, dims=("c", "y", "x"),
                                                    transformations={"global": Identity()})})
        if f"{fov}_labels" in raw.labels:
            lab = np.asarray(scale0(raw.labels[f"{fov}_labels"]).data)
            ids = np.unique(lab)
            ids = ids[ids > 0]
            if len(ids):
                cen = np.asarray(ndimage.center_of_mass(lab > 0, lab, ids))          # (row, col) px
                # shapes, not points: a second points element would compete with the transcripts
                els["shapes"] = {"cell_centroids": ShapesModel.parse(
                    np.c_[cen[:, 1] * px_um, cen[:, 0] * px_um], geometry=0, radius=np.ones(len(ids)),
                    index=pd.Index([str(i) for i in ids], name="cell_label"), transformations=to_px)}
        sd.SpatialData(**els).write(str(out))
        done.append(unit)
        log(f"{unit}: {len(df):,} transcripts -> {out}")
    run.set_units(units)
    for unit in units:
        add_mask_image(run, read_sdata(run.udir(unit) / "spatialdata.zarr"))
    return {"units": units, "written": done, "px_um": px_um}


# ---------------------------------------------------------------------------------- SpatialData zarr
def load_zarr(run, unit):
    import spatialdata as sd
    from ssaft import Analysis
    from ssaft_images import keep_planes
    cfg, Z = run.cfg, run.cfg["zarr"]
    src = cfg["dataset_path"]
    if cfg["platform"] == "cosmx" and Z.get("per_fov_cosmx"):
        log(f"per-FOV CosMx zarr: one unit per FOV from {src}")
        return load_cosmx_fovs(run, sd.read_zarr(src))
    if cfg["platform"] == "cosmx":
        unit = "sample"
        run.set_units(["sample"])
    out = run.udir(unit) / "spatialdata.zarr"
    fresh(out)
    sdata = sd.read_zarr(src)
    stale = set(Analysis._OWN_ELEMENTS) | {"ssaft_mask_image", "ssaft_nuclei", ZKEY}
    for kind in ("images", "labels", "points", "shapes", "tables"):
        for k in [k for k in getattr(sdata, kind) if k in stale]:
            log(f"dropping {kind}[{k!r}]: an SSAF-T output of an earlier run, not input")
            del getattr(sdata, kind)[k]
    if cfg["platform"] == "visiumhd" and not Z.get("cellman_table"):
        # the bin size chosen (load.bin_size, visiumhd_bin_um): the raw reader loads only that one, and
        # the count backend takes the smallest bin table it finds
        keep = f"square_{int(cfg['load'].get('bin_size', 2)):03d}um"
        for kind in ("tables", "shapes"):
            for k in [k for k in getattr(sdata, kind) if re.search(r"square_\d{3}um", k) and keep not in k]:
                log(f"dropping {kind}[{k!r}]: Visium HD bins other than {keep} (visiumhd_bin_um)")
                del getattr(sdata, kind)[k]
    if cfg["regime"] == "count" and not Z.get("cellman_table") and not sdata.attrs.get("spatialdata_io_reader"):
        sdata.attrs["spatialdata_io_reader"] = {"visiumhd": "visium_hd", "stereoseq": "stereoseq"}[cfg["platform"]]
    if cfg["dimension"] == 3 and Z.get("zstack_element"):
        el = sdata.images.pop(Z["zstack_element"])
        planes = (cfg.get("zstack") or {}).get("planes")
        if planes:
            el = keep_planes(el, [int(v) for v in planes])
            log(f"z-stack: planes {planes} of images[{Z['zstack_element']!r}] (z_planes)")
        sdata.images[ZKEY] = el
        log(f"z-stack: images[{Z['zstack_element']!r}] -> images[{ZKEY!r}]")
    if cfg["platform"] == "merscope" and cfg["dimension"] == 3:
        key = transcripts_key(sdata, "merscope")
        zz = sdata.points[key]["z"].drop_duplicates().compute()
        if len(zz) <= 64 and (zz.round() == zz).all():
            merscope_z_to_um(sdata, float(cfg["load"].get("z_plane_um", 1.0)))
    log("elements", elements(sdata))
    sdata.write(str(out))
    sdata = read_sdata(out)
    add_mask_image(run, sdata)
    return {"zarr": str(out), "source": src, "images": sorted(sdata.images), "points": sorted(sdata.points),
            "shapes": sorted(sdata.shapes), "tables": sorted(sdata.tables)}


# ------------------------------------------------------------------------- Stereo-seq SAW in one folder
def load_saw_flat(run, unit):
    """A Stereo-seq SAW output in one folder: saw_flat.py builds the SpatialData from
    the GEM bins, the registered image and the cellbin.gef cells (cells_polygons)."""
    import saw_flat
    out = run.udir(unit) / "spatialdata.zarr"
    fresh(out)
    sdata, summ = saw_flat.load(run.cfg, log)
    log("elements", elements(sdata))
    sdata.write(str(out))
    return dict(summ, zarr=str(out), source=run.cfg["dataset_path"], images=sorted(sdata.images),
                shapes=sorted(sdata.shapes), tables=sorted(sdata.tables))
