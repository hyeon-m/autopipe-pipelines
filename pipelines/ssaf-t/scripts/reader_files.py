"""The files the load stage opens in a raw output folder, checked at dry_run so that a run does not fail
at load because one is missing.

The lists follow the readers as installed in the image (spatialdata_io 0.7.0) called the way the
load stage calls them (ssaft_load.st_load), plus what the load reads itself:

  xenium     spatialdata_io.xenium(cells_boundaries=False, cells_labels=False, nucleus_labels=False,
             cells_as_circles=False, morphology_mip=False, aligned_images=False): experiment.xenium,
             cells.zarr.zip (opened in every call), cells.parquet + cell_feature_matrix.h5 (the cell
             table), nucleus_boundaries.parquet (nucleus_boundaries stays on), transcripts.parquet,
             the morphology focus image (morphology_focus.ome.tif before Xenium Analyzer 2.0, else
             morphology_focus/morphology_focus_0000.ome.tif or morphology_focus/ch0000_dapi.ome.tif);
             downsampling on the z-stack: the z-stack file (morphology.ome.tif)
  merscope   spatialdata_io.merscope(mosaic_images=False, transcripts=False, cells_boundaries=False):
             images/micron_to_mosaic_pixel_transform.csv; merscope_load: detected_transcripts.csv and
             images/mosaic_<stain>_z<layer>.tif of every stain loaded at the mask layer; downsampling
             on the z-stack: the merscope_stain mosaic of every plane chosen. The cell table
             (cell_by_gene.csv, cell_metadata.csv) is optional: the reader skips it when absent
  cosmx      the whole-slide loader (ssaft_load.load_cosmx_slide; it reads the flat files, not
             spatialdata_io.cosmx): *tx_file.csv, *fov_positions_file.csv, CellComposite/ images;
             with segmentation: *metadata_file.csv (the cell centroids)
  visiumhd   spatialdata_io.visium_hd(bin_size=<chosen>): [<id>_]feature_slice.h5,
             binned_outputs/square_<NNN>um/filtered_feature_bc_matrix.h5, .../spatial/
             scalefactors_json.json, .../spatial/tissue_positions.parquet
  stereoseq  the numbered SAW layout, spatialdata_io.stereoseq(read_square_bin=True):
             03.register/*_regist.tif (the image), 04.tissuecut/*.gef (the square bins),
             041.cellcut/*.cellbin.gef, 051.cellcluster/*cell.cluster.h5ad.
             A SAW output in one folder, the pipeline's saw_flat.py:
             <chip>.tissue.gem.gz (the counts), <chip>_regist.tif (the image); with segmentation
             <chip>.cellbin.gef (the cells)

A folder that shows a platform's files but lacks one the inspection (inspect_raw.py) needs to recognise it
(e.g. MERSCOPE mosaics without detected_transcripts.csv, a Stereo-seq GEM table without the SAW .gef
files) is reported by near_miss with the files missing, instead of as an unknown layout.

cosmx_camera reads the CosMx pixel size the way the load does (ssaft_load.cosmx_px_um): from
the size of the first CellComposite image.

Used by plan.py only (no stage imports it).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

READER = {"xenium": "spatialdata_io.xenium", "merscope": "spatialdata_io.merscope and the MERSCOPE loader",
          "cosmx": "the CosMx slide loader", "visiumhd": "spatialdata_io.visium_hd",
          "stereoseq": "spatialdata_io.stereoseq",
          "stereoseq_flat": "the flat SAW loader (saw_flat.py)"}

# the dry_run `loader :` line: what reads a raw folder into the run's SpatialData zarr
LOADER = {
    "xenium": "spatialdata_io.xenium (transcripts, morphology_focus, nucleus boundaries, cell table)",
    "merscope": "spatialdata_io.merscope (the cell table, if any) and the MERSCOPE loader: the z_layers mosaic "
                "of every stain, the transcripts read with pyarrow",
    "cosmx": "the CosMx slide loader (it reads *tx_file.csv, *fov_positions_file.csv, CellComposite/ "
             "and *metadata_file.csv directly and joins the FOVs into one slide; not spatialdata_io.cosmx)",
    "visiumhd": "spatialdata_io.visium_hd (the bin size chosen, all images)",
    "stereoseq": "spatialdata_io.stereoseq (square bins, cell bins, the registered image)",
    "stereoseq_flat": "the pipeline's flat SAW loader (saw_flat.py): the DNBs of the tissue GEM summed into bins of "
                      "vf_bin_um, the registered image (1 px = 1 DNB = 0.5 µm), and the cells of the "
                      "cellbin.gef as cells_polygons; spatialdata_io.stereoseq reads the numbered layout only",
}

# CosMx FOV image size (width, height) -> µm per pixel, as the load decides it
COSMX_CAMERAS = {(5472, 3648): (0.18, "the earlier CosMx camera, 5472 x 3648 px (e.g. the public NSCLC, "
                                      "breast and colon flat files)"),
                 (4256, 4256): (0.12028, "the newer CosMx camera, 4256 x 4256 px")}


def _any(p: Path, pattern: str) -> bool:
    return p.is_dir() and any(p.glob(pattern))


def _xenium_version(p: Path):
    """(major, minor) of the Xenium analysis software in experiment.xenium, or None."""
    try:
        v = json.loads((p / "experiment.xenium").read_text()).get("analysis_sw_version", "")
    except (OSError, ValueError):
        return None
    m = re.search(r"(\d+)\.(\d+)", str(v))
    return (int(m.group(1)), int(m.group(2))) if m else None


def cosmx_camera(p: Path):
    """(width, height, µm per px or None, which camera) of the first CellComposite image, the one the
    load sizes the pixel from; (None, None, None, why) when there is no image to read."""
    d = Path(p) / "CellComposite"
    f = next(iter(sorted(d.glob("*"))), None) if d.is_dir() else None
    if f is None:
        return None, None, None, "no CellComposite image"
    try:
        from PIL import Image
        Image.MAX_IMAGE_PIXELS = None
        with Image.open(f) as im:
            w, h = im.size
    except Exception as e:  # noqa: BLE001
        return None, None, None, f"{f.name} could not be read ({type(e).__name__})"
    px, cam = COSMX_CAMERAS.get((w, h), (None, f"{w} x {h} px, not a known CosMx camera"))
    return w, h, px, cam


def required(plat: str, p: Path, cfg: dict, planes=None) -> list[tuple[str, bool, str]]:
    """[(file, present, what it is for)] for this run of this raw folder. `planes`: the z-stack planes
    chosen, when the run downsamples on a z-stack."""
    p = Path(p)
    stack = (cfg.get("downsample_image") or {}).get("kind") == "zstack"
    seg = "segment" in cfg.get("stages", [])
    rows = []

    def add(name, ok, why):
        rows.append((name, bool(ok), why))

    if plat == "xenium":
        for name, why in (("experiment.xenium", "the run's specifications"),
                          ("cells.zarr.zip", "the cell index, opened by the reader in every run"),
                          ("cells.parquet", "the cell table (metadata)"),
                          ("cell_feature_matrix.h5", "the cell table (counts)"),
                          ("nucleus_boundaries.parquet", "nucleus polygons, read by the reader in every run "
                                                         "(and the segmentation seeds)"),
                          ("transcripts.parquet", "the transcripts")):
            add(name, (p / name).exists(), why)
        ver = _xenium_version(p)
        if ver is not None and ver < (2, 0):
            add("morphology_focus.ome.tif", (p / "morphology_focus.ome.tif").exists(),
                "the 2D image (Xenium Analyzer < 2.0)")
        else:
            d = p / "morphology_focus"
            ok = (d / "morphology_focus_0000.ome.tif").exists() or (d / "ch0000_dapi.ome.tif").exists()
            add("morphology_focus/morphology_focus_0000.ome.tif (or ch0000_dapi.ome.tif)", ok,
                "the 2D image (Xenium Analyzer >= 2.0)")
        if stack:
            src = Path(str(cfg.get("zstack", {}).get("source") or p / "morphology.ome.tif"))
            add(str(src.relative_to(p)) if str(src).startswith(str(p)) else str(src), src.exists(),
                "the z-stack downsampling uses")
    elif plat == "merscope":
        add("detected_transcripts.csv", (p / "detected_transcripts.csv").exists(),
            "the transcripts (read by name by the MERSCOPE loader; a .csv.gz is not read)")
        add("images/micron_to_mosaic_pixel_transform.csv",
            (p / "images" / "micron_to_mosaic_pixel_transform.csv").exists(), "µm to mosaic pixel transform")
        L = cfg.get("load", {})
        z = L.get("z_layers")
        z = z[0] if isinstance(z, list) else z
        for s in L.get("stains") or []:
            add(f"images/mosaic_{s}_z{z}.tif", (p / "images" / f"mosaic_{s}_z{z}.tif").exists(),
                "the 2D mosaic")
        if stack and planes:
            st = cfg.get("zstack", {}).get("stain", "DAPI")
            for k in planes:
                add(f"images/mosaic_{st}_z{k}.tif", (p / "images" / f"mosaic_{st}_z{k}.tif").exists(),
                    "a z-stack plane chosen in z_planes")
    elif plat == "cosmx":
        add("*tx_file.csv", _any(p, "*tx_file*"), "the transcripts")
        add("*fov_positions_file.csv", _any(p, "*fov_positions*"), "the FOV positions")
        add("CellComposite/ (images)", _any(p / "CellComposite", "*"), "the FOV images")
        if seg:
            add("*metadata_file.csv", _any(p, "*metadata_file*"), "the cell centroids (segmentation seeds)")
    elif plat == "visiumhd":
        ok = (p / "feature_slice.h5").exists() or _any(p, "*_feature_slice.h5")
        add("feature_slice.h5 (or <id>_feature_slice.h5)", ok, "the slide layout")
        b = int(cfg.get("load", {}).get("bin_size", 2))
        d = f"binned_outputs/square_{b:03d}um"
        for name, why in ((f"{d}/filtered_feature_bc_matrix.h5", "the bin counts"),
                          (f"{d}/spatial/scalefactors_json.json", "the bin scale factors"),
                          (f"{d}/spatial/tissue_positions.parquet", "the bin positions")):
            add(name, (p / name).exists(), why)
    elif plat == "stereoseq" and cfg.get("load", {}).get("layout") == "saw_flat":
        L = cfg["load"]
        chip = L.get("dataset_id") or "<chip>"

        def there(key, fallback):
            v = L.get(key)
            return (Path(v).name if v else fallback), bool(v) and Path(v).exists()
        name, ok = there("gem", f"{chip}.tissue.gem.gz")
        add(name, ok, "the counts, one row per gene per DNB, summed into bins (SAW writes it; `saw convert "
                      "gef2gem` makes it from the .tissue.gef)")
        name, ok = there("image", f"{chip}_regist.tif")
        add(name, ok, "the image, registered to the DNB grid (1 px = 1 DNB)")
        if seg:
            name, ok = there("cellbin", f"{chip}.cellbin.gef")
            add(name, ok, "the cells (cellBin/cell and cellBorder): the segmentation seeds; or segment: false")
    elif plat == "stereoseq":
        gef = re.compile(r".*(\.raw\.gef|\.cellbin\.gef|\.tissuecut\.gef|\.gef)$")

        def has(sub, rx):
            d = p / sub
            return d.is_dir() and any(rx.match(f.name) for f in d.iterdir())
        add("03.register/*_regist.tif", has("03.register", re.compile(r".*_regist\.tif$")), "the registered image")
        add("04.tissuecut/*.gef", has("04.tissuecut", gef), "the square bins")
        add("041.cellcut/*.cellbin.gef", has("041.cellcut", gef), "the cell bins")
        add("051.cellcluster/*cell.cluster.h5ad", has("051.cellcluster", re.compile(r".*cell\.cluster\.h5ad$")),
            "the cell table")
    return rows


def missing_message(plat: str, p: Path, rows, reader: str | None = None) -> str | None:
    gone = [(n, w) for n, ok, w in rows if not ok]
    if not gone:
        return None
    return (f"{READER.get(reader or plat, plat)} cannot load {p}: these files it opens are missing:\n"
            + "\n".join(f"  - {n} ({w})" for n, w in gone)
            + "\nAdd them to the folder (or, where the message offers one, change the setting that needs "
              "them), then run dry_run again.")


# a platform's traces in a folder, and the files that recognising and loading it then needs
_TRACES = {
    "xenium": (("experiment.xenium", "transcripts.parquet", "morphology.ome.tif", "morphology_focus.ome.tif",
                "morphology_focus", "cells.zarr.zip"),
               ("experiment.xenium", "cells.zarr.zip", "cells.parquet", "cell_feature_matrix.h5",
                "nucleus_boundaries.parquet", "transcripts.parquet")),
    "merscope": (("detected_transcripts.csv", "images/mosaic_*_z*.tif", "cell_by_gene.csv"),
                 ("detected_transcripts.csv", "images/micron_to_mosaic_pixel_transform.csv",
                  "images/mosaic_*_z*.tif")),
    "cosmx": (("*tx_file*", "*fov_positions*", "*exprMat_file*", "CellComposite"),
              ("*tx_file.csv", "*fov_positions_file.csv", "CellComposite/")),
    "visiumhd": (("binned_outputs", "*feature_slice.h5"),
                 ("feature_slice.h5", "binned_outputs/square_002um/filtered_feature_bc_matrix.h5")),
    # a GEM table (and maybe the image) without the SAW .gef files: not a SAW output
    "stereoseq": (("*.gem.gz", "*.gem", "*regist.tif"),
                  ("*.tissue.gef (the SAW output's square-bin GEF, by which a SAW output is recognised)",
                   "*.tissue.gem.gz", "*_regist.tif")),
}


def _present(p: Path, pattern: str) -> bool:
    pattern = pattern.rstrip("/")
    if any(ch in pattern for ch in "*?["):
        return any(p.glob(pattern))
    return (p / pattern).exists()


def near_miss(p: Path) -> str | None:
    """For a folder no layout matched: the platform whose files it shows, and what is missing."""
    p = Path(p)
    best = None
    for plat, (traces, needed) in _TRACES.items():
        found = [t for t in traces if _present(p, t)]
        if found and (best is None or len(found) > len(best[1])):
            best = (plat, found, [n for n in needed if not _present(p, n.replace("(", "").split(" ")[0])])
    if best is None or not best[2]:
        return None
    plat, found, gone = best
    if plat == "stereoseq":
        return (f"{p} holds Stereo-seq files ({', '.join(found)}) but is not a SAW output: a GEM table alone is "
                f"not read. Missing:\n" + "\n".join(f"  - {n}" for n in gone)
                + "\nGive the SAW output folder (in one folder: <chip>.tissue.gef, <chip>.tissue.gem.gz, "
                  "<chip>_regist.tif, <chip>.cellbin.gef; or the numbered layout 03.register/, 04.tissuecut/, "
                  "...), or a SpatialData zarr.")
    return (f"{p} looks like a {plat} output folder (it has {', '.join(found)}), but these files "
            f"{READER[plat]} needs are missing:\n" + "\n".join(f"  - {n}" for n in gone)
            + "\nAdd them (a folder still being copied may not be complete yet), then run dry_run again.")
