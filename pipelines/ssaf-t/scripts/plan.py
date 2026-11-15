#!/usr/bin/env python
"""Plan one SSAF-T run: inspect the input, turn config.yaml into answers, and build run.json with
run_config.build / validate (the platform defaults in platform_defaults.json, the answers over them).

Called by the Snakefile at parse time (so dry_run shows it), with the SSAF-T environment's Python:

    plan.py --config-json '<config as JSON>' --run-dir ssaft

stdout: one JSON object {cfg, inspect, answers, summary, zplanes, images}
stderr: the human-readable summary (platform, loader, files, elements, images, choices, z planes,
        stage plan)
exit 2: the configuration cannot run; stdout holds the message, with the valid choices taken from
        the input wherever there are some. A value that has to come from the user (species or organ,
        the z plane spacing of a 3D map on MERSCOPE or CosMx, the z planes of a z-stack that lacks
        some) does not stop the plan where it is found: the whole plan made from the data is printed
        first (platform, files, the images and the z plane table by plane number, ...), then every
        value still needed.

Two choices are made apart (SSAF-T 23dd889): the map's dimension (`dimension`, 2D or 3D) and the
image downsampling uses (`image_key`: a 2D image, or a z-stack). A 2D map uses a 2D image; a 3D map
uses either, a 2D image being placed by SSAF-T at the middle of the transcripts' z range (no z-stack
step, no planes to choose). The candidates, the recommendation and the choice come from
ssaft_images.py.

Raw dataset folders are inspected with inspect_raw.py and loaded by ssaft_load.py with the
platform's spatialdata_io reader (CosMx: from its flat files); the files the load opens are checked
here (reader_files.py). A Stereo-seq SAW output with every file in one folder, which
spatialdata_io.stereoseq does not read, is inspected and loaded by saw_flat.py. SpatialData .zarr inputs are inspected here from the zarr metadata
and the element names, without loading any data. A run downsampling on a Xenium / MERSCOPE z-stack also gets the molecules
per z plane and a recommended plane range (zplanes.py).

The config may carry `_cores`, the cores the run was given (the Snakefile adds it): a blank `threads`
uses that.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import inspect_raw as idt  # noqa: E402
import run_config as mc  # noqa: E402
import reader_files as RF  # noqa: E402
import saw_flat as SF  # noqa: E402
import ssaft_images as IM  # noqa: E402
import zplanes as ZP  # noqa: E402

DEFAULTS = mc.DEFAULTS
PLATFORMS = tuple(DEFAULTS["platforms"])
# default nuclei element per platform (ssaft_steps.nuclei_element uses the same names)
NUCLEI_DEFAULT = {"xenium": "nucleus_boundaries", "stereoseq": "cells_polygons", "cosmx": "cell_centroids"}
# elements ssaft (or this pipeline) writes itself; never offered as input choices
OWN = {"ssaft_cells", "ssaft_downsampled_vectors", "ssaft_cell_boundaries", "ssaft_cell_boundaries_2d",
       "ssaft_cell_boundaries_3d", "ssaft_cell_labels_3d", "ssaft_celltype_map", "ssaft_celltype_map_rgb",
       "ssaft_mask_image", "ssaft_nuclei", "ssaft_zstack"}
# Visium HD bin sizes the pipeline offers (Space Ranger's square_002um / 008um / 016um), the default first
VISIUMHD_BINS = (2, 8, 16)
# Stereo-seq DNB pitch (µm): the bin1 table's bin size, the blur raster when no SSAF-T zarr records one
STEREOSEQ_DNB_UM = SF.DNB_UM
# smallest sampling distance (grid spacing), µm: the transcript platforms; the floor of Stereo-seq's,
# which is its bin size (vf_bin_um) or this, whichever is larger. Visium HD: the bin size
MIN_SAMPLING_UM = 1.0
# anchors sampled per 10 µm block on the transcript platforms (n_samples blank), by the map's
# dimension: 2D / 3D (per imaging plane). SSAF-T's own defaults
N_SAMPLES = {2: 10, 3: 5}
# anchors per 10 µm block on the count platforms, always 2D (n_samples_per_block blank)
N_SAMPLES_PER_BLOCK = 10
KERNELS = ("gaussian", "tophat", "epanechnikov", "exponential", "linear", "cosine")   # ssaft._count.KDE_KERNELS
# platforms whose z-stack is stacked from image planes, so that z planes can be chosen
ZSTACK_PLATFORMS = ("xenium", "merscope")
# config keys that only a 3D map uses
ZKEYS = ("zstack_source", "z_step", "z_planes", "merscope_z_plane_um", "cosmx_z_plane_um")
# config keys that only downsampling on a z-stack uses
STACK_KEYS = ("zstack_source", "z_step", "z_planes")
# run_config.validate's lines about a z-stack, not wanted when a 3D map downsamples on a 2D image
STACK_ERRORS = ("3D needs a z-stack source", "z-stack source not found", "3D needs the z plane spacing")
# label-list files CeLLMan's annotate-reference-mode reads
REF_LIST_SUFFIXES = (".json", ".csv", ".tsv")
# the MERSCOPE mosaic z layer spatialdata_io loads as the 2D image (merscope_mask_z blank)
_zl = DEFAULTS["platforms"]["merscope"]["load"]["z_layers"]
MERSCOPE_MASK_Z = int(_zl[0] if isinstance(_zl, list) else _zl)


class PlanError(Exception):
    pass


# ------------------------------------------------------------------------------------------ helpers
def blank(v) -> bool:
    return v is None or (isinstance(v, str) and v.strip() == "") or (isinstance(v, (list, tuple)) and not v)


def choices(C, key, allowed):
    if not blank(C.get(key)) and str(C[key]).strip() not in allowed:
        raise PlanError(f"{key} must be one of {', '.join(allowed)}, got {C[key]!r}")


def to_bool(v, key):
    if isinstance(v, bool):
        return v
    s = str(v).strip().lower()
    if s in ("true", "yes", "1", "on"):
        return True
    if s in ("false", "no", "0", "off"):
        return False
    raise PlanError(f"{key}: expected true or false, got {v!r}")


def to_num(v, key, kind=float):
    try:
        return kind(v)
    except (TypeError, ValueError):
        raise PlanError(f"{key}: expected a number, got {v!r}")


def input_path(v) -> str:
    s = str(v).strip()
    return s if s.startswith("/") else "/input/" + s


def read_json(p: Path):
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def zattrs(d: Path) -> dict:
    """Attributes of a zarr group/array, v3 (zarr.json) or v2 (.zattrs)."""
    j = read_json(d / "zarr.json")
    if j is not None:
        return j.get("attributes") or {}
    return read_json(d / ".zattrs") or {}


def zshape(d: Path):
    j = read_json(d / "zarr.json")
    if j is not None and "shape" in j:
        return j["shape"]
    j = read_json(d / ".zarray")
    return j.get("shape") if j else None


def gpus_visible() -> list[str]:
    try:
        r = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30)
        return [l.strip() for l in r.stdout.splitlines() if l.strip()] if r.returncode == 0 else []
    except Exception:  # noqa: BLE001
        return []


# --------------------------------------------------------------------------------- zarr inspection
def is_zarr(p: Path) -> bool:
    if not p.is_dir():
        return False
    if p.suffix == ".zarr":
        return True
    has_meta = (p / "zarr.json").exists() or (p / ".zgroup").exists()
    return has_meta and any((p / k).is_dir() for k in ("images", "points", "shapes", "tables", "labels"))


def image_info(d: Path) -> dict:
    a = zattrs(d)
    ome = a.get("ome") or a
    ms = ome.get("multiscales") or [{}]
    axes = [x.get("name") if isinstance(x, dict) else x for x in ms[0].get("axes", [])]
    ds = ms[0].get("datasets") or [{}]
    shape = zshape(d / ds[0]["path"]) if ds and ds[0].get("path") else None
    chans = [c.get("label") for c in (ome.get("omero") or {}).get("channels", []) if isinstance(c, dict)]
    return {"axes": axes, "shape": shape, "channels": chans, "is_3d": "z" in axes}


def points_info(d: Path) -> dict:
    cols = []
    try:
        import pyarrow.parquet as pq
        pp = d / "points.parquet"
        parts = sorted(pp.rglob("*.parquet")) if pp.is_dir() else [pp]
        if parts and parts[0].exists():
            cols = list(pq.read_schema(parts[0]).names)
    except Exception:  # noqa: BLE001
        pass
    gene = next((c for c in ("gene", "gene_symbol", "feature_name", "target") if c in cols), None)
    fk = (zattrs(d).get("spatialdata_attrs") or {}).get("feature_key")
    return {"columns": cols, "has_z": "z" in cols, "gene_column": fk or gene}


def table_info(d: Path) -> dict:
    info = {"shape": zattrs(d / "X").get("shape"), "cellman": (d / "uns" / "cellman").is_dir(), "meta": None}
    if info["cellman"]:
        try:
            import zarr
            try:
                from anndata.io import read_elem
            except ImportError:  # anndata < 0.11
                from anndata.experimental import read_elem
            meta = read_elem(zarr.open_group(str(d / "uns" / "cellman"), mode="r"))
            info["meta"] = json.loads(json.dumps(meta, default=lambda o: o.tolist()
                                                 if hasattr(o, "tolist") else str(o)))
        except Exception as e:  # noqa: BLE001
            info["meta_error"] = f"{type(e).__name__}: {e}"
    return info


def inspect_zarr(p: Path) -> dict:
    el = {}
    for kind in ("images", "labels", "points", "shapes", "tables"):
        d = p / kind
        el[kind] = sorted(x.name for x in d.iterdir() if x.is_dir() and not x.name.startswith(".")) \
            if d.is_dir() else []
    root = zattrs(p)
    return {"elements": el,
            "root_attrs": {k: v for k, v in root.items() if k != "spatialdata_attrs"},
            "images": {k: image_info(p / "images" / k) for k in el["images"]},
            "points": {k: points_info(p / "points" / k) for k in el["points"]},
            "tables": {k: table_info(p / "tables" / k) for k in el["tables"]}}


def detect_zarr_platform(z: dict):
    """Platform of a SpatialData zarr from its element names; the reader attribute only breaks ties.

    Returns (platform or None, how it was decided, candidate platforms)."""
    el = z["elements"]
    I, P, S, T = el["images"], el["points"], el["shapes"], el["tables"]
    cellman = [k for k, t in z["tables"].items() if t["cellman"]]
    if cellman:
        meta = z["tables"][cellman[0]].get("meta") or {}
        plat = str(meta.get("platform", "")).lower().replace("_", "").replace("-", "")
        if plat in ("visiumhd", "stereoseq"):
            return plat, f"SSAF-T count zarr: tables[{cellman[0]!r}].uns['cellman']['platform'] = {plat!r}", [plat]
        return None, (f"SSAF-T count zarr (tables[{cellman[0]!r}] carries uns['cellman']) but its platform "
                      f"is {meta.get('platform')!r}"), ["visiumhd", "stereoseq"]
    hits = {}
    if any(re.fullmatch(r"square_\d{3}um", t) for t in T) or \
            any(i.endswith(("_full_image", "_hires_image", "_lowres_image")) for i in I):
        hits["visiumhd"] = "a square_###um bin table / *_full_image or *_hires_image image"
    if any(i.endswith("_regist") for i in I) or any(re.fullmatch(r"bin\d+_table", t) for t in T) \
            or "cells_polygons" in S:
        hits["stereoseq"] = "a *_regist image / bin<N>_table / cells_polygons"
    if any(i.startswith("morphology_") for i in I) or "nucleus_boundaries" in S:
        hits["xenium"] = "morphology_* images / nucleus_boundaries shapes"
    if any(k.endswith("_transcripts") for k in P) and any(re.search(r"_z\d+$", i) for i in I):
        hits["merscope"] = "a <dataset>_transcripts points element with <dataset>_z<N> mosaic images"
    if any(re.fullmatch(r"\d+_points", k) for k in P) or any(re.fullmatch(r"\d+_image", i) for i in I):
        hits["cosmx"] = "per-FOV <fov>_points / <fov>_image elements"
    elif "transcripts" in P and "composite" in I:
        hits["cosmx"] = "a transcripts points element with a composite image (whole-slide layout)"
    if len(hits) == 1:
        k = next(iter(hits))
        return k, f"element names: {hits[k]}", [k]
    reader = {"xenium": "xenium", "merscope": "merscope", "cosmx": "cosmx", "visium_hd": "visiumhd",
              "stereoseq": "stereoseq"}.get(z["root_attrs"].get("spatialdata_io_reader"))
    if reader and (not hits or reader in hits):
        return reader, f"root attribute spatialdata_io_reader = {z['root_attrs'].get('spatialdata_io_reader')!r}", [reader]
    how = ("element names match " + ", ".join(f"{k} ({v})" for k, v in hits.items())) if hits \
        else "no element names of a known platform"
    return None, how + ", and the zarr names no reader", sorted(hits) or list(PLATFORMS)


def stack_planes(kind, cfg, p: Path) -> list[int]:
    """Plane numbers of the run's z-stack, from its metadata only (no molecules read)."""
    element = (cfg.get("zarr") or {}).get("zstack_element") if kind == "zarr" else None
    if element:
        return [k for k, _ in ZP.element_planes(p, element)[0]]
    return [k for k, _ in ZP.file_planes(cfg["zstack"]["source"], cfg["zstack"].get("z_step"))[0]]


def merscope_stains(kind, insp, z) -> list[str]:
    """The stains a MERSCOPE dataset has: from the mosaic file names, or the channels of a zarr's z mosaic."""
    if kind == "raw":
        return list((insp.get("hints") or {}).get("stains") or [])
    for k, v in z["images"].items():
        if re.search(r"_z\d+$", k) and not v["is_3d"] and v.get("channels"):
            return [str(c) for c in v["channels"]]
    return []


def merscope_layers(kind, p: Path, z) -> dict:
    """The z layers that have a 2D mosaic. raw: {layer: stains with a mosaic_<stain>_z<layer>.tif};
    zarr: {layer: the <dataset>_z<layer> 2D image element}."""
    out = {}
    if kind == "raw":
        for f in (p / "images").glob("mosaic_*_z*.tif") if (p / "images").is_dir() else []:
            m = re.match(r"mosaic_(.+)_z(\d+)\.tif$", f.name)
            if m:
                out.setdefault(int(m.group(2)), set()).add(m.group(1))
        return out
    for k, v in z["images"].items():
        m = re.search(r"_z(\d+)$", k)
        if m and not v["is_3d"] and k not in OWN:
            out[int(m.group(1))] = k
    return out


def count_zarr_image(z, meta) -> str | None:
    """The image an SSAF-T count zarr names in uns['cellman'] (the one SSAF-T reads), 'staining' if unnamed."""
    if not any(t["cellman"] for t in z["tables"].values()):
        return None
    return str(meta.get("image_key") or (meta.get("meta") or {}).get("image_key") or "staining")


def label_row(table: str, row: str, text: str) -> str:
    """The summary table with a note on one row (display only: run.json is unchanged)."""
    out = []
    for line in table.splitlines():
        if line.split("|")[0].strip() == row:
            line += f"  -- {text}"
        out.append(line)
    return "\n".join(out)


def file_stack_planes(source):
    """Plane numbers of a z-stack file / glob, from its metadata; None when it cannot be read."""
    try:
        return [k for k, _ in ZP.file_planes(source, None)[0]]
    except Exception:  # noqa: BLE001
        return None


# -------------------------------------------------------------------------------------- main plan
def plan(C: dict, run_dir: str) -> dict:
    # pending: values the user has to give, asked for after the whole plan made from the data
    notes, summary, pending = [], [], []
    z_um_from = None
    if blank(C.get("input")):
        inp = Path("/input")
        raise PlanError("input is required: the dataset folder or SpatialData .zarr under /input. /input holds: "
                        f"{sorted(x.name for x in inp.iterdir())[:40] if inp.is_dir() else '(not mounted)'}")
    p = Path(input_path(C["input"])).expanduser()
    if not p.exists():
        parent = p.parent if p.parent.exists() else Path("/input")
        listing = sorted(x.name for x in parent.iterdir())[:40] if parent.is_dir() else []
        raise PlanError(f"input {str(p)!r} does not exist. {parent} holds: {listing}")
    p = p.resolve()
    override = None if blank(C.get("platform")) else str(C["platform"]).strip().lower()
    if override and override not in PLATFORMS:
        raise PlanError(f"platform must be one of {list(PLATFORMS)} (or blank to detect), got {C['platform']!r}")

    # ---- inspection
    z, cands = None, []
    if is_zarr(p):
        kind = "zarr"
        z = inspect_zarr(p)
        plat, how, cands = detect_zarr_platform(z)
        if override:
            plat, how = override, f"set in config (detected: {how})"
        if plat is None:
            raise PlanError(f"could not tell the platform of {p}: {how}. Set platform to one of {cands}. "
                            f"Elements: {json.dumps(z['elements'])}")
        P = DEFAULTS["platforms"][plat]
        img3d = [k for k, v in z["images"].items() if v["is_3d"] and k not in OWN]
        insp = {"platform": plat, "supported": True, "path": str(p), "input_kind": "zarr",
                "dimensions": P["dimensions"], "experimental": P.get("experimental", False),
                "zstack": bool(img3d), "files": {}, "hints": {}, "zarr": z, "detected_by": how}
    else:
        kind = "raw"
        rep = None
        for fn in (idt.xenium, idt.merscope, idt.cosmx, idt.visiumhd, idt.stereoseq):
            rep = fn(p)
            if rep:
                break
        how = "inspect_raw.py"
        if rep is not None and rep.get("platform") == "stereoseq" and not rep.get("supported", True) \
                and SF.is_flat(p):
            # a SAW output in one folder: spatialdata_io.stereoseq does not read it, saw_flat.py does
            rep = SF.inspect(p)
            how = ("inspect_raw.py: a Stereo-seq SAW output with every file in one folder, which "
                   "spatialdata_io.stereoseq does not read; read by the pipeline's flat SAW loader (saw_flat.py)")
        if rep is None:
            miss = RF.near_miss(p)
            if miss:
                raise PlanError(miss)
            raise PlanError(f"{p}: no Xenium / MERSCOPE / CosMx / Visium HD / Stereo-seq layout recognised and it "
                            f"is not a SpatialData .zarr. It holds: {sorted(x.name for x in p.iterdir())[:60]}")
        if not rep.get("supported", True):
            raise PlanError(f"{p} ({rep['platform']}): {rep.get('problem')}")
        rep["path"] = str(p)
        insp = dict(rep, input_kind="raw", detected_by=how)
        plat = rep["platform"]
        if override and override != plat:
            raise PlanError(f"platform is set to {override!r} but {p} is a {plat} output folder")
        cands = [plat]
        P = DEFAULTS["platforms"][plat]
    flat = insp.get("layout") == "saw_flat"
    loader = (RF.LOADER["stereoseq_flat" if flat else plat] if kind == "raw" else
              "spatialdata.read_zarr: the SpatialData zarr is copied into the run (no spatialdata_io reader)")

    # ---- dimension = the map's, 2D or 3D. Visium HD / Stereo-seq: always 2. Xenium, MERSCOPE, CosMx:
    # blank = 3. The image downsampling uses is chosen apart from it (image_key, below)
    regime, dims = P["regime"], P["dimensions"]
    if blank(C.get("dimension")):
        dim = 2 if regime == "count" else 3
    else:
        dim = to_num(C["dimension"], "dimension", int)
    if dim not in dims:
        if regime == "count":
            # Visium HD and Stereo-seq: counts on a 2D bin lattice, no z (2D/3D is never asked)
            raise PlanError(f"dimension {dim} is not possible for {plat}: Visium HD and Stereo-seq are count "
                            f"platforms (bins on a 2D lattice, no z) and map in 2D only. Set dimension: 2 or "
                            f"leave it blank.")
        raise PlanError(f"{plat} maps in {dims}D only; dimension {dim} is not possible")
    if dim == 3 and kind == "zarr" and regime == "transcript":
        tx = {k: v for k, v in z["points"].items() if v.get("gene_column")}
        if tx and not any(v["has_z"] for v in tx.values()):
            raise PlanError(f"a 3D map needs the transcripts' z, and the zarr's transcripts {sorted(tx)} have no z "
                            f"column. Set dimension: 2 (a 2D map).")
    if dim == 2:
        unused = [k for k in ZKEYS if not blank(C.get(k)) and not (plat == "merscope" and k == "z_planes")]
        if unused:
            notes.append(f"{', '.join(unused)} set but not used: this map is 2D (no z-stack, no z)")
    try:
        zsel = ZP.parse(C.get("z_planes"))
    except ZP.SelectionError as e:
        raise PlanError(str(e))

    # ---- answers (run_config.build input); only what config.yaml sets, the rest stays default
    cores = int(C.get("_cores") or DEFAULTS["common"]["ncores"])
    ncores = min(to_num(C["threads"], "threads", int), cores) if not blank(C.get("threads")) else cores
    A: dict = {"run_dir": run_dir, "python": sys.executable, "dimension": dim,
               "condition": "normal" if blank(C.get("condition")) else str(C["condition"]).strip(),
               "ncores": ncores}

    def sec(s):
        return A.setdefault(s, {})

    def put(section, key, ckey, conv=None):
        v = C.get(ckey)
        if not blank(v):
            sec(section)[key] = conv(v) if conv else (v.strip() if isinstance(v, str) else v)

    def num(k):
        return lambda v: to_num(v, k)

    def integer(k):
        return lambda v: to_num(v, k, int)

    def boolean(k):
        return lambda v: to_bool(v, k)

    # tissue: required; an SSAF-T count zarr records it, so it may come from there. When it is not
    # known, the rest of the plan is still made from the data and shown before the request for it
    species, organ = C.get("species"), C.get("organ")
    meta = {}
    if kind == "zarr":
        meta = next((t["meta"] for t in z["tables"].values() if t.get("meta")), None) or {}
        if blank(species) and meta.get("species") not in (None, "", "unknown"):
            species = meta["species"]
            notes.append(f"species {species!r} taken from the zarr's uns['cellman']")
        if blank(organ) and meta.get("organ") not in (None, "", "unknown"):
            organ = meta["organ"]
            notes.append(f"organ {organ!r} taken from the zarr's uns['cellman']")
    tissue_missing = None
    if blank(species) or blank(organ):
        try:
            mods = [f"{w['species']}/{w['organ']}" for w in
                    mc.modules(DEFAULTS["cellman_server"], None if blank(species) else str(species).strip())]
        except Exception as e:  # noqa: BLE001
            mods = [f"(module list unavailable: {e})"]
        tissue_missing = (f"species and organ are required (a Cell-Ontology module, e.g. mouse / kidney); "
                          f"{'species' if blank(species) else 'organ'}{' and organ' if blank(species) and blank(organ) else ''}"
                          f" {'are' if blank(species) and blank(organ) else 'is'} blank. Decide them with the user "
                          f"from what the data shows above and set them. Modules served: {mods}")
        A["species"] = "" if blank(species) else str(species).strip()
        A["organ"] = "" if blank(organ) else str(organ).strip()
    else:
        A["species"], A["organ"] = str(species).strip(), str(organ).strip()
        try:
            served = {w["organ"] for w in mc.modules(DEFAULTS["cellman_server"], A["species"])}
        except Exception as e:  # noqa: BLE001
            served = None
            notes.append(f"could not check the organ against the CeLLMan server ({e}); it is checked again later")
        if served is not None:
            if not served:
                raise PlanError(f"no Cell-Ontology modules for species {A['species']!r} (use mouse or human)")
            if A["organ"] not in served:
                raise PlanError(f"organ {A['organ']!r} is not a {A['species']} Cell-Ontology module. Served: {sorted(served)}")

    # compute: the device is resolved at run time by the stage driver; only the request is recorded
    gpu_req = "auto" if blank(C.get("gpu")) else str(C["gpu"]).strip().lower()
    if gpu_req not in ("auto", "cpu", "cuda"):
        raise PlanError(f"gpu must be auto, cpu or cuda, got {C['gpu']!r}")
    A["gpu_device"] = None if gpu_req == "cpu" else "cuda:0"

    # platform input
    if plat == "visiumhd" and kind == "raw" and not blank(C.get("visiumhd_fullres_image")):
        v = str(C["visiumhd_fullres_image"]).strip()
        sec("load")["fullres_image_file"] = None if v.lower() == "none" else input_path(v)
    bin_um = None
    if plat == "visiumhd":
        # the bin size the run uses: the raw reader loads only that square_###um table, and a zarr keeps
        # only that one at load (ssaft_load.py); an SSAF-T count zarr has a single bin size, recorded
        recorded = (meta.get("meta") or {}).get("bin_um") if kind == "zarr" else None
        cellman_zarr = kind == "zarr" and any(t["cellman"] for t in z["tables"].values())
        asked = None if blank(C.get("visiumhd_bin_um")) else to_num(C["visiumhd_bin_um"], "visiumhd_bin_um")
        if recorded:
            bin_um = float(recorded)
            if asked is not None and asked != bin_um:
                raise PlanError(f"visiumhd_bin_um is {asked:g} but this SSAF-T zarr holds {bin_um:g} µm bins only; "
                                f"leave visiumhd_bin_um blank (or {bin_um:g})")
            notes.append(f"Visium HD bin size {bin_um:g} µm, recorded in the zarr's uns['cellman']")
        else:
            b = VISIUMHD_BINS[0] if asked is None else asked
            if b not in VISIUMHD_BINS:
                raise PlanError(f"visiumhd_bin_um must be one of {list(VISIUMHD_BINS)} (µm), got {C['visiumhd_bin_um']!r}")
            b = int(b)
            want = f"square_{b:03d}um"
            if not cellman_zarr:
                have = (sorted(t for t in z["elements"]["tables"] if re.search(r"square_\d{3}um", t)) if kind == "zarr"
                        else sorted(x.name for x in (p / "binned_outputs").glob("square_*um") if x.is_dir()))
                if not any(want in h for h in have):
                    raise PlanError(f"visiumhd_bin_um {b}: this dataset has no {want} bins; it has {have or 'none'}. "
                                    f"Set visiumhd_bin_um to one of the sizes it has.")
            bin_um = float(b)
            sec("load")["bin_size"] = b
    elif not blank(C.get("visiumhd_bin_um")):
        notes.append(f"visiumhd_bin_um is set but applies to Visium HD only; ignored for {plat}")
    px_line = px_used = px_src = None
    cosmx_z_missing = False
    if plat == "cosmx":
        # µm per image pixel. A raw folder: the load takes it from the FOV image size
        # (ssaft_load.cosmx_px_um) unless load.px_um_override is set; cosmx_px_um sets that override
        # when it differs from what the image size gives (so a value equal to it changes nothing)
        asked = None if blank(C.get("cosmx_px_um")) else to_num(C["cosmx_px_um"], "cosmx_px_um")
        if kind == "raw":
            w, h, px_img, cam = RF.cosmx_camera(p)
            if asked is not None:
                sec("load")["px_um"] = asked
                px_used, px_src = asked, "cosmx_px_um"
                if px_img is not None and abs(asked - px_img) < 1e-9:
                    px_line = f"{asked:g} µm/px (cosmx_px_um; the same as the FOV image size gives: {cam})"
                else:
                    sec("load")["px_um_override"] = asked
                    px_line = (f"{asked:g} µm/px (cosmx_px_um, which overrides the FOV image size"
                               + (f"; the {w} x {h} images would give {px_img:g} µm/px)" if px_img is not None
                                  else f"; {cam})"))
            elif px_img is not None:
                px_used, px_src = px_img, f"from the FOV image size, {w} x {h} px"
                px_line = f"{px_img:g} µm/px, from the FOV image size: {cam} (cosmx_px_um blank)"
            elif w is not None:
                raise PlanError(f"cosmx_px_um: the CellComposite images are {w} x {h} px, not a known CosMx camera "
                                f"(5472 x 3648 px is 0.18 µm/px, 4256 x 4256 px is 0.12028 µm/px), so the pixel "
                                f"size cannot be told from them. Ask the user for the instrument's µm per pixel and "
                                f"set cosmx_px_um.")
        else:
            put("load", "px_um", "cosmx_px_um", num("cosmx_px_um"))
            px_line = (f"{asked:g} µm/px (cosmx_px_um)" if asked is not None else
                       f"{DEFAULTS['platforms']['cosmx']['load']['px_um']:g} µm/px (cosmx_px_um blank: the "
                       f"default for a zarr, which carries no FOV image file to size)")
        if not blank(C.get("cosmx_fovs")) and str(C["cosmx_fovs"]).strip().lower() != "all":
            sec("load")["fovs"] = [int(x) for x in re.split(r"[,\s]+", str(C["cosmx_fovs"]).strip()) if x]
        if dim == 3:
            if blank(C.get("cosmx_z_plane_um")):
                # asked for after the plan, as the MERSCOPE spacing is
                cosmx_z_missing = True
                pending.append("cosmx_z_plane_um is required for a 3D map on CosMx: the transcripts' z is a plane "
                               "number, and the µm between planes is not recorded in the flat files. Ask the user "
                               "for this experiment's z plane spacing and set cosmx_z_plane_um; do not propose a "
                               "value of your own. A 2D map needs no spacing.")
            else:
                sec("load")["z_plane_um"] = to_num(C["cosmx_z_plane_um"], "cosmx_z_plane_um")
    elif not blank(C.get("cosmx_px_um")):
        notes.append(f"cosmx_px_um is set but applies to CosMx only; ignored for {plat}")

    # MERSCOPE: every stain is loaded; the one chosen (merscope_stain) makes the tissue mask and the z-stack
    stains, stain, mask_z, mask_layers, mname = [], None, None, [], None
    z_um_missing = False
    if plat == "merscope":
        stains = merscope_stains(kind, insp, z)
        stain = None if blank(C.get("merscope_stain")) else str(C["merscope_stain"]).strip()
        if stain is None:
            if not stains or "DAPI" in stains:
                stain = "DAPI"
            else:
                raise PlanError(f"merscope_stain: choose the stain for the tissue mask and the z-stack, one of "
                                f"{stains} (this dataset has no DAPI)")
        if stains and stain not in stains:
            raise PlanError(f"merscope_stain {stain!r} is not in this dataset; its stains are {stains}")
        # the z layer spatialdata_io loads as the 2D mosaic, the MERSCOPE folder's 2D image
        mz_set = not blank(C.get("merscope_mask_z"))
        layers = merscope_layers(kind, p, z)
        mask_layers = sorted(k for k, s in layers.items() if stain in s) if kind == "raw" else sorted(layers)
        if kind == "zarr":
            # a zarr's images are chosen by name (image_key); no layer has to exist
            mz = to_num(C["merscope_mask_z"], "merscope_mask_z", int) if mz_set else MERSCOPE_MASK_Z
            mask_z = mz if mz in mask_layers else None
            if mz_set and mask_z is None:
                notes.append(f"merscope_mask_z {mz}: the zarr has no z{mz} mosaic; not used by this run")
        else:
            mask_z = to_num(C["merscope_mask_z"], "merscope_mask_z", int) if mz_set else MERSCOPE_MASK_Z
            if mask_z not in mask_layers:
                given = "merscope_mask_z" if mz_set else "merscope_mask_z is blank, so the default z layer"
                raise PlanError(
                    f"{given} {mask_z}: {p / 'images'} has no mosaic_{stain}_z{mask_z}.tif. The z layers with a "
                    f"{stain} mosaic are {mask_layers or 'none'}. Set merscope_mask_z to one of them.")
            sec("load")["z_layers"] = mask_z
            lacking = [s for s in stains if s not in layers[mask_z]]
            if lacking:
                stains = [s for s in stains if s in layers[mask_z]]
                notes.append(f"stains without a z{mask_z} mosaic are not loaded: {lacking}")
            if not mz_set:
                notes.append(f"merscope_mask_z blank: the z{mask_z} mosaic (the layer spatialdata_io loads by default)")
            # the elements are named after the data folder (ssaft_images.rename_merscope)
            mname = IM.folder_name(insp["path"])
            sec("image")["key"] = f"{mname}_z{mask_z}"
            if stains:
                sec("load")["stains"] = stains
        sec("zstack")["stain"] = stain
        if stains:
            sec("image")["channel"] = stain
        else:
            notes.append(f"the zarr's mosaic records no channel names, so the tissue mask uses all its channels "
                         f"(merscope_stain {stain!r} is used for nothing else)")
        if dim == 3:
            # a 3D map needs the transcripts' z in µm, whichever image is used
            header_um, header_from = (ZP.header_spacing(p / "images", stain)
                                      if blank(C.get("merscope_z_plane_um")) and kind == "raw" else (None, None))
            if header_um:
                sec("load")["z_plane_um"] = header_um
                z_um_from = header_from
            elif blank(C.get("merscope_z_plane_um")):
                z_um_missing = True
                pending.append(f"merscope_z_plane_um is required for a 3D map on MERSCOPE: the transcripts' z is a "
                               f"plane number and {header_from or 'a zarr records no spacing'} (the z plane table "
                               f"above, if any, is by plane number). Offer the user the choice: the spacing of "
                               f"their experiment, or the usual MERSCOPE spacing of 1.5 µm. Then set "
                               f"merscope_z_plane_um.")
            else:
                sec("load")["z_plane_um"] = to_num(C["merscope_z_plane_um"], "merscope_z_plane_um")
    else:
        for k in ("merscope_stain", "merscope_mask_z"):
            if not blank(C.get(k)):
                notes.append(f"{k} is set but applies to MERSCOPE only; ignored for {plat}")

    # ---- the image downsampling uses (transcript platforms): candidates, recommendation, choice
    src = None if blank(C.get("zstack_source")) else str(C["zstack_source"]).strip()
    if plat == "cosmx" and kind == "raw" and src:
        raise PlanError("zstack_source is set for a CosMx folder: CosMx has no z-stack here; a 3D map uses the 2D "
                        "composite image, and only the cell-type map is 3D. Leave zstack_source blank")
    per_fov = kind == "zarr" and plat == "cosmx" and any(re.fullmatch(r"\d+_image", i) for i in z["elements"]["images"])
    images, used, uses_stack, file_stack_name = None, None, False, None
    if regime == "transcript" and not per_fov:
        if kind == "zarr":
            rows = IM.zarr_images(z, OWN)
            if src is not None and src not in z["elements"]["images"]:
                if plat == "merscope":
                    raise PlanError("zstack_source names a file for a MERSCOPE zarr; choose a z-stack of the zarr "
                                    "with image_key instead")
                file_stack_name = IM.zstack_key(src)
                pl = file_stack_planes(input_path(src))
                rows.append({"name": file_stack_name, "dims": 3, "axes": "czyx", "shape": None, "channels": [],
                             "planes": pl, "origin": f"z-stack: zstack_source {input_path(src)}"})
        else:
            if plat == "merscope":
                if src is not None:
                    raise PlanError("zstack_source is set for a MERSCOPE folder: its z-stack is built from its own "
                                    "mosaics (the planes chosen in z_planes); leave zstack_source blank")
                stack_planes_raw, source = mask_layers, f"images/mosaic_{stain}_z*.tif"
                file_stack_name = IM.zstack_key(None, stain)
            else:
                path = input_path(src) if src is not None else (
                    str(p / "morphology.ome.tif") if plat == "xenium" and (p / "morphology.ome.tif").exists() else None)
                source = (src if src is not None else "morphology.ome.tif")
                stack_planes_raw = file_stack_planes(path) if path else None
                file_stack_name = IM.zstack_key(path) if path else None
            rows = IM.raw_image_rows(plat, insp["path"], insp, stack_planes_raw, zsel, name=mname, stains=stains,
                                     stain=stain, mask_z=mask_z, source=source, stack_name=file_stack_name)
        rec, why = IM.recommend(plat, dim, rows, mask_z, stain)
        asked = C.get("image_key")
        if kind == "zarr" and src is not None and src in z["elements"]["images"]:
            # zstack_source naming an element of the zarr is that element chosen
            if not blank(asked) and str(asked).strip() != src:
                raise PlanError(f"image_key {str(asked).strip()!r} and zstack_source {src!r} name different "
                                f"images; set only image_key")
            asked = src
        try:
            used, how = IM.choose(asked, dim, rows, rec)
        except IM.ImageError as e:
            raise PlanError(str(e))
        if used is None:
            raise PlanError(f"there is no image to make the tissue mask from (images: {[r['name'] for r in rows]})")
        urow = next(r for r in rows if r["name"] == used)
        uses_stack = dim == 3 and urow["dims"] == 3
        images = {"rows": rows, "used": used, "kind": "zstack" if uses_stack else "2d", "recommended": rec,
                  "why": why, "how": how, "dims": dim}
        if kind == "zarr":
            if uses_stack:
                two = IM.recommend(plat, 2, rows, mask_z, stain)[0]
                if two:
                    sec("image")["key"] = two       # a 2D image is kept too: tiles, pixel grid
            else:
                sec("image")["key"] = used
    elif per_fov and not blank(C.get("image_key")):
        pick = str(C["image_key"]).strip()
        if pick not in z["elements"]["images"]:
            raise PlanError(f"image_key {pick!r} is not in the zarr; images: {z['elements']['images']}")
        sec("image")["key"] = pick
    elif regime == "count" and not blank(C.get("image_key")):
        notes.append(f"image_key is set but not used: on {plat} SSAF-T reads the count data's own image")

    # z-stack: only when the image used is one (a 3D map on a z-stack)
    zstack_element = None
    if dim == 3:
        if kind == "zarr" and plat == "cosmx" and (per_fov or uses_stack):
            raise PlanError("CosMx 3D on a zarr's z-stack (or a per-FOV zarr) is not supported; use the zarr's 2D "
                            "composite image (image_key) for a 3D map")
        if uses_stack:
            if kind == "zarr" and used != file_stack_name:
                zstack_element = used
                sec("zstack")["source"] = str(p / "images" / used)
                sec("zstack")["element"] = used
            elif src is not None:
                sec("zstack")["source"] = input_path(src)
            put("zstack", "z_step", "z_step", num("z_step"))
            if plat == "merscope" and blank(C.get("z_step")) and "z_plane_um" in A.get("load", {}):
                # the mosaic planes are the transcripts' z planes: same spacing
                sec("zstack")["z_step"] = A["load"]["z_plane_um"]
        elif not per_fov:
            if zsel:
                raise PlanError(f"z_planes is set, which chooses the planes of a z-stack, but the image this 3D map "
                                f"downsamples on is the 2D image {used!r} (one plane, placed at the middle of the "
                                f"transcripts' z range). Set image_key to the z-stack to use its planes, or clear "
                                f"z_planes.")
            unused = [k for k in STACK_KEYS if not blank(C.get(k)) and k != "z_planes"]
            if unused:
                notes.append(f"{', '.join(unused)} set but not used: this 3D map downsamples on the 2D image {used!r}")
    elif plat == "merscope" and zsel:
        raise PlanError("z_planes is set, which chooses the planes of a z-stack, but the map is 2D (dimension 2), "
                        "which downsamples on a 2D image. Set dimension: 3 and image_key to the z-stack to use those "
                        "planes, or clear z_planes.")

    # tiling: Xenium and MERSCOPE; default tiled with a
    # z-stack. The tiles carry the z-stack, so a 3D map on a 2D image is not tiled
    if plat in ("xenium", "merscope"):
        tile_on = uses_stack if blank(C.get("tile")) else to_bool(C["tile"], "tile")
        if tile_on and dim == 3 and not uses_stack:
            raise PlanError(f"tile: a tiled 3D map carries a z-stack in every tile, and this map downsamples on the 2D "
                            f"image {used!r}. Leave tile blank (or false), or downsample on the z-stack.")
        sec("tile")["enabled"] = tile_on
        put("tile", "tile_size", "tile_size", num("tile_size"))
    elif not blank(C.get("tile")) and to_bool(C["tile"], "tile"):
        raise PlanError(f"tiling is available for Xenium and MERSCOPE only, not {plat}")

    # image / anchors
    mask_line = None
    flat_bin = None
    if regime == "transcript":
        if images:
            if uses_stack:
                mask_line = (f"{used} ({images['how']}): a z-stack, one tissue mask per plane (the planes chosen in "
                             f"z_planes)")
            elif dim == 3:
                mask_line = (f"{used} ({images['how']}): a 2D image in a 3D map, one tissue mask, placed at the "
                             f"middle of the transcripts' z range")
            else:
                mask_line = f"{used} ({images['how']}): a 2D image"
        elif per_fov:
            mask_line = "each FOV's own <fov>_image" + (f" ({sec('image')['key']}, image_key)"
                                                        if not blank(C.get("image_key")) else "")
        if plat == "merscope":
            mask_line += f"; stain {stain} (merscope_stain; stains loaded: {stains or 'unknown'})"
            if mask_layers:
                mask_line += f"; z layers with a mosaic: {mask_layers}"
        else:
            put("image", "channel", "image_channel")
            if not blank(C.get("image_channel")):
                mask_line += f"; channel {C['image_channel']} (image_channel)"
        put("image", "revert_image", "revert_image", boolean("revert_image"))
        # anchors per 10 µm block of each imaging plane: 2D map 10, 3D map 5, unless n_samples is set
        sec("downsample")["n_samples"] = (N_SAMPLES[dim] if blank(C.get("n_samples"))
                                          else to_num(C["n_samples"], "n_samples", int))
        # the Otsu threshold as computed unless otsu_mult is set (no platform-specific factor)
        sec("downsample")["otsu_mult"] = None if blank(C.get("otsu_mult")) else to_num(C["otsu_mult"], "otsu_mult")
        put("downsample", "bandwidth", "bandwidth", num("bandwidth"))
        put("downsample", "kernel", "kernel")
    else:
        if not blank(C.get("image_channel")):
            notes.append(f"image_channel is set but not used: on {plat} SSAF-T reads the count data's own image")
        own = count_zarr_image(z, meta) if kind == "zarr" else None
        flat_imgs = SF.images(p) if flat else []
        if flat:
            fimg = flat_imgs[0] if flat_imgs else None
            if fimg is None:
                mask_line = "the registered image <chip>_regist.tif (missing: see the files check)"
            else:
                mask_line = (f"{fimg.name}, registered to the DNB grid (1 px = 1 DNB = 0.5 µm): "
                             + ("H&E, dark tissue on bright (he_revert true)" if SF.is_he(fimg) else
                                "ssDNA, bright tissue on dark (he_revert false)")
                             + (f"; other registered images not used: {[f.name for f in flat_imgs[1:]]}"
                                if len(flat_imgs) > 1 else ""))
        elif own:
            there = "an image element of this zarr" if own in z["elements"]["images"] else \
                f"not among this zarr's images {z['elements']['images']}"
            mask_line = f"{own}, the staining image this SSAF-T zarr names in uns['cellman'] ({there})"
        elif plat == "visiumhd":
            mask_line = "chosen by SSAF-T from the count data: the full-resolution H&E (or the hires image)"
        else:
            mask_line = "chosen by SSAF-T from the count data: the registered ssDNA image (*_regist)"
        # anchors per 10 µm block (always 2D): 10 unless n_samples_per_block is set
        sec("downsample")["n_samples_per_block"] = (
            N_SAMPLES_PER_BLOCK if blank(C.get("n_samples_per_block"))
            else to_num(C["n_samples_per_block"], "n_samples_per_block", int))
        # the Otsu threshold of the tissue image as computed unless otsu_mult is set (SSAF-T 23dd889
        # takes it on the count path too); only written when set, so blank changes nothing
        if not blank(C.get("otsu_mult")):
            sec("downsample")["otsu_mult"] = to_num(C["otsu_mult"], "otsu_mult")
            mask_line += f"; Otsu threshold x {sec('downsample')['otsu_mult']:g} (otsu_mult)"
        put("vf", "bandwidth_um", "vf_bandwidth_um", num("vf_bandwidth_um"))
        put("vf", "kernel", "vf_kernel")
        put("panel", "hvg_n", "hvg_n", integer("hvg_n"))
        gem_head = ((insp.get("hints") or {}).get("gem_header") or {}) if flat else {}
        data_bin = SF.data_bin_um(gem_head) if flat else None
        if not blank(C.get("vf_bin_um")):
            sec("vf")["bin_um"] = to_num(C["vf_bin_um"], "vf_bin_um")
        elif plat == "visiumhd":
            sec("vf")["bin_um"] = int(bin_um) if bin_um.is_integer() else bin_um     # the bins' own size
        elif flat:
            sec("vf")["bin_um"] = data_bin
            notes.append(f"vf_bin_um {data_bin:g} µm: the GEM's own bin size (BinSize {gem_head.get('BinSize') or 1} "
                         f"x {STEREOSEQ_DNB_UM:g} µm DNB)")
        else:
            rec_bin = (meta.get("meta") or {}).get("bin_um") if kind == "zarr" else None
            if rec_bin:
                sec("vf")["bin_um"] = float(rec_bin)
                notes.append(f"vf_bin_um {float(rec_bin):g} µm: the bin size recorded in this SSAF-T zarr "
                             f"(uns['cellman']['meta']['bin_um'])")
            else:
                sec("vf")["bin_um"] = STEREOSEQ_DNB_UM
                notes.append(f"vf_bin_um {STEREOSEQ_DNB_UM:g} µm: the bin1 table's bin size (the {STEREOSEQ_DNB_UM:g} µm "
                             f"DNB pitch)")
        if flat:
            # the load sums the DNBs into bins of vf_bin_um, and the vector field blurs those bins
            flat_bin = float(A["vf"]["bin_um"])
            steps = flat_bin / STEREOSEQ_DNB_UM
            if abs(steps - round(steps)) > 1e-9 or flat_bin < data_bin - 1e-9:
                raise PlanError(f"vf_bin_um {flat_bin:g}: the load sums the GEM's DNBs ({STEREOSEQ_DNB_UM:g} µm) into "
                                f"bins of vf_bin_um, so it must be a multiple of {STEREOSEQ_DNB_UM:g} µm and at least "
                                f"the GEM's own bin size, {data_bin:g} µm. Leave it blank for {data_bin:g} µm.")
            fimg = flat_imgs[0] if flat_imgs else None
            sec("load").update(layout="saw_flat", gem=insp["files"].get("gem"),
                               image=str(fimg) if fimg else None, cellbin=insp["files"].get("cellbin"),
                               bin_um=flat_bin, he_revert=bool(fimg is not None and SF.is_he(fimg)))
            offs = {k: gem_head.get(k) for k in ("OffsetX", "OffsetY") if str(gem_head.get(k) or "0") != "0"}
            if offs:
                notes.append(f"GEM {offs}: not applied; the GEM's x, y are used as they are, as the registered "
                             f"image and the cellbin.gef are in that frame")

    # sampling distance: the grid spacing of the vector field and of the cell-type map. It cannot be finer
    # than the bins the counts are measured on (Visium HD: the bin size; Stereo-seq: the blur raster bin,
    # vf_bin_um), nor finer than 1 µm (Stereo-seq whose bins are smaller, and the transcript platforms)
    if plat == "visiumhd":
        sd_min = bin_um
    elif plat == "stereoseq":
        sd_min = max(float(A["vf"]["bin_um"]), MIN_SAMPLING_UM)
    else:
        sd_min = MIN_SAMPLING_UM
    if blank(C.get("sampling_distance")):
        sdist = sd_min
    else:
        sdist = to_num(C["sampling_distance"], "sampling_distance")
        if sdist < sd_min:
            if plat == "visiumhd":
                why_sd = (f"it cannot be smaller than the Visium HD bin size used, {sd_min:g} µm (visiumhd_bin_um): "
                          f"the counts are measured per bin, so a finer grid has nothing to resolve")
            elif plat == "stereoseq":
                vb = float(A["vf"]["bin_um"])
                why_sd = (f"the minimum for Stereo-seq is the bin size the counts are on, {vb:g} µm (vf_bin_um), or "
                          f"1 µm when the bins are smaller: max({vb:g}, 1) = {sd_min:g} µm. A finer grid has nothing "
                          f"to resolve below the bins, and below 1 µm multiplies the memory and run time")
            else:
                why_sd = (f"the minimum for {plat} is {sd_min:g} µm: a finer grid multiplies the memory and run time "
                          f"of the vector field and the map")
            raise PlanError(f"sampling_distance {sdist:g} µm is too small: {why_sd}. Set sampling_distance to "
                            f"{sd_min:g} or more, or leave it blank ({sd_min:g} µm).")
    sec("vf" if regime == "count" else "downsample")["sampling_distance"] = sdist
    sd_how = {"visiumhd": "the Visium HD bin size",
              "stereoseq": f"max(vf_bin_um {float((A.get('vf') or {}).get('bin_um') or 0):g}, 1)"}.get(plat, "1 µm")

    # clustering, mapping, filter
    put("cluster", "resolution", "resolution", num("resolution"))
    put("cluster", "n_neighbors", "n_neighbors", integer("n_neighbors"))
    put("map", "method", "map_method")
    put("filter", "norm_f", "norm_f", num("norm_f"))
    put("filter", "show_unknown", "show_unknown", boolean("show_unknown"))

    # annotation
    host = "claude" if blank(C.get("annotate_host")) else str(C["annotate_host"]).strip().lower()
    if host not in ("claude", "codex"):
        raise PlanError(f"annotate_host must be claude or codex, got {C['annotate_host']!r}")
    sec("annotate")["host"] = host
    put("annotate", "disease", "disease")
    put("annotate", "note", "annotate_note")
    replicates = 3 if blank(C.get("replicates")) else to_num(C["replicates"], "replicates", int)
    if replicates < 1:
        raise PlanError(f"replicates must be 1 or more, got {C['replicates']!r}")
    if replicates % 2 == 0:
        notes.append(f"replicates {replicates}: a majority needs an odd number, using {replicates + 1}")
        replicates += 1
    sec("annotate")["replicates"] = replicates
    # reference mode (labels limited to given cell types): an annotated h5ad, or a list of cell types
    ref_line = None
    if not blank(C.get("reference_h5ad")) and not blank(C.get("reference_celltypes")):
        raise PlanError("give the reference either as reference_h5ad (+ reference_column) or as "
                        "reference_celltypes, not both")
    if not blank(C.get("reference_h5ad")):
        ref = input_path(C["reference_h5ad"])
        if not Path(ref).exists():
            raise PlanError(f"reference_h5ad {ref!r} does not exist")
        sec("annotate")["reference_h5ad"] = ref
        put("annotate", "reference_column", "reference_column")
        ref_line = f"reference mode: {ref} : {C.get('reference_column')}"
    elif not blank(C.get("reference_celltypes")):
        v = C["reference_celltypes"]
        if isinstance(v, (list, tuple)):
            labels = [str(x).strip() for x in v if str(x).strip()]
        elif str(v).strip().lower().endswith(REF_LIST_SUFFIXES):
            f = input_path(v)
            if not Path(f).exists():
                raise PlanError(f"reference_celltypes {f!r} does not exist")
            labels = None
            sec("annotate")["reference_labels_file"] = f
            ref_line = f"reference mode: cell types listed in {f}"
        else:
            labels = [x.strip() for x in str(v).split(",") if x.strip()]
        if labels is not None:
            if not labels:
                raise PlanError("reference_celltypes is set but lists no cell type")
            sec("annotate")["reference_labels"] = labels
            ref_line = f"reference mode: {len(labels)} cell types ({', '.join(labels[:12])}{', ...' if len(labels) > 12 else ''})"
    if ref_line:
        if not blank(C.get("reference_scope")):
            sc = str(C["reference_scope"]).strip()
            if sc not in ("strict", "refine"):
                raise PlanError(f"reference_scope must be strict or refine, got {sc!r}")
            sec("annotate")["reference_scope"] = sc
        ref_line += f"; scope {A.get('annotate', {}).get('reference_scope', 'strict')}"
        if not blank(C.get("replicates")):
            notes.append("replicates is set but a reference run is one judgment (CeLLMan has no consensus for it); ignored")
    elif not blank(C.get("reference_column")) or not blank(C.get("reference_scope")):
        notes.append("reference_column / reference_scope are set but no reference is given; ignored")
    sec("annotate")["model"] = None if blank(C.get("llm_model")) else str(C["llm_model"]).strip()
    sb = "danger-full-access" if blank(C.get("codex_sandbox")) else str(C["codex_sandbox"]).strip()
    if sb not in ("workspace-write", "danger-full-access"):
        raise PlanError(f"codex_sandbox must be workspace-write or danger-full-access, got {sb!r}")
    sec("annotate")["codex_sandbox"] = sb

    # segmentation
    if not blank(C.get("segment")):
        sec("segment")["enabled"] = to_bool(C["segment"], "segment")
    if not blank(C.get("nuclei")):
        sec("segment")["nuclei"] = str(C["nuclei"]).strip()
    if not blank(C.get("nuclei_file")):
        sec("segment")["nuclei_file"] = input_path(C["nuclei_file"])
        sec("segment").setdefault("nuclei", "file")
        if blank(C.get("segment")):
            sec("segment")["enabled"] = True        # a nuclei file turns the segmentation on
    put("segment", "nuclei_frame", "nuclei_frame")
    put("segment", "seed_rule", "seed_rule")
    choices(C, "map_method", ("lda", "gnb", "corr"))
    choices(C, "kernel", KERNELS)
    choices(C, "vf_kernel", KERNELS)
    choices(C, "nuclei_frame", ("image", "um"))
    choices(C, "seed_rule", ("majority", "all"))

    # ---- build + validate (run_config.py)
    cfg = mc.build(insp, A)
    errs = mc.validate(cfg, insp)
    if kind == "zarr":
        # this check reads the raw folder layout; a zarr carries its images inside
        errs = [e for e in errs if not e.startswith("Visium HD: no full-resolution")]
    if tissue_missing:
        # asked for at the end, after the plan made from the data
        errs = [e for e in errs if e not in ("species is required", "organ is required")]
    if cosmx_z_missing:
        errs = [e for e in errs if not e.startswith("CosMx 3D needs cosmx_z_plane_um")]
    # the image downsampling uses: a 3D map on a 2D image has no z-stack step (SSAF-T places the image
    # at the middle of the transcripts' z range)
    if regime == "transcript":
        cfg["downsample_image"] = {"key": used if images else cfg["image"].get("key"),
                                   "kind": "zstack" if uses_stack else "2d"}
        if dim == 3 and not uses_stack:
            cfg["stages"] = [s for s in cfg["stages"] if s != "zstack"]
            errs = [e for e in errs if not e.startswith(STACK_ERRORS)]

    # nuclei: a SpatialData input names its own candidates; the user picks one, nothing is switched off
    if "segment" in cfg["stages"]:
        key = cfg["segment"].get("nuclei", "auto")
        if key == "auto":
            key = NUCLEI_DEFAULT.get(plat, "file")
        if kind == "zarr" and key != "file":
            tx = [k for k, v in z["points"].items() if v.get("gene_column")]
            cand_s = [k for k in z["elements"]["shapes"] if k not in OWN]
            cand_p = [k for k in z["elements"]["points"] if k not in OWN and k not in tx]
            if key not in cand_s + cand_p:
                errs.append(
                    f"segmentation needs nuclei, and {key!r} is not in this zarr. Set nuclei to one of its "
                    f"shapes {cand_s} or points {cand_p}, or nuclei: file with nuclei_file, or segment: false "
                    f"(map only)")
        if flat and key not in ("cells_polygons", "file"):
            errs.append(f"nuclei {key!r}: a flat SAW folder's cells are loaded as cells_polygons (from the "
                        f"cellbin.gef); leave nuclei blank, or nuclei: file with nuclei_file")
        cfg["segment"]["nuclei"] = key

    # ---- the files the load opens in a raw folder (reader_files.py): all there, or stop now
    files_line = None
    if kind == "raw":
        frows = RF.required(plat, p, cfg, zsel if uses_stack else None)
        gone = RF.missing_message(plat, p, frows, reader="stereoseq_flat" if flat else None)
        if gone:
            raise PlanError(gone)
        files_line = f"{len(frows)} files the load opens, all present: " + ", ".join(n for n, _, _ in frows)
    if errs:
        raise PlanError("the configuration cannot run:\n  - " + "\n  - ".join(errs))

    # ---- what the pipeline adds to run.json (keys run_config.build does not set)
    cfg["input_kind"] = kind
    cfg["gpu_request"] = gpu_req
    cfg["memory_gate"] = True if blank(C.get("memory_gate")) else to_bool(C["memory_gate"], "memory_gate")
    if kind == "zarr":
        cfg["zarr"] = {"elements": z["elements"], "zstack_element": zstack_element,
                       "cellman_table": next((k for k, t in z["tables"].items() if t["cellman"]), None),
                       "per_fov_cosmx": plat == "cosmx" and any(re.fullmatch(r"\d+_points", k)
                                                                for k in z["elements"]["points"])}
    if plat == "merscope" and kind == "raw":
        # element names after the data folder (ssaft_images.rename_merscope, right after loading)
        cfg["merscope_images"] = {"reader_id": IM.merscope_reader_id(insp["path"]), "name": mname,
                                  "mosaic_key": f"{mname}_z{mask_z}"}

    # ---- z planes (a Xenium / MERSCOPE z-stack): molecules per plane, a recommended range, the selection
    zp, sel, every = None, None, None
    if uses_stack and plat in ZSTACK_PLATFORMS:
        sel = zsel
        try:
            zp = ZP.census(kind, plat, insp, cfg, p, Path(run_dir))
        except Exception as e:  # noqa: BLE001
            notes.append(f"the molecules per z plane could not be counted ({type(e).__name__}: {e})")
        try:
            have = [r["plane"] for r in zp["planes"]] if zp else stack_planes(kind, cfg, p)
        except Exception as e:  # noqa: BLE001
            have = None
            if sel:
                raise PlanError(f"z_planes is set, but the planes of the z-stack could not be read: "
                                f"{type(e).__name__}: {e}")
        if sel:
            problem = ZP.check(sel, have)
            if problem:
                raise PlanError(problem)
            cfg["zstack"]["planes"] = sel
        elif have is not None and ZP.complete(have):
            # blank: every plane, allowed because the stack has every plane from 0. MERSCOPE's z-stack is
            # a glob of per-plane files: the plane numbers are passed as z=, so their order does not matter
            every = list(have)
            if plat == "merscope":
                cfg["zstack"]["planes"] = every
        elif have is not None:
            # asked for after the plan, which shows the molecules per plane and the recommended range
            pending.append(f"z_planes: the z-stack has only planes {have}, not every plane from 0, so z_planes "
                           f"cannot be left blank: the planes would be placed one z_step apart from 0 instead of "
                           f"at their own numbers, away from their molecules. Show the user the molecules per "
                           f"plane and the recommended range (the z planes table above) and ask which planes to "
                           f"use; set z_planes to a continuous range (e.g. 4-6) or a range at a constant step "
                           f"(e.g. 0-3:3), or image_key to the 2D image.")
        elif plat == "merscope":
            pending.append("z_planes: the planes of the MERSCOPE z-stack could not be read, so every plane cannot "
                           "be taken for a blank z_planes; set z_planes to the planes to use (e.g. 0-6).")
    elif uses_stack and zsel:
        cfg["zstack"]["planes"] = zsel

    # ---- the images block: candidates, the recommendation and the one used
    image_block = None
    if images:
        if uses_stack and (sel or every):
            for r in images["rows"]:
                if r["dims"] == 3 and r["name"] == used and r.get("shape") and len(r["shape"]) == 4:
                    r["shape"][1] = len(sel or every)
        image_block = IM.lines(images["rows"], used, dim, images["recommended"], images["why"], images["how"])

    # ---- summary
    summary.append(f"input      : {p} ({'SpatialData zarr' if kind == 'zarr' else 'raw output folder'}"
                   f"{'; a SAW output in one folder' if flat else ''})")
    summary.append(f"platform   : {plat}{' (experimental)' if P.get('experimental') else ''} -- {insp['detected_by']}"
                   + ("" if override else f"; candidates {cands}: confirm with the user, or set platform"))
    summary.append(f"loader     : {loader}")
    if files_line:
        summary.append(f"files      : {files_line}")
    if regime == "count":
        dline = " (count platform: 2D only)"
    elif images:
        dline = (f" (dimension{' blank: 3' if blank(C.get('dimension')) else ''}); downsampling on {used}, "
                 f"{'a z-stack' if uses_stack else 'a 2D image'}")
    else:
        dline = ""
    summary.append(f"map        : possible {dims}D; this run a {cfg['dimension']}D map{dline}")
    if kind == "zarr":
        for k in ("images", "points", "shapes", "tables", "labels"):
            if not z["elements"][k] or (k == "images" and image_block):
                continue
            detail = []
            for e in z["elements"][k]:
                if k == "images":
                    i = z["images"][e]
                    detail.append(f"{e} ({''.join(map(str, i['axes']))} {i['shape']})")
                elif k == "points":
                    detail.append(f"{e} (gene={z['points'][e]['gene_column']}, z={'yes' if z['points'][e]['has_z'] else 'no'})")
                elif k == "tables":
                    detail.append(f"{e} ({z['tables'][e]['shape']}{', uns[cellman]' if z['tables'][e]['cellman'] else ''})")
                else:
                    detail.append(e)
            summary.append(f"  {k:<8} : " + "; ".join(detail))
    else:
        for k, v in (insp.get("files") or {}).items():
            summary.append(f"  {k:<24}: {v}")
        for k, v in (insp.get("hints") or {}).items():
            if v in (None, [], "", {}):
                continue
            if plat == "merscope" and k in ("z_plane_um", "dataset_id_guess", "slide_name"):
                continue        # inspection guesses; the spacing and the element names used are below
            if plat == "cosmx" and k in ("px_um_candidates", "note"):
                continue        # the pixel size is on the `pixel size` line
            summary.append(f"  hint {k:<19}: {v}")
        if mname:
            summary.append(f"  element names           : {mname}_z{mask_z}, {mname}_transcripts "
                           f"(after the data folder {mname!r})")
    if px_line:
        summary.append(f"pixel size : {px_line}")
    summary.append(f"compute    : gpu={gpu_req}; GPUs visible in this container: {gpus_visible() or 'none'}; "
                   f"cores per stage {ncores} ({'threads' if not blank(C.get('threads')) else 'the cores the run is given'})")
    if image_block:
        summary.extend(image_block)
    if mask_line:
        summary.append(f"mask image : {mask_line}")
    if plat == "visiumhd":
        summary.append(f"bin size   : {bin_um:g} µm (visiumhd_bin_um; choices {list(VISIUMHD_BINS)})")
    if flat:
        summary.append(f"bin size   : {flat_bin:g} µm: the load sums the DNBs of {Path(cfg['load']['gem'] or '<chip>.tissue.gem.gz').name} "
                       f"into bins of this size ({'vf_bin_um' if not blank(C.get('vf_bin_um')) else 'vf_bin_um blank: the GEM' + chr(39) + 's own'}; "
                       f"a multiple of 0.5 µm), each at its centre in the chip's µm frame")
    summary.append(f"grid       : sampling distance {sdist:g} µm "
                   f"({'blank: the minimum' if blank(C.get('sampling_distance')) else 'sampling_distance'}; "
                   f"minimum {sd_min:g} µm = {sd_how}; any value from the minimum up)"
                   + (f"; blur raster bin {cfg['vf']['bin_um']:g} µm" if regime == "count" else ""))
    if regime == "transcript":
        ns = cfg["downsample"]["n_samples"]
        summary.append(f"anchors    : {ns} per 10 µm block of each imaging plane "
                       f"({'n_samples' if not blank(C.get('n_samples')) else 'n_samples blank: 10 for a 2D map, 5 for a 3D map'})")
    else:
        ns = cfg["downsample"]["n_samples_per_block"]
        summary.append(f"anchors    : {ns} per 10 µm block "
                       f"({'n_samples_per_block' if not blank(C.get('n_samples_per_block')) else f'n_samples_per_block blank: {N_SAMPLES_PER_BLOCK}'})")
    if cfg["dimension"] != 3 and not blank(C.get("seed_rule")):
        notes.append("seed_rule applies to a 3D map only; ignored")
    if cfg["dimension"] == 3 and plat in ("merscope", "cosmx"):
        key_um = "merscope_z_plane_um" if plat == "merscope" else "cosmx_z_plane_um"
        if z_um_missing or cosmx_z_missing:
            summary.append(f"z spacing  : not set ({key_um} blank): the µm between the transcripts' z planes is not "
                           f"recorded in the data" + (", so the planes are shown by number" if zp else ""))
        else:
            summary.append(f"z spacing  : {cfg['load']['z_plane_um']:g} µm per transcript z plane ({z_um_from or key_um})"
                           + (f"; z-stack spacing {cfg['zstack']['z_step']:g} µm (z_step)"
                              if uses_stack and cfg['zstack'].get('z_step') else ""))
    if zp:
        summary.extend(ZP.text(zp, sel, um=not z_um_missing))
    elif uses_stack and sel:
        summary.append(f"z planes   : {ZP.words(sel)} (z_planes: \"{ZP.spec(sel)}\")")
    ann = cfg["annotate"]
    summary.append(f"tissue     : species {A['species'] or '(blank)'}, organ {A['organ'] or '(blank)'}, "
                   f"condition {cfg['condition']}")
    summary.append(f"annotation : {host}, model {ann['model'] or '(CLI default; the name used is recorded in 00_check/env.json)'}"
                   f"; {ref_line or 'no reference'}"
                   + ("" if ref_line or ann.get("replicates", 3) == 1
                      else f"; {ann.get('replicates', 3)} judgments merged (CeLLMan consensus)"))
    if "segment" in cfg["stages"]:
        summary.append(f"nuclei     : {cfg['segment']['nuclei']}"
                       f"{' -> ' + str(cfg['segment'].get('nuclei_file')) if cfg['segment']['nuclei'] == 'file' else ''}"
                       + (f" (the cells of {Path(cfg['load']['cellbin']).name}, loaded as cells_polygons)"
                          if flat and cfg["segment"]["nuclei"] == "cells_polygons" else ""))
    else:
        summary.append("nuclei     : segmentation off")
    for n in notes:
        summary.append(f"note       : {n}")
    summary.append("")
    table = mc.summary(cfg)
    if cfg["dimension"] == 3 and regime == "transcript" and images:
        if uses_stack:
            if kind == "zarr" and zstack_element:
                stack_name = zstack_element
            else:
                stack_name = (f"the z-stack {used}, built by Analysis.zstack_image from "
                              f"{cfg.get('zstack', {}).get('source') or 'the z-stack files'}")
            if sel or every:
                stack_name += f" (planes {ZP.spec(sel or every)})"
            table = label_row(table, "image", f"2D reference grid for tiling and pixel coordinates; downsample uses "
                                              f"{stack_name}")
        else:
            table = label_row(table, "image", "the 2D image downsample uses, in this 3D map placed at the middle of "
                                              "the transcripts' z range (no z-stack)")
    if px_used is not None:
        table = label_row(table, "load", f"this run uses {px_used:g} µm/px ({px_src})"
                          + ("; px_um is used only for an unknown camera" if px_src.startswith("from") else ""))
    if z_um_missing:
        table = label_row(table, "load", "z_plane_um is the default until merscope_z_plane_um is set")
    summary.append(table)
    requests = ([tissue_missing] if tissue_missing else []) + pending
    if requests:
        # everything the data shows first, then every value still needed
        raise PlanError("the plan made from the data, before the values it still needs are set:\n"
                        + "=" * 77 + "\n" + "\n".join(summary) + "\n" + "=" * 77 + "\nstill needed:\n"
                        + "\n".join(f"  {i}. {r}" for i, r in enumerate(requests, 1)))
    return {"cfg": cfg, "inspect": insp, "answers": A, "summary": "\n".join(summary),
            "zplanes": dict(zp, selected=sel or every) if zp else None,
            "images": images}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config-json", required=True)
    ap.add_argument("--run-dir", required=True)
    a = ap.parse_args()
    try:
        out = plan(json.loads(a.config_json), a.run_dir)
    except PlanError as e:
        print(f"SSAF-T configuration problem: {e}")
        sys.exit(2)
    print("=" * 30 + " SSAF-T run plan " + "=" * 30, file=sys.stderr)
    print(out["summary"], file=sys.stderr)
    print("=" * 77, file=sys.stderr, flush=True)
    print(json.dumps(out, default=str))


if __name__ == "__main__":
    main()
