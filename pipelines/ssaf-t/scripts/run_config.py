"""run.json from the platform defaults (platform_defaults.json) and the answers plan.py derives from
config.yaml (standard library only).

  build(inspect, answers)   the run.json dict: every section's platform default, the answers merged over
  validate(cfg, inspect)    the combinations that would fail late (hours in), as error strings
  summary(cfg)              the stage plan and values as a table
  modules(server, species)  the Cell-Ontology modules the CeLLMan server serves

answers hold only what config.yaml sets; anything absent keeps its platform default. Keys: run_dir,
python, dimension (2|3), species, organ, condition, gpu_device, ncores and per-section overrides, each a
dict merged over the default section: load, zstack, tile, image, panel, vf, downsample, cluster,
annotate, map, filter, segment.
"""
from __future__ import annotations

import copy
import json
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULTS = json.loads((HERE / "platform_defaults.json").read_text())
SECTIONS = ("load", "zstack", "tile", "image", "panel", "vf", "downsample", "cluster", "annotate",
            "map", "filter", "segment")


def merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = merge(out[k], v)
        else:
            out[k] = v
    return out


def build(inspect: dict, answers: dict) -> dict:
    plat = inspect["platform"]
    if plat not in DEFAULTS["platforms"]:
        raise ValueError(f"unsupported platform {plat!r}")
    P = DEFAULTS["platforms"][plat]
    dim = int(answers.get("dimension", P["default_dimension"]))
    cfg = {
        "platform": plat,
        "regime": P["regime"],
        "dimension": dim,
        "dataset_path": inspect["path"],
        "run_dir": str(Path(answers["run_dir"]).expanduser().resolve()),
        "python": answers.get("python", ""),
        "species": answers.get("species", ""),
        "organ": answers.get("organ", ""),
        "condition": answers.get("condition", "normal"),
        "gpu_device": answers.get("gpu_device", "cuda:0"),
        "ncores": int(answers.get("ncores", DEFAULTS["common"]["ncores"])),
        "cellman_server": answers.get("cellman_server", DEFAULTS["cellman_server"]),
    }
    for s in SECTIONS:
        base = copy.deepcopy(DEFAULTS["common"].get(s, {}))
        if s == "tile" and P["regime"] == "transcript" and not P.get("per_fov"):
            base["enabled"] = False               # the whole section at once; tile only when asked
        pv = P.get(s, {})
        if s == "downsample":
            pv = pv.get(str(dim)) or pv.get(str(P["default_dimension"]), {})
        base = merge(base, {k: v for k, v in pv.items() if not (s == "tile" and k == "enabled" and not P.get("per_fov"))})
        cfg[s] = merge(base, answers.get(s, {}))
    # hints from the inspection fill 'auto' / null load fields
    h = inspect.get("hints", {})
    if plat in ("cosmx", "stereoseq") and not cfg["load"].get("dataset_id"):
        cfg["load"]["dataset_id"] = h.get("dataset_id")
    if plat == "visiumhd" and cfg["load"].get("fullres_image_file") == "auto":
        cfg["load"]["fullres_image_file"] = h.get("fullres_image_file")
    if plat == "visiumhd" and "bin_um" not in answers.get("vf", {}):
        # vf.bin_um is the raster the bins are blurred on, so it follows the bin level actually
        # loaded. Finer than the data leaves the raster mostly empty -- 8 µm bins on a 2 µm
        # raster fill 3.4% of its nodes, and a 2.5 µm sigma cannot bridge the 8 µm gaps, so the
        # field comes out as isolated blobs instead of a continuous one.
        cfg["vf"]["bin_um"] = float(cfg["load"].get("bin_size", 2))
    if plat == "xenium" and dim == 3 and not str(cfg["zstack"].get("source", "")).startswith("/"):
        cfg["zstack"]["source"] = str(Path(inspect["path"]) / cfg["zstack"]["source"])
    if plat == "merscope" and dim == 3 and not str(cfg["zstack"].get("source", "")).startswith("/"):
        cfg["zstack"]["source"] = str(Path(inspect["path"]) /
                                      cfg["zstack"]["source"].format(stain=cfg["zstack"].get("stain", "DAPI")))
    stages = list(P["stages"].get(str(dim)) or P["stages"][str(P["default_dimension"])])
    if "tile" in stages and not cfg["tile"].get("enabled"):
        stages.remove("tile")
    if "segment" in stages and not cfg["segment"].get("enabled", True):
        stages.remove("segment")
    cfg["stages"] = stages
    return cfg


def validate(cfg: dict, inspect: dict) -> list[str]:
    """VALIDATION: the combinations that would fail late (hours in) are refused here."""
    err = []
    P = DEFAULTS["platforms"][cfg["platform"]]
    if not inspect.get("supported", False):
        err.append(f"dataset not readable: {inspect.get('problem')}")
    if cfg["dimension"] not in P["dimensions"]:
        err.append(f"{cfg['platform']} runs in {P['dimensions']}D only, not {cfg['dimension']}D")
    if cfg["dimension"] == 3:
        src = cfg["zstack"].get("source")
        if not src:
            err.append("3D needs a z-stack source (zstack_source)")
        elif "*" not in str(src) and not Path(src).exists():
            err.append(f"z-stack source not found: {src}")
        if cfg["platform"] in ("merscope", "cosmx") and not cfg["zstack"].get("z_step"):
            err.append("3D needs the z plane spacing in µm (merscope_z_plane_um / cosmx_z_plane_um)")
        if cfg["platform"] == "cosmx" and not cfg["load"].get("z_plane_um"):
            err.append("CosMx 3D needs cosmx_z_plane_um to convert the transcripts' plane index to µm")
    if cfg["platform"] == "visiumhd" and not cfg["load"].get("fullres_image_file"):
        # the reader still gives the spatial/ hires image, whose scalefactor puts it in the same frame.
        # At 0.18-0.28 x full resolution that is ~1-1.5 µm/px, at or below the 2 µm bin lattice, so it
        # serves as the tissue mask and the portal's background; only a bundle with no image at all
        # cannot run.
        spatial = Path(inspect["path"]) / "binned_outputs" / f"square_{int(cfg['load'].get('bin_size', 2)):03d}um" / "spatial"
        if not any((spatial / f"tissue_{r}_image.png").exists() for r in ("hires", "lowres")):
            err.append("Visium HD: no full-resolution H&E and no spatial/tissue_*_image.png; "
                       "set visiumhd_fullres_image")
    if cfg["platform"] == "cosmx" and not cfg["load"].get("px_um"):
        err.append("CosMx needs cosmx_px_um (µm per pixel)")
    seg = cfg["segment"]
    if "segment" in cfg["stages"] and seg.get("nuclei") == "file" and not seg.get("nuclei_file"):
        err.append("segmentation on this platform needs a nuclei file (nuclei_file), or set segment: false")
    if seg.get("nuclei_file") and not Path(seg["nuclei_file"]).exists():
        err.append(f"nuclei file not found: {seg['nuclei_file']}")
    ref = cfg["annotate"].get("reference_h5ad")
    if ref and not cfg["annotate"].get("reference_column"):
        err.append("reference_h5ad given without reference_column")
    if not cfg["python"] or not Path(cfg["python"]).exists():
        err.append(f"python not found: {cfg['python']}")
    for k in ("species", "organ"):
        if not cfg.get(k):
            err.append(f"{k} is required")
    return err


def summary(cfg: dict) -> str:
    rows = [("platform", cfg["platform"] + (" (experimental)" if DEFAULTS["platforms"][cfg["platform"]].get("experimental") else "")),
            ("dimension", f"{cfg['dimension']}D"), ("dataset", cfg["dataset_path"]), ("run dir", cfg["run_dir"]),
            ("tissue", f"{cfg['species']}, {cfg['organ']}, {cfg['condition']}"),
            ("device", cfg["gpu_device"] or "CPU"), ("stages", " -> ".join(cfg["stages"]))]
    for s in SECTIONS:
        if s in ("zstack",) and cfg["dimension"] != 3:
            continue
        if s in ("panel", "vf") and cfg["regime"] != "count":
            continue
        if s in ("tile",) and "tile" not in cfg["stages"]:
            continue
        if s == "segment" and "segment" not in cfg["stages"]:
            rows.append(("segment", "off"))
            continue
        vals = ", ".join(f"{k}={v}" for k, v in cfg[s].items() if v is not None and not k.endswith("_unit"))
        rows.append((s, vals or "defaults"))
    w = max(len(r[0]) for r in rows)
    return "\n".join(f"{k.ljust(w)} | {v}" for k, v in rows)


def modules(server: str, species: str | None) -> list:
    with urllib.request.urlopen(server.rstrip("/") + "/modules", timeout=30) as r:
        data = json.loads(r.read())
    wings = data.get("wings", data if isinstance(data, list) else [])
    return [w for w in wings if species is None or w.get("species") == species]
