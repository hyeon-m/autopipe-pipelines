#!/usr/bin/env python
"""Run one SSAF-T stage for every unit.

    ssaft_stage.py --config .ssaft/run.json --stage downsample [--cores 16]
    ssaft_stage.py --config .ssaft/run.json --stage final

The stages are ssaft_load.py (load) and ssaft_steps.py (the others); ssaft_core.run_one runs each
for one unit with its log and state.json record. This script adds around them:

  * output layout (layout.py): after a stage, the files it made are moved from .ssaft/units/<unit>/
    into the stage's numbered folder, and a relative link is left behind, so later stages and
    re-runs read and write them where the stage code expects. `final` moves the SpatialData into
    09_final/;
  * the device: resolved here at run time (GPU when one is attached and gpu is not cpu), so
    changing the GPU never invalidates finished stages.
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import ssaft_core as core  # noqa: E402
from layout import FINAL, STAGE_DIR, STORE  # noqa: E402
from ssaft_core import log  # noqa: E402
from ssaft_steps import stages  # noqa: E402


def resolve_device(cfg) -> str | None:
    req = cfg.get("gpu_request", "auto")
    if req == "cpu":
        return None
    try:
        import torch
        ok = torch.cuda.is_available()
        name = torch.cuda.get_device_name(0) if ok else None
    except Exception as e:  # noqa: BLE001
        ok, name = False, None
        log(f"torch CUDA check failed: {e}")
    if not ok:
        if req == "cuda":
            raise SystemExit("gpu: cuda was requested but no CUDA device is usable in this container "
                             "(no GPU attached, or the driver does not support CUDA 13)")
        log("no usable CUDA device: running on the CPU")
        return None
    log(f"device cuda:0 ({name})")
    return "cuda:0"


# ------------------------------------------------------------------------- output layout
def root(run) -> Path:
    return Path(run.cfg["run_dir"]).parent


def multi(run) -> bool:
    return run.units() not in ([], ["sample"])


def stage_dir(run, folder, unit) -> Path:
    d = root(run) / folder
    if unit not in ("sample", "*") and multi(run):
        d = d / unit
    d.mkdir(parents=True, exist_ok=True)
    return d


def remove(p: Path):
    if p.is_symlink() or p.is_file():
        p.unlink()
    elif p.is_dir():
        shutil.rmtree(p)


def relocate(src: Path, dst: Path):
    """Move src into the layout at dst and leave a relative link at src. Idempotent."""
    if src.is_symlink():
        if src.resolve() == dst.resolve():
            return
        target = src.resolve()
        src.unlink()
        if target.exists() and target != dst.resolve():
            dst.parent.mkdir(parents=True, exist_ok=True)
            remove(dst)
            shutil.move(str(target), str(dst))
    elif src.exists():
        dst.parent.mkdir(parents=True, exist_ok=True)
        remove(dst)
        shutil.move(str(src), str(dst))
    else:
        return
    src.symlink_to(os.path.relpath(dst, src.parent))


def drop_linked(p: Path):
    """A link the layout left: remove it and what it points to (the stage re-creates both)."""
    if p.is_symlink():
        target = p.resolve()
        p.unlink()
        if target.exists():
            remove(target)


def prepare(run, unit, stage):
    """Before a (re-)run: a stage rmtree()s its own partial outputs, which fails on a link."""
    units_dir = Path(run.cfg["run_dir"]) / "units"
    if stage == "load" and units_dir.is_dir():
        for z in units_dir.glob("*/spatialdata.zarr"):
            drop_linked(z)
    if stage == "tile" and unit != "*":
        drop_linked(run.udir(unit) / "crops")


def cluster_png(h5ad: Path, png: Path):
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    a = ad.read_h5ad(h5ad, backed="r")
    xy = np.asarray(a.obsm["spatial"])
    lab = a.obs["leiden"].astype("category")
    cmap = plt.get_cmap("tab20")
    fig, ax = plt.subplots(figsize=(10, 10))
    for i, c in enumerate(lab.cat.categories):
        m = (lab == c).to_numpy()
        ax.scatter(xy[m, 0], xy[m, 1], s=0.2, color=cmap(i % 20), label=str(c), linewidths=0)
    ax.set_aspect("equal")
    ax.invert_yaxis()
    ax.set_axis_off()
    ax.legend(markerscale=20, fontsize=6, ncol=2, loc="upper left", bbox_to_anchor=(1.01, 1), frameon=False)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)


SECTIONS = {"load": ("load", "image", "merscope_images"), "zstack": ("zstack",), "tile": ("tile",),
            "panel_vf": ("panel", "vf"), "downsample": ("downsample", "image", "downsample_image"),
            "cluster": ("cluster",), "map": ("map",), "filter": ("filter",), "segment": ("segment",)}


def publish(run, unit, stage, seconds):
    d = stage_dir(run, STAGE_DIR[stage], unit)
    if unit != "*":
        u = run.udir(unit)
        if stage == "tile":
            relocate(u / "crops", d / "crops")
        if stage in ("panel_vf", "downsample"):
            relocate(u / "store", root(run) / STORE / (unit if multi(run) else ""))
        if stage == "panel_vf":
            relocate(u / "panel.txt", d / "panel.txt")
        if stage == "downsample":
            relocate(u / "anchors.h5ad", d / "anchors.h5ad")
        if stage == "cluster":
            cls = run.cls_path(unit)
            relocate(cls, d / cls.name)
            try:
                cluster_png(d / cls.name, d / "clusters_spatial.png")
            except Exception as e:  # noqa: BLE001
                log(f"clusters_spatial.png not drawn: {type(e).__name__}: {e}")
        if stage == "segment":
            relocate(u / "ssaft_cells.h5ad", d / "ssaft_cells.h5ad")
    rec = run.state()["units"].get(unit, {}).get(stage, {})
    out = {"stage": stage, "unit": unit, "status": rec.get("status"), "seconds": round(seconds),
           "device": run.cfg["gpu_device"] or "cpu", "ncores": run.cfg["ncores"],
           "params": {s: run.cfg.get(s) for s in SECTIONS.get(stage, ())},
           "summary": rec.get("summary"), "finished": rec.get("finished")}
    (d / "summary.json").write_text(json.dumps(out, indent=1, default=str))


def final(run):
    """09_final: the SpatialData with every input and SSAF-T element, and what it holds."""
    units = run.units()
    out = {}
    for unit in units:
        d = stage_dir(run, FINAL, unit)
        dst = d / "spatialdata.zarr"
        relocate(run.udir(unit) / "spatialdata.zarr", dst)
        el = {k: sorted(x.name for x in (dst / k).iterdir() if x.is_dir()) if (dst / k).is_dir() else []
              for k in ("images", "labels", "points", "shapes", "tables")}
        out[unit] = {"zarr": str(dst), "elements": el,
                     "segment": "done" if "segment" in run.cfg["stages"] else "not run"}
        (d / "summary.json").write_text(json.dumps(dict(out[unit], unit=unit,
                                                        finished=datetime.now().isoformat(timespec="seconds")),
                                                   indent=1))
    print(json.dumps({"ok": True, "stage": "final", "units": out}, indent=1))


# ---------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--cores", type=int, default=None)
    a = ap.parse_args()
    run = core.Run(a.config)
    if a.stage == "final":
        final(run)
        return
    if a.cores:
        run.cfg["ncores"] = a.cores
    run.cfg["gpu_device"] = resolve_device(run.cfg)
    STAGES = stages()
    if a.stage not in STAGES:
        raise SystemExit(f"unknown stage {a.stage!r}")
    units = run.stage_units(a.stage)
    if not units:
        raise SystemExit("no units yet: the load stage has not run")
    for unit in units:
        prepare(run, unit, a.stage)
        t0 = time.time()
        if not core.run_one(run, unit, a.stage, STAGES[a.stage]):
            rec = run.state()["units"].get(unit, {}).get(a.stage, {})
            raise SystemExit(f"stage {a.stage} failed for unit {unit}: {rec.get('error')} "
                             f"(log: {rec.get('log')})")
        publish(run, unit, a.stage, time.time() - t0)
        with open(Path(run.cfg["run_dir"]) / "devices.tsv", "a") as fh:
            fh.write(f"{a.stage}\t{unit}\t{run.cfg['gpu_device'] or 'cpu'}\t{run.cfg['ncores']}\t{time.time() - t0:.0f}s\n")
    print(json.dumps({"ok": True, "stage": a.stage, "units": run.units(), "device": run.cfg["gpu_device"]}))


if __name__ == "__main__":
    main()
