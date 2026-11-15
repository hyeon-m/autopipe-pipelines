#!/usr/bin/env python
"""Final report: per unit the cell types on the map, the kept voxel
fraction and the cells segmented; where the outputs are (the numbered stage folders, layout.py);
total run time; a PNG of the map in 07_filter/.

    report.py --config .ssaft/run.json
writes report.md and report.json next to the stage folders.

The device line names the GPU both ways: inside the container the attached GPUs are numbered from 0
(nvidia-smi there says "GPU 0"), and the server GPU it is comes from the GPUs this run was given
(AutoPipe's run record .autopipe-run.json, `gpus`, written when the run starts; 00_check/env.json
keeps the one of the run that made it, used when there is no record). The GPU names are read again
here, since a run resumed on another GPU keeps the 00_check of the first start.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from layout import CHECK, FINAL, STAGE_DIR, STORE  # noqa: E402


def rgb_preview(zarr_path: Path, store: Path, png: Path, max_side=3000):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    import spatialdata as sd
    import zarr
    from matplotlib.patches import Patch

    sdata = sd.read_zarr(str(zarr_path), selection=("images",))
    if "ssaft_celltype_map_rgb" not in sdata.images:
        return None
    el = sdata.images["ssaft_celltype_map_rgb"]
    if hasattr(el, "children") or not hasattr(el, "dims"):
        scales = sorted(el.keys(), key=lambda s: int("".join(c for c in s if c.isdigit()) or 0))
        pick = scales[-1]
        for s in scales:
            node = el[s]
            var = node[list(node.data_vars)[0]] if hasattr(node, "data_vars") else node["image"]
            if max(var.sizes["y"], var.sizes["x"]) <= max_side:
                pick = s
                break
        node = el[pick]
        da = node[list(node.data_vars)[0]] if hasattr(node, "data_vars") else node["image"]
    else:
        da = el
    if "z" in da.dims:
        da = da.isel(z=da.sizes["z"] // 2)
    arr = np.asarray(da.transpose("y", "x", "c").values)[..., :3]
    if arr.dtype != np.uint8:
        arr = np.clip(arr * (255 if arr.max() <= 1 else 1), 0, 255).astype(np.uint8)
    g = zarr.open_group(str(store), mode="r")
    names = [str(x) for x in g["celltype_names"][:]] if "celltype_names" in g else []
    colors = list(g.attrs.get("cellman_celltype_colors", []))
    png.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(12, 12 * arr.shape[0] / max(arr.shape[1], 1) + 0.5))
    ax.imshow(arr, interpolation="nearest")
    ax.set_axis_off()
    if names and len(colors) == len(names):
        ax.legend(handles=[Patch(color=c, label=n) for n, c in zip(names, colors)], loc="upper left",
                  bbox_to_anchor=(1.01, 1), fontsize=7, frameon=False)
    fig.savefig(png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return str(png)


def read_json(p: Path):
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return {}


def server_gpus(root: Path) -> dict | None:
    """`gpus` of AutoPipe's run record in the output folder (or a batch's above it)."""
    for d in [root, *root.parents]:
        j = read_json(d / ".autopipe-run.json")
        if j:
            return {"gpus": j.get("gpus"), "gpu_flags": j.get("gpu_flags"), "record": str(d / ".autopipe-run.json")}
    return None


def attached_gpus() -> list[str]:
    """`nvidia-smi -L` of this container: the GPUs of this run, which a resumed run may have changed."""
    try:
        out = subprocess.run(["nvidia-smi", "-L"], capture_output=True, text=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return []
    return [l.strip() for l in out.stdout.splitlines() if l.startswith("GPU ")] if out.returncode == 0 else []


def device_text(env: dict) -> str:
    """'cuda:0 = server GPU 1, NVIDIA ... (container GPU 0: ...; the run was given gpus 1)'."""
    dev = env.get("device") or "cpu"
    smi = env.get("nvidia_smi") or []
    sg = (env.get("server_gpus") or {}).get("gpus")
    sgs = str(sg).strip().lower() if sg not in (None, "") else ""
    ids = [s.strip() for s in str(sg).split(",") if s.strip()] if sgs and sgs not in ("all", "none") else None
    gpus = []
    for i, line in enumerate(smi):
        m = re.match(r"GPU (\d+): (.*?)(?: \(UUID: [^)]*\))?$", line.strip())
        k, name = (int(m.group(1)), m.group(2)) if m else (i, line.strip())
        srv = ids[k] if ids and k < len(ids) else (str(k) if sgs == "all" else None)
        gpus.append((k, srv, name))
    if not gpus:
        return f"{dev} (no GPU attached)"
    given = f"the run was given gpus {sg}" if sgs else "the server GPU number was not recorded"
    others = gpus
    if dev.startswith("cuda"):
        k0 = int(dev.split(":")[-1]) if ":" in dev else 0
        g = next((x for x in gpus if x[0] == k0), gpus[0])
        head = (f"{dev} = server GPU {g[1]}, {g[2]}" if g[1] is not None else f"{dev}, {g[2]}")
        head += f" (container GPU {g[0]}: inside the container the attached GPUs are numbered from 0; {given})"
        others = [x for x in gpus if x is not g]
    else:
        head = f"{dev} ({given})"
    if others:
        head += "; also attached: " + ", ".join(
            f"container GPU {k}" + (f" = server GPU {s}" if s is not None else "") + f" ({n})" for k, s, n in others)
    return head


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()
    cfg = json.loads(a.config.read_text())
    run_dir = Path(cfg["run_dir"])
    root = run_dir.parent
    st = read_json(run_dir / "state.json") or {"units": {}}
    ann = read_json(root / STAGE_DIR["annotate"] / "cellman" / "annotation_summary.json")
    env = read_json(root / CHECK / "env.json")
    # this run's record first: env.json is made once and may be from an earlier run on another GPU
    env["server_gpus"] = server_gpus(root) or env.get("server_gpus")
    env["nvidia_smi"] = attached_gpus() or env.get("nvidia_smi")

    units = st.get("unit_list") or ["sample"]
    many = units != ["sample"]
    sub = lambda folder, u: (root / folder / u) if many else (root / folder)  # noqa: E731
    rel = lambda p: str(Path(p).relative_to(root)) if p and str(p).startswith(str(root)) else p  # noqa: E731
    rep = {"platform": cfg["platform"], "dimension": cfg["dimension"], "dataset": cfg["dataset_path"],
           "tissue": f"{cfg['species']}, {cfg['organ']}, {cfg.get('condition')}", "stages": cfg["stages"],
           "stage_folders": {s: STAGE_DIR[s] for s in cfg["stages"]},
           "annotation": {k: ann.get(k) for k in ("host", "cli_version", "model", "model_used", "models_called",
                                                  "seconds")},
           "environment": {k: env.get(k) for k in ("device", "nvidia_smi", "server_gpus", "cpu_cores", "avx512f",
                                                   "ssaft")},
           "units": {}}
    rep["environment"]["device_text"] = device_text(env)
    started, finished = [], []
    used = ann.get("model_used")
    lines = ["# SSAF-T run report", "", f"- dataset: `{cfg['dataset_path']}`",
             f"- platform: {cfg['platform']}, {cfg['dimension']}D; tissue: {rep['tissue']}",
             f"- stages: {' -> '.join(cfg['stages'])}",
             f"- annotation: {ann.get('host')} {ann.get('cli_version') or ''}, model {ann.get('model')}"
             + (f" (used: {used})" if used and used != ann.get("model") else ""),
             f"- device: {rep['environment']['device_text']}", ""]
    for u in units:
        # a stage run once for every unit (e.g. the CosMx load) is recorded under the unit '*'
        rec = {**st["units"].get("*", {}), **st["units"].get(u, {})}
        for s in rec.values():
            if s.get("started"):
                started.append(s["started"])
            if s.get("finished"):
                finished.append(s["finished"])
        summ = {s: rec.get(s, {}).get("summary") for s in cfg["stages"]}
        labels = (ann.get("labels") or {}).get(u, {})
        cells = (summ.get("segment") or {}).get("cells")
        kept = (summ.get("filter") or {}).get("kept_voxel_fraction")
        zarr_final = sub(FINAL, u) / "spatialdata.zarr"
        store = (root / STORE / u) if many else (root / STORE)
        png = sub(STAGE_DIR["filter"], u) / "celltype_map.png"
        try:
            png = rgb_preview(zarr_final, store, png)
        except Exception as e:  # noqa: BLE001
            png = f"(preview failed: {type(e).__name__}: {e})"
        consensus = sub(STAGE_DIR["annotate"], u) / "annotation_out" / "consensus_report.md"
        outputs = {"spatialdata": rel(zarr_final), "store": rel(store),
                   "consensus_report": rel(consensus) if consensus.exists() else None,
                   "annotation_out": rel(sub(STAGE_DIR["annotate"], u) / "annotation_out"),
                   "anchors": rel(sub(STAGE_DIR["downsample"], u) / "anchors.h5ad"),
                   "clusters": rel(sub(STAGE_DIR["cluster"], u) / f"{cfg['species']}_{cfg['organ']}_cls.h5ad"),
                   "cells_h5ad": rel(sub(STAGE_DIR["segment"], u) / "ssaft_cells.h5ad") if cells is not None else None,
                   "celltype_map_png": rel(png)}
        rep["units"][u] = {"celltypes_on_map": (summ.get("map") or {}).get("celltypes"),
                           "unknown_groups_on_map": (summ.get("map") or {}).get("unknown_groups"),
                           "annotation_labels": labels, "kept_voxel_fraction": kept, "cells": cells,
                           "stage_seconds": {s: rec.get(s, {}).get("seconds") for s in cfg["stages"]},
                           "stage_summaries": summ, "outputs": outputs}
        lines += [f"## {u}" if many else "## Results", "",
                  f"- cell types on the map: {rep['units'][u]['celltypes_on_map']}"
                  f" (plus {rep['units'][u]['unknown_groups_on_map']} unknown groups)",
                  f"- kept voxel fraction: {kept}",
                  f"- cells segmented: {cells if cells is not None else 'segmentation not run'}",
                  "- annotation labels (anchors): " + ", ".join(f"{k} {v}" for k, v in labels.items()), "",
                  "| stage | folder | seconds | summary |", "|---|---|---|---|"]
        for s in cfg["stages"]:
            lines.append(f"| {s} | `{STAGE_DIR[s]}` | {rec.get(s, {}).get('seconds')} | "
                         f"{json.dumps(summ.get(s), default=str)[:300]} |")
        lines += ["", "Outputs:"] + [f"- {k}: `{v}`" for k, v in outputs.items() if v] + [""]
    if started and finished:
        t = (datetime.fromisoformat(max(finished)) - datetime.fromisoformat(min(started))).total_seconds()
        rep["wall_time_hours"] = round(t / 3600, 2)
        lines.insert(6, f"- wall time: {t / 3600:.2f} h")
    lines += [f"`{FINAL}/spatialdata.zarr` holds labels['ssaft_celltype_map'] and images['ssaft_celltype_map_rgb'] "
              "(+ shapes['ssaft_cell_boundaries_2d'], tables['ssaft_cells'] with segmentation), in the same "
              "coordinate system as the platform's own image. View with spatialdata_plot "
              "`render_images('ssaft_celltype_map_rgb')`, or the cells h5ad in scanpy."]
    (root / "report.json").write_text(json.dumps(rep, indent=1, default=str))
    (root / "report.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
