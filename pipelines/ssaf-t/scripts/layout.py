"""Output layout: one numbered folder per stage, shared by the Snakefile and the scripts.

    00_check/        env.json, gpu_info.txt, plan.txt, inspect.json, answers.json, run.json, images.json,
                     zstack_planes.json
    01_load/         summary.json
    02_zstack/       summary.json                       (a 3D map on a z-stack)
    02a_tile/        crops/, summary.json               (tiled Xenium / MERSCOPE)
    02_panel_vf/     panel.txt, summary.json            (Visium HD, Stereo-seq)
    03_downsample/   anchors.h5ad, summary.json
    04_cluster/      <species>_<organ>_cls.h5ad, clusters_spatial.png, summary.json
    05_annotation/   annotation_out/ (CeLLMan), cellman/ (prompt, session events), summary.json
    06_map/          summary.json
    07_filter/       celltype_map.png, summary.json
    08_segment/      ssaft_cells.h5ad, summary.json     (when segmentation runs)
    09_final/        spatialdata.zarr (every input and SSAF-T element), summary.json
    ssaft_store/     the SSAF-T zarr store (written by several stages)
    report.md, report.json
    .ssaft/          the run directory (run.json, state.json, units/ with links
                     to the files above, logs)

A CosMx per-FOV run has one subfolder per unit inside each stage folder.
"""
import json

WORK = ".ssaft"
CHECK = "00_check"
FINAL = "09_final"
STORE = "ssaft_store"
STAGE_DIR = {
    "load": "01_load",
    "zstack": "02_zstack",
    "tile": "02a_tile",
    "panel_vf": "02_panel_vf",
    "downsample": "03_downsample",
    "cluster": "04_cluster",
    "annotate": "05_annotation",
    "map": "06_map",
    "filter": "07_filter",
    "segment": "08_segment",
}
DIR_STAGE = {v: k for k, v in STAGE_DIR.items()}

# 05_annotation/cellman/<this>: what the annotation in annotation_out/ was made from
ANNOTATION_INPUTS = "annotation_inputs.json"


def annotation_inputs(cfg: dict, state: dict, units) -> dict:
    """What an annotation depends on: the run.json values the annotate stage reads (its params in
    rules/single.smk) and, per unit, the clustering its anchors come from (the cluster stage's record
    in state.json, which a re-clustering replaces)."""
    clusters = {}
    for u in units:
        rec = (state.get("units") or {}).get(u, {}).get("cluster", {})
        clusters[u] = {"finished": rec.get("finished"), "summary": rec.get("summary")}
    values = {k: cfg.get(k) for k in ("species", "organ", "condition", "annotate")}
    return json.loads(json.dumps({"values": values, "clusters": clusters}, sort_keys=True, default=str))
