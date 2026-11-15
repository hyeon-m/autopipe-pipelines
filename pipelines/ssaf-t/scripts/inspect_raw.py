"""What a raw dataset folder is and holds, without loading it (standard library only): one detector
per platform returns {platform, supported, problem, dimensions, zstack, files, hints}, or None when
the folder is not that platform's layout. plan.py runs them in order and reports the first hit.
"""
from __future__ import annotations

import csv
import json
import re
from pathlib import Path


def _size_gb(p: Path) -> float:
    tot = 0
    for f in p.rglob("*"):
        try:
            if f.is_file() and not f.is_symlink():
                tot += f.stat().st_size
        except OSError:
            pass
    return round(tot / 2**30, 2)


def _first(p: Path, *patterns):
    for pat in patterns:
        hit = sorted(p.glob(pat))
        if hit:
            return hit[0]
    return None


def xenium(p: Path) -> dict | None:
    if not (p / "experiment.xenium").exists():
        return None
    stack = p / "morphology.ome.tif"
    focus = [str(x.relative_to(p)) for x in sorted(p.glob("morphology_focus*"))]
    exp = {}
    try:
        exp = json.loads((p / "experiment.xenium").read_text())
    except Exception:  # noqa: BLE001
        pass
    return dict(platform="xenium", supported=True, dimensions=[2, 3],
                zstack=stack.exists(),
                files=dict(transcripts=str(_first(p, "transcripts.parquet", "transcripts.csv.gz") or ""),
                           morphology_focus=focus, zstack=str(stack) if stack.exists() else None,
                           nucleus_boundaries=(p / "nucleus_boundaries.parquet").exists()),
                hints=dict(panel=exp.get("panel_name"), region=exp.get("region_name"),
                           analysis_sw_version=exp.get("analysis_sw_version")))


def merscope(p: Path) -> dict | None:
    tx = _first(p, "detected_transcripts.csv", "detected_transcripts.csv.gz")
    if tx is None:
        return None
    imgs = sorted((p / "images").glob("mosaic_*_z*.tif")) if (p / "images").is_dir() else []
    stains, planes = set(), set()
    for f in imgs:
        m = re.match(r"mosaic_(.+)_z(\d+)\.tif$", f.name)
        if m:
            stains.add(m.group(1)); planes.add(int(m.group(2)))
    return dict(platform="merscope", supported=True, dimensions=[2, 3],
                zstack=len(planes) > 1,
                files=dict(transcripts=str(tx), images=[f.name for f in imgs][:40],
                           cell_boundaries=bool(_first(p, "cell_boundaries.parquet"))),
                hints=dict(stains=sorted(stains), planes=sorted(planes),
                           region_name=p.name, slide_name=p.parent.name,
                           dataset_id_guess=f"{p.parent.name}_{p.name}",
                           z_plane_um=1.5))


def cosmx(p: Path) -> dict | None:
    tx = _first(p, "*_tx_file.csv", "*_tx_file.csv.gz")
    fovf = _first(p, "*fov_positions_file.csv", "*fov_positions_file.csv.gz")
    if tx is None or fovf is None:
        return None
    dataset_id = tx.name.split("_tx_file")[0]
    n_fov, fovs = None, []
    try:
        with open(fovf, newline="") as fh:
            rows = list(csv.DictReader(fh))
        fovs = sorted({int(float(r.get("fov") or r.get("FOV"))) for r in rows})
        n_fov = len(fovs)
    except Exception:  # noqa: BLE001
        pass
    header = []
    try:
        with open(tx, newline="") as fh:
            header = next(csv.reader(fh))
    except Exception:  # noqa: BLE001
        pass
    return dict(platform="cosmx", supported=True, experimental=True, dimensions=[2, 3],
                zstack=False,
                files=dict(transcripts=str(tx), fov_positions=str(fovf),
                           cell_composite=(p / "CellComposite").is_dir(),
                           cell_labels=(p / "CellLabels").is_dir()),
                hints=dict(dataset_id=dataset_id, n_fov=n_fov, fovs=fovs[:50],
                           transcript_columns=header, px_um_candidates=[0.12028, 0.18]))


def visiumhd(p: Path) -> dict | None:
    b = p / "binned_outputs"
    if not b.is_dir():
        if (p / "outs" / "binned_outputs").is_dir():
            return dict(platform="visiumhd", supported=False,
                        problem=f"point the dataset path at {p / 'outs'}")
        return None
    bins = sorted(x.name for x in b.glob("square_*um"))
    fullres = sorted(str(x) for x in list(p.glob("*.btf")) + list(p.glob("*tissue_image*.tif*"))
                     + list(p.parent.glob("*.btf")))
    return dict(platform="visiumhd", supported="square_002um" in bins, dimensions=[2], zstack=False,
                problem=None if "square_002um" in bins else "no binned_outputs/square_002um",
                files=dict(bins=bins, fullres_image_candidates=fullres),
                hints=dict(fullres_image_file=fullres[0] if fullres else None,
                           nuclei="none shipped: a nucleus segmentation (e.g. StarDist/ENACT) is needed for segmentation"))


def stereoseq(p: Path) -> dict | None:
    legacy = (p / "03.register").is_dir() and (p / "04.tissuecut").is_dir()
    flat_gef = sorted(p.glob("*.tissue.gef")) + sorted(p.glob("*.cellbin.gef"))
    if not legacy and not flat_gef:
        return None
    if not legacy:
        return dict(platform="stereoseq", supported=False,
                    problem="SAW output with every file in one directory (*.tissue.gef, *_regist.tif). "
                            "spatialdata_io.stereoseq reads only the numbered SAW layout "
                            "(02.count/, 03.register/, 04.tissuecut/, 041.cellcut/, 051.cellcluster/).",
                    files=dict(gef=[x.name for x in flat_gef]))
    need = {"03.register": "*_regist.tif", "04.tissuecut": "*.gef", "041.cellcut": "*.cellbin.gef",
            "051.cellcluster": "*cell.cluster.h5ad"}
    missing = [f"{d}/{pat}" for d, pat in need.items() if not list((p / d).glob(pat))]
    regist = sorted(x.name for x in (p / "03.register").glob("*_regist.tif"))
    gef = sorted((p / "04.tissuecut").glob("*.gef"))
    return dict(platform="stereoseq", supported=not missing,
                problem=("missing " + ", ".join(missing)) if missing else None,
                dimensions=[2], zstack=False,
                files=dict(regist=regist, tissuecut=[x.name for x in gef]),
                hints=dict(dataset_id=(gef[0].name.split(".")[0] if gef else None),
                           ssdna=[r for r in regist if "HE_" not in r],
                           he=[r for r in regist if "HE_" in r]))
