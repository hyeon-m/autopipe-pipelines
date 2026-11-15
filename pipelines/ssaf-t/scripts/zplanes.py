"""z planes of a 3D run (Xenium, MERSCOPE): how many molecules each imaging plane of the z-stack
holds, which range of planes to keep, and the selection config.yaml makes (z_planes).

Used by plan.py (the census and the recommendation are printed at dry_run, and kept in
00_check/zstack_planes.json).

A plane is named by its own number: the page of a z-stack file (Xenium morphology.ome.tif), the
_z<N> of a per-plane file (MERSCOPE images/mosaic_<stain>_z<N>.tif), or the position along z of a 3D
image element of a SpatialData zarr. The kept planes must be evenly spaced (a range, or every k-th
plane of a range), because the z-stack element places its planes with one z scale and translation.

z_planes may stay blank (every plane) only when the z-stack has every plane from 0 on (0, 1, ...,
n-1), on Xenium and MERSCOPE alike: ssaft's zstack_image given no z= places the files' planes one
z_step apart in list order, which is right only then (a blank MERSCOPE z_planes is passed on as every
plane number, so the glob's order does not matter). A bundle shipping some planes only (z0 and z3, or
z4 to z6) would have them placed at 0, 1, 2 ... z_step instead of their own numbers, away from their
molecules; plan.py stops such a run and asks for z_planes, which passes the numbers as z= so each
plane sits at its own number.

The census of a file z-stack is ssaft's Analysis.zstack_planes, called with the plane numbers as z=
so that the z shown is the plane's own. A z-stack that is an image element of the input zarr has no
file for it, so its planes are taken from the element (their z through its 'global' transformation,
as ssaft reads them) and counted the same way: molecules whose z is a plane number are matched to
that plane exactly; molecules in µm are counted within one plane spacing centred on the plane.
Reading the molecules' z column is the slow part (a whole MERSCOPE CSV), so the census is cached in
the run directory, keyed by the files it was made from.

Before the MERSCOPE plane spacing is known (merscope_z_plane_um blank) the table is shown by plane
number only (text(um=False)): MERSCOPE molecules carry the plane number as their z, so the counts do
not depend on the spacing.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

RECOMMEND_SHARE = 0.95          # the range recommended holds at least this share of the molecules on planes
UM = ("µm", "um", "\u03bcm")
CENSUS_VERSION = 3              # part of the cache key: a census cached by an earlier version is made again


def header_spacing(images_dir, stain):
    """(µm between planes, how it was read) from the OME headers of the MERSCOPE mosaics of a stain, or
    (None, why not). The planes' own PositionZ decide when at least two planes carry one; otherwise
    PhysicalSizeZ. Planes that are not evenly spaced give None: one spacing cannot place them."""
    import tifffile
    pos, size = {}, set()
    for f in Path(images_dir).glob(f"mosaic_{stain}_z*.tif*"):
        m = re.search(r"_z(\d+)\.", f.name)
        if not m:
            continue
        try:
            with tifffile.TiffFile(str(f)) as t:
                xml = t.ome_metadata
        except Exception:  # noqa: BLE001 -- an unreadable header is a missing one
            xml = None
        if not xml:
            continue
        z = re.search(r'PositionZ="([0-9.eE+-]+)"', xml)
        unit = re.search(r'PositionZUnit="([^"]+)"', xml)
        if z and (unit is None or unit.group(1) in UM):
            pos[int(m.group(1))] = float(z.group(1))
        s = re.search(r'PhysicalSizeZ="([0-9.eE+-]+)"', xml)
        unit = re.search(r'PhysicalSizeZUnit="([^"]+)"', xml)
        if s and (unit is None or unit.group(1) in UM):
            size.add(float(s.group(1)))
    if len(pos) >= 2:
        planes = sorted(pos)
        steps = [(pos[b] - pos[a]) / (b - a) for a, b in zip(planes, planes[1:])]
        if max(steps) - min(steps) > 1e-3:
            return None, f"the image headers place the planes unevenly ({steps})"
        return round(abs(steps[0]), 4), f"PositionZ of {len(pos)} planes in the OME headers"
    if len(size) == 1:
        return size.pop(), "PhysicalSizeZ in the OME headers"
    return None, "the mosaics of this stain have no z in their headers"


class SelectionError(ValueError):
    pass


# ---------------------------------------------------------------------------------- the selection
def even(nums) -> bool:
    """Evenly spaced plane numbers (one or two planes always are)."""
    nums = list(nums)
    return len(nums) < 3 or len({b - a for a, b in zip(nums, nums[1:])}) == 1


def complete(nums) -> bool:
    """Every plane from 0 on (0, 1, ..., n-1): the only stack z_planes may leave blank for."""
    nums = list(nums)
    return nums == list(range(len(nums)))


def parse(text) -> list[int] | None:
    """config z_planes -> sorted plane numbers, or None (blank: every plane).

    "3-8" (planes 3 to 8), "2-10:2" (every second plane from 2 to 10), "2,4,6" or "5"."""
    if text is None or (isinstance(text, str) and not text.strip()):
        return None
    if isinstance(text, (list, tuple)):
        vals = [int(v) for v in text]
    else:
        s = str(text).strip().replace(" ", "")
        m = re.fullmatch(r"(\d+)-(\d+)(?::(\d+))?", s)
        if m:
            a, b, k = int(m.group(1)), int(m.group(2)), int(m.group(3) or 1)
            if b < a or k < 1:
                raise SelectionError(f"z_planes {text!r}: write a range low-high, e.g. 3-8, or low-high:k")
            vals = list(range(a, b + 1, k))
        elif re.fullmatch(r"\d+(,\d+)*", s):
            vals = [int(v) for v in s.split(",")]
        else:
            raise SelectionError(f"z_planes {text!r}: write a range (3-8), every k-th plane of a range "
                                 f"(2-10:2), or a list of evenly spaced planes (2,4,6)")
    vals = sorted(set(vals))
    if not even(vals):
        raise SelectionError(f"z_planes {text!r} picks planes {vals}, which are not evenly spaced. Only a "
                             f"continuous range (e.g. 3-8) or a range at a constant step (e.g. 2-10:2) can "
                             f"be used: the z-stack places its planes at one constant spacing")
    return vals


def check(sel: list[int], have: list[int]) -> str | None:
    """Why the selection cannot be used with these planes, or None."""
    missing = [v for v in sel if v not in have]
    if missing:
        return f"z_planes: no plane numbered {missing}; this z-stack has planes {have}"
    return None


def spec(sel: list[int]) -> str:
    """The z_planes value that writes this selection."""
    if len(sel) == 1:
        return str(sel[0])
    k = sel[1] - sel[0]
    return f"{sel[0]}-{sel[-1]}" + (f":{k}" if k != 1 else "")


def words(sel: list[int]) -> str:
    """planes 4 to 6 / planes 0 to 6, every 3rd / plane 5"""
    if len(sel) == 1:
        return f"plane {sel[0]}"
    k = sel[1] - sel[0]
    nth = {2: "2nd", 3: "3rd"}.get(k, f"{k}th")
    return f"planes {sel[0]} to {sel[-1]}" + (f", every {nth}" if k != 1 else "")


# ------------------------------------------------------------------------------- reading molecules
def _stat(p) -> list:
    try:
        st = Path(p).stat()
        return [str(p), st.st_size, st.st_mtime_ns]
    except OSError:
        return [str(p), None, None]


def _read_column(path: Path, names: list[str]):
    """The first of these columns of a parquet file / dataset or a (gzipped) CSV, as float64."""
    import numpy as np
    import pyarrow as pa
    if path.is_dir() or path.suffix == ".parquet":
        import pyarrow.dataset as ds
        d = ds.dataset(str(path), format="parquet")
        col = next((c for c in names if c in d.schema.names), None)
        if col is None:
            raise ValueError(f"{path}: none of the columns {names} (has {d.schema.names})")
        return np.asarray(d.to_table(columns=[col]).column(col).to_numpy(), dtype=np.float64), col
    import pyarrow.csv as pacsv
    with pacsv.open_csv(str(path)) as r:
        header = r.schema.names
    col = next((c for c in names if c in header), None)
    if col is None:
        raise ValueError(f"{path}: none of the columns {names} (has {header[:20]})")
    t = pacsv.read_csv(str(path), read_options=pacsv.ReadOptions(use_threads=True, block_size=1 << 26),
                       convert_options=pacsv.ConvertOptions(include_columns=[col],
                                                            column_types={col: pa.float64()}))
    return np.asarray(t.column(col).to_numpy(), dtype=np.float64), col


def molecules_file(kind: str, plat: str, insp: dict, zarr_path: Path):
    """(file or folder holding the molecules' z, column names to try)."""
    if kind == "raw":
        tx = (insp.get("files") or {}).get("transcripts")
        if not tx:
            raise ValueError("no transcripts file found in the dataset folder")
        return Path(tx), (["z_location", "z"] if plat == "xenium" else ["global_z", "z"])
    pts = (insp.get("zarr") or {}).get("points") or {}
    name = next((k for k, v in pts.items() if v.get("has_z") and v.get("gene_column")), None)
    if name is None:
        raise ValueError("the zarr has no points element with a z column")
    return zarr_path / "points" / name / "points.parquet", ["z"]


# --------------------------------------------------------------------------------------- planes
def file_planes(source, z_step):
    """[(number, z)] and the spacing of a z-stack file / glob, from its metadata and file names;
    a plane's z is its number times the spacing."""
    import glob
    from ssaft import Analysis
    src = sorted(glob.glob(source)) if "*" in str(source) else source
    planes, step = Analysis._image_planes(src, z_step=z_step, z_origin=0.0)
    return [(int(p["index"]), float(p["z"])) for p in planes], float(step), src


def element_planes(zarr_path: Path, key: str):
    """[(position, z)] and the spacing of a 3D image element of the input zarr."""
    import spatialdata as sd
    from ssaft import Analysis
    sdata = sd.read_zarr(str(zarr_path), selection=("images",))
    planes = Analysis._sdata_image_planes(sdata, key)
    zs = [float(p["z"]) for p in planes]
    step = (zs[1] - zs[0]) if len(zs) > 1 else 1.0
    return [(int(p["index"]), z) for p, z in zip(planes, zs)], float(step)


def _count(planes, step, tz, z_unit=None):
    """Molecules per plane, the way Analysis.zstack_planes counts them (for element planes)."""
    import numpy as np
    if z_unit is None:
        uniq = np.unique(tz) if tz.size else tz
        z_unit = ("plane" if (uniq.size and uniq.size <= 100 and np.all(uniq == np.rint(uniq))
                              and uniq.min() >= 0) else "um")
    half = float(step) / 2.0
    rows, inside = [], np.zeros(len(tz), dtype=bool)
    for k, zc in planes:
        m = (tz == k) if z_unit == "plane" else ((tz >= zc - half) & (tz < zc + half))
        inside |= m
        rows.append([k, zc, int(m.sum())])
    return rows, z_unit, int((~inside).sum())


def census(kind, plat, insp, cfg, zarr_path: Path, cache_dir: Path) -> dict:
    """{planes: [{plane, z, molecules, fraction, cumulative}], unit, outside, total, step, source,
    recommended: {planes, spec, share_planes, share_all}} for the run's z-stack."""
    Zc = cfg["zstack"]
    element = (cfg.get("zarr") or {}).get("zstack_element") if kind == "zarr" else None
    tx, names = molecules_file(kind, plat, insp, zarr_path)
    if element:
        src_desc = f"images[{element!r}] of {zarr_path}"
        key = ["element", str(zarr_path), element, _stat(zarr_path / "images" / element / "zarr.json"),
               _stat(zarr_path / "images" / element / ".zattrs"), _stat(tx)]
    else:
        src_desc = str(Zc["source"])
        key = ["file", src_desc, Zc.get("z_step"), _stat(tx)]
    key = json.dumps([CENSUS_VERSION] + key, sort_keys=True, default=str)
    cache_p = cache_dir / "zstack_census.json"
    try:
        cache = json.loads(cache_p.read_text())
    except (OSError, ValueError):
        cache = {}
    if key in cache:
        return cache[key]

    tz, col = _read_column(tx, names)
    if element:
        planes, step = element_planes(zarr_path, element)
        rows, unit, outside = _count(planes, step, tz)
        how = src_desc
    else:
        from ssaft import Analysis
        planes, step, src = file_planes(Zc["source"], Zc.get("z_step"))
        nums = [k for k, _ in planes]
        # the plane numbers as z=, so each plane is shown (and counted) at its own number
        zsel = nums if len(nums) > 1 and even(nums) else None
        df = Analysis.zstack_planes(src, z=zsel, transcripts_z=tz, z_step=Zc.get("z_step"), verbose=False)
        rows = [[int(k), float(r["z"]), int(r["molecules"])] for k, r in df.iterrows()]
        unit, outside = df.attrs["unit"], int(df.attrs["outside"])
        how = f"{src_desc} (Analysis.zstack_planes)"
    total = int(len(tz))
    out, cum = [], 0.0
    for k, zc, n in rows:
        frac = n / max(total, 1)
        cum += frac
        out.append({"plane": k, "z": zc, "molecules": n, "fraction": frac, "cumulative": cum})
    res = {"planes": out, "unit": unit, "outside": outside, "total": total, "step": float(step),
           "source": how, "molecules_from": f"{tx} [{col}]", "recommended": recommend(out, total)}
    cache = {k: v for k, v in cache.items() if k.startswith(f"[{CENSUS_VERSION},")}
    cache[key] = res
    try:
        cache_dir.mkdir(parents=True, exist_ok=True)
        tmp = cache_p.with_suffix(".tmp")
        tmp.write_text(json.dumps(cache))
        tmp.replace(cache_p)
    except OSError:
        pass
    return res


def recommend(rows: list[dict], total: int, share: float = RECOMMEND_SHARE) -> dict | None:
    """The shortest run of consecutive, evenly spaced planes holding at least `share` of the
    molecules on planes (the most molecules among equally short runs)."""
    on = sum(r["molecules"] for r in rows)
    if not rows or on == 0:
        return None
    best = None
    for i in range(len(rows)):
        for j in range(i, len(rows)):
            sel = [r["plane"] for r in rows[i:j + 1]]
            if not even(sel):
                break
            n = sum(r["molecules"] for r in rows[i:j + 1])
            if n >= share * on:
                cand = (len(sel), -n, sel, n)
                if best is None or cand[:2] < best[:2]:
                    best = cand
                break
    if best is None:
        return None
    _, _, sel, n = best
    return {"planes": sel, "spec": spec(sel), "share_planes": n / on, "share_all": n / max(total, 1)}


# ------------------------------------------------------------------------------------------ text
def text(res: dict, sel: list[int] | None, um: bool = True) -> list[str]:
    """The census as lines of the dry_run plan. um=False: by plane number only (the spacing in µm is
    not known yet)."""
    if res["unit"] == "plane":
        unit = "molecule z is the plane number: matched to each plane exactly"
    else:
        unit = f"molecule z is in µm: counted within {res['step']:g} µm around each plane"
    rows = res["planes"]
    spacing = (f"{res['step']:g} µm per plane number" if um else
               "by plane number (the µm between planes is not set yet)")
    lines = [f"z planes   : {len(rows)} planes of {res['source']}, {spacing}; "
             f"{res['total']:,} molecules ({res['molecules_from']}); {unit}",
             f"  {'':1} {'plane':>5}" + (f"  {'z (µm)':>8}" if um else "")
             + f"  {'molecules':>12}  {'fraction':>8}  {'cumulative':>10}"]
    for r in rows:
        mark = "*" if sel and r["plane"] in sel else " "
        lines.append(f"  {mark} {r['plane']:>5}" + (f"  {r['z']:>8.2f}" if um else "")
                     + f"  {r['molecules']:>12,}  {r['fraction']:>8.1%}  {r['cumulative']:>10.1%}")
    if res["outside"]:
        lines.append(f"  on no plane of the z-stack: {res['outside']:,} molecules "
                     f"({res['outside'] / max(res['total'], 1):.1%}; e.g. on planes the bundle has no image of)")
    rec = res.get("recommended")
    if rec:
        lines.append(f"  recommended: {words(rec['planes'])} (z_planes: \"{rec['spec']}\") hold "
                     f"{rec['share_planes']:.1%} of the molecules on the planes ({rec['share_all']:.1%} of all), "
                     f"the shortest range with at least {RECOMMEND_SHARE:.0%}")
    have = [r["plane"] for r in rows]
    if sel:
        n = sum(r["molecules"] for r in rows if r["plane"] in sel)
        lines.append(f"  this run   : {words(sel)} (z_planes: \"{spec(sel)}\"; {len(sel)} of {len(rows)}, "
                     f"marked *), holding {n / max(res['total'], 1):.1%} of all molecules")
    elif complete(have):
        lines.append(f"  this run   : all {len(rows)} planes (z_planes blank: the z-stack has every plane from 0)")
    else:
        lines.append(f"  this run   : none yet: the z-stack has only planes {have}, not every plane from 0, so "
                     f"z_planes has to name the planes to use")
    lines.append("  z_planes takes a continuous range (3-8) or a range at a constant step (2-10:2) only")
    return lines
