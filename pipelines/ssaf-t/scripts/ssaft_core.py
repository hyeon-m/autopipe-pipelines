"""The run record and what every stage shares: run.json, state.json, the unit layout, logging, the
SpatialData helpers, the SSAF-T Analysis, and the memory / GPU gates a stage waits behind.

A run lives in <output>/.ssaft/ (run.json "run_dir"):
  state.json                       status per unit and stage (running / done / failed), timings, a
                                   short summary, the log path
  units/<unit>/logs/<stage>.log    the stage's stdout/stderr (run_one tees into it)
  units/<unit>/spatialdata.zarr    the SpatialData (inputs + every SSAF-T element written back)
  units/<unit>/store/              the SSAF-T zarr store
A unit is one independent analysis: 'sample', or one per FOV for a per-FOV CosMx zarr ('fov<N>').
"""
from __future__ import annotations

import contextlib
import datetime as _dt
import fcntl
import json
import os
import re
import sys
import time
import traceback
from pathlib import Path

os.environ.setdefault("NUMPY_MADVISE_HUGEPAGE", "0")   # large numpy allocations can stall for minutes otherwise

ORDER = ["load", "zstack", "tile", "panel_vf", "downsample", "cluster", "annotate", "map", "filter", "segment"]


# ----------------------------------------------------------------------------------------------- state
class Run:
    def __init__(self, cfg_path: Path):
        self.cfg_path = cfg_path.resolve()
        self.cfg = json.loads(self.cfg_path.read_text())
        self.dir = Path(self.cfg["run_dir"])
        self.dir.mkdir(parents=True, exist_ok=True)
        self.state_path = self.dir / "state.json"

    # state.json is read by report.py and annotate.py too; lock around every read-modify-write
    @contextlib.contextmanager
    def _locked_state(self):
        lock = open(self.dir / ".state.lock", "w")
        fcntl.flock(lock, fcntl.LOCK_EX)
        try:
            st = json.loads(self.state_path.read_text()) if self.state_path.exists() else {"units": {}}
            yield st
            tmp = self.state_path.with_suffix(".tmp")
            tmp.write_text(json.dumps(st, indent=1))
            tmp.replace(self.state_path)
        finally:
            fcntl.flock(lock, fcntl.LOCK_UN)
            lock.close()

    def state(self) -> dict:
        return json.loads(self.state_path.read_text()) if self.state_path.exists() else {"units": {}}

    def units(self) -> list[str]:
        st = self.state()
        if st.get("unit_list"):
            return st["unit_list"]
        return ["sample"] if self.cfg["platform"] != "cosmx" else []

    def set_units(self, units):
        with self._locked_state() as st:
            st["unit_list"] = list(units)

    def mark(self, unit, stage, **kw):
        with self._locked_state() as st:
            rec = st["units"].setdefault(unit, {}).setdefault(stage, {})
            rec.update(kw)

    def stage_units(self, stage) -> list[str]:
        """A CosMx load decides the units itself (one slide, or one per FOV); that load is unit '*'."""
        if self.cfg["platform"] == "cosmx" and stage == "load":
            return ["*"]
        return self.units()

    def udir(self, unit) -> Path:
        if unit == "*":
            (self.dir / "logs").mkdir(parents=True, exist_ok=True)
            return self.dir
        d = self.dir / "units" / unit
        (d / "logs").mkdir(parents=True, exist_ok=True)
        return d

    def cls_path(self, unit) -> Path:
        return self.udir(unit) / f"{self.cfg['species']}_{self.cfg['organ']}_cls.h5ad"


class Tee:
    def __init__(self, path):
        self.fh = open(path, "a", buffering=1)
        self.out = sys.__stdout__

    def write(self, s):
        self.fh.write(s)
        self.out.write(s)

    def flush(self):
        self.fh.flush()
        self.out.flush()


def log(*a):
    print(f"[{_dt.datetime.now():%H:%M:%S}]", *a, flush=True)


# ------------------------------------------------------------------------------------------- helpers
def scale0(el):
    """Full-resolution DataArray of an image/labels element (single-scale DataArray or multiscale DataTree)."""
    if hasattr(el, "children") or not hasattr(el, "shape"):
        keys = list(el.keys())
        node = el["scale0"] if "scale0" in keys else el[keys[0]]
        names = list(node.data_vars) if hasattr(node, "data_vars") else list(node.keys())
        return node["image"] if "image" in names else node[names[0]]
    return el


def read_sdata(path):
    import spatialdata as sd
    return sd.read_zarr(str(path))


def persist(sdata, keys):
    """Write elements the analysis added to the in-memory SpatialData back into its zarr."""
    on_disk = set()
    try:
        on_disk = set(sdata.elements_paths_on_disk())
    except Exception:  # noqa: BLE001
        pass
    for k in keys:
        present = any(k in getattr(sdata, kind) for kind in ("images", "labels", "points", "shapes", "tables"))
        if not present:
            continue
        with contextlib.suppress(Exception):
            if any(p.endswith("/" + k) for p in on_disk) or k in on_disk:
                sdata.delete_element_from_disk(k)
        try:
            sdata.write_element(k, overwrite=True)
        except TypeError:
            sdata.write_element(k)
        log(f"wrote element {k!r} to {sdata.path}")


def device(run):
    return run.cfg.get("gpu_device") or None


def transcripts_key(sdata, platform):
    pts = list(sdata.points)
    if "transcripts" in pts:
        return "transcripts"
    cand = [k for k in pts if k.endswith("_transcripts")] or pts
    if not cand:
        raise ValueError("no transcripts points element in the SpatialData")
    return cand[0]


def image_key_2d(run, sdata):
    img = run.cfg["image"]
    key = img.get("key", "auto")
    if img.get("channel") is not None:
        return "ssaft_mask_image"
    if key and key != "auto":
        if key not in sdata.images:
            raise KeyError(f"image {key!r} not in {sorted(sdata.images)}")
        return key
    plat = run.cfg["platform"]
    if plat == "merscope":
        z = run.cfg["load"].get("z_layers", 3)
        z = z[0] if isinstance(z, list) else z
        hit = [k for k in sdata.images if k.endswith(f"_z{z}")]
        if hit:
            return hit[0]
    cands = [k for k in sdata.images if not k.startswith("ssaft_")]
    if not cands:
        raise KeyError("no image element")
    return cands[0]


# ----------------------------------------------------------------- SSAF-T 23dd889 workaround 2 (remove
# when SSAF-T fixes it): the tiled (chunked) 3D downsample reads each tile with selection=('images',)
# and then asks _transcripts_z_mid of that object, which has no transcripts; its log line and errors
# name a variable `images` that the function does not have (NameError under verbose).
def _zmid_with_points(an):
    orig = an._transcripts_z_mid

    def zmid(sdata, element_name=None):
        if not getattr(sdata, "points", None) and getattr(sdata, "path", None):
            import spatialdata as sd
            sdata = sd.read_zarr(str(sdata.path), selection=("points",))
        return orig(sdata, element_name)
    an._transcripts_z_mid = zmid


def analysis(run, unit):
    """ssaft.Analysis over the unit's store, its downsample_vector made safe for the tiled 3D path."""
    import ssaft
    an = ssaft.Analysis(ssaft.Dataset(str(run.udir(unit) / "store")), ncores=run.cfg["ncores"], verbose=True)
    _zmid_with_points(an)
    dv = an.downsample_vector

    def downsample_vector(*a, **k):
        import ssaft._analysis as _A
        had = hasattr(_A, "images")
        _A.images = k.get("image_key")                 # the name its messages use (workaround 2)
        try:
            return dv(*a, **k)
        finally:
            if not had:
                del _A.images
    an.downsample_vector = downsample_vector
    return an


# ------------------------------------------------------------------------------------- memory gate
# bytes of host memory per map voxel a stage needs at its peak (measured on 3D Xenium maps: filter
# peaked at ~34 B/voxel on 8.69 G voxels, segment at ~30 B/voxel on 8.0 G); 40 keeps a margin
MEM_PER_VOXEL = {"filter": 40, "segment": 40, "map": 30}   # map: 57-126 GiB on 2-8 G-voxel 3D maps
MEM_RESERVE_GIB = 80           # kept free for everything else on the host
# stages whose peak does not follow the map size (full sections): a MERSCOPE downsample holds the
# z-stack plane pipeline and ~0.5 billion transcripts (~220 GiB), a MERSCOPE map the transcripts plus
# ~10 M anchor vectors (143 GiB), a 45-FOV CosMx whole-slide load the mosaic and its pyramid
MEM_PLATFORM_GIB = {("merscope", "downsample"): 230, ("merscope", "map"): 160, ("cosmx", "load"): 200,
                    ("cosmx", "downsample"): 150}


def _meminfo_gib(key):
    for line in Path("/proc/meminfo").read_text().splitlines():
        if line.startswith(key + ":"):
            return int(line.split()[1]) / 2**20
    return float("inf")


_HELD = {}            # lock path -> open handle this process still holds after its stage started
_HELD_BY_STAGE = []   # gate locks the running stage owns; freed as soon as it finishes


def _hold(lockp, fh, seconds):
    """Keep a gate lock for `seconds` after a stage starts, so the next stage measures after it allocated."""
    import threading
    _HELD[str(lockp)] = fh
    t = threading.Timer(seconds, _drop, args=(lockp, fh))
    t.daemon = True
    t.start()


def _drop(lockp, only=None):
    """Release a gate lock this process holds. The next stage of the same process calls this first:
    otherwise it would queue behind its own lock until the timer fires."""
    if only is not None and _HELD.get(str(lockp)) is not only:
        return                                  # a later stage already re-took the lock: its own timer frees it
    fh = _HELD.pop(str(lockp), None)
    if fh is not None:
        try:
            fcntl.flock(fh, fcntl.LOCK_UN)
            fh.close()
        except (OSError, ValueError):
            pass


def estimated_need(run, unit, stage):
    """GiB of host memory the stage is expected to need, or None when it is not a heavy stage."""
    need = MEM_PLATFORM_GIB.get((run.cfg.get("platform"), stage))
    per = MEM_PER_VOXEL.get(stage)
    if need is not None or not per:
        return need
    import numpy as np
    shape = None
    try:
        import zarr
        g = zarr.open_group(str(run.udir(unit) / "store"), mode="r")
        shape = g["celltype_maps"].shape if "celltype_maps" in g else g["vf_norm"].shape if "vf_norm" in g else None
    except Exception:  # noqa: BLE001
        pass
    if shape is None:          # before the first map: the lattice downsample reported
        logp = run.udir(unit) / "logs" / "downsample.log"
        m = re.findall(r"vf grid \((\d+), (\d+), (\d+)\)", logp.read_text(errors="ignore")) if logp.exists() else []
        if m:
            shape = tuple(int(v) for v in m[-1])
    return None if shape is None else float(np.prod(shape)) * per / 2**30


def wait_for_memory(run, unit, stage):
    """Before a memory-heavy stage, wait until the host has room for it, so two heavy stages on one
    host (a batch, or another job) do not start together. The check-and-start is serialised by a
    lock next to the run directory, held a few minutes after the start so the stage has allocated
    before the next one measures. A need above what this host can ever free does not wait (it would
    wait forever); memory_gate: false turns the gate off."""
    if not run.cfg.get("memory_gate", True):
        return None
    need = estimated_need(run, unit, stage)
    if need is None:
        return None
    cap = _meminfo_gib("MemTotal") - MEM_RESERVE_GIB
    if need > cap:
        log(f"{stage} is estimated at ~{need:.0f} GiB (full-section measurement), more than this host can "
            f"free ({cap:.0f} GiB = total - {MEM_RESERVE_GIB} reserve): starting without waiting")
        return None
    lockp = Path(run.cfg["run_dir"]).parent / ".heavy_stage.lock"
    _drop(lockp)
    t0, noted = time.time(), False
    while True:
        fh = open(lockp, "a")
        fcntl.flock(fh, fcntl.LOCK_EX)
        if _meminfo_gib("MemAvailable") - MEM_RESERVE_GIB >= need:
            break
        fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()
        if not noted:
            log(f"{stage} needs ~{need:.0f} GiB; {_meminfo_gib('MemAvailable'):.0f} GiB available: waiting")
            noted = True
        time.sleep(60)
    if noted:
        log(f"memory free after {(time.time() - t0) / 60:.0f} min; starting {stage}")
    _hold(lockp, fh, 600)                     # a 3D segment takes ~7 min to reach its peak
    return need


# ---------------------------------------------------------------------------------------- GPU gate
# GPU-heavy stages start only when their card has this much free (GiB): a map started on a nearly
# full card dies with a CUDA out-of-memory an hour in
GPU_FREE_GIB = {"panel_vf": 10, "downsample": 10, "map": 12}
# MERSCOPE carries ~0.5 billion transcripts: downsample uploads them whole; map keeps them on the
# host and uploads one x-slab at a time
GPU_FREE_GIB_PLATFORM = {("merscope", "map"): 16, ("merscope", "downsample"): 26, ("cosmx", "map"): 18}


def wait_for_gpu(run, stage):
    need = GPU_FREE_GIB_PLATFORM.get((run.cfg.get("platform"), stage), GPU_FREE_GIB.get(stage))
    dev = device(run)
    if not need or not dev or not str(dev).startswith("cuda"):
        return
    try:
        import torch
        idx = int(str(dev).split(":")[-1]) if ":" in str(dev) else 0
        free_gib = lambda: torch.cuda.mem_get_info(idx)[0] / 2**30  # noqa: E731
        total_gib = torch.cuda.mem_get_info(idx)[1] / 2**30
    except Exception:  # noqa: BLE001
        return
    need = min(need, 0.9 * total_gib)          # a need above the card's size would wait forever
    lockp = Path(run.cfg["run_dir"]).parent / f".gpu{idx}_stage.lock"
    _drop(lockp)
    t0, noted = time.time(), False
    while True:
        # hold the lock only to check and start, never while waiting
        fh = open(lockp, "a")
        fcntl.flock(fh, fcntl.LOCK_EX)
        if free_gib() >= need:
            break
        fcntl.flock(fh, fcntl.LOCK_UN)
        fh.close()
        if not noted:
            log(f"{stage} waits for {need:.0f} GiB free on {dev} ({free_gib():.1f} GiB now)")
            noted = True
        # a waiter must not sit on the memory reservation it took a moment ago
        _drop(Path(run.cfg["run_dir"]).parent / ".heavy_stage.lock")
        time.sleep(60)
    if noted:
        log(f"{dev} has room after {(time.time() - t0) / 60:.0f} min; starting {stage}")
    # heavy GPU stages keep the card for as long as they run (run_one frees it); light ones only
    # need the allocation window
    _HELD_BY_STAGE.append(lockp)
    _hold(lockp, fh, 24 * 3600 if stage in ("map", "downsample") else 180)


# ------------------------------------------------------------------------------------------ run_one
def run_one(run, unit, stage, fn):
    """Run one stage function for one unit: wait behind the gates, tee its output into the unit log,
    and record running / done / failed with timings and the summary it returns in state.json."""
    logf = run.udir(unit) / "logs" / f"{stage}.log"
    wait_for_memory(run, unit, stage)
    wait_for_gpu(run, stage)
    t0 = time.time()
    run.mark(unit, stage, status="running", started=_dt.datetime.now().isoformat(timespec="seconds"),
             log=str(logf), pid=os.getpid(), error=None)
    tee = Tee(logf)
    ok = True
    with contextlib.redirect_stdout(tee), contextlib.redirect_stderr(tee):
        log(f"== {stage} [{unit}] ==")
        try:
            summary = fn(run, unit)
            log(f"== {stage} [{unit}] done in {time.time() - t0:.0f}s: {summary}")
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            summary, ok = None, False
            run.mark(unit, stage, status="failed", error=f"{type(e).__name__}: {e}"[:2000],
                     finished=_dt.datetime.now().isoformat(timespec="seconds"), seconds=round(time.time() - t0))
    while _HELD_BY_STAGE:                       # the stage is over: another stage may take the card now
        _drop(_HELD_BY_STAGE.pop())
    if ok:
        run.mark(unit, stage, status="done", summary=summary, seconds=round(time.time() - t0),
                 finished=_dt.datetime.now().isoformat(timespec="seconds"))
    return ok
