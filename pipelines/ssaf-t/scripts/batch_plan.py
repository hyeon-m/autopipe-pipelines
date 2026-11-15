#!/usr/bin/env python
"""Plan a batch of samples for this pipeline's config.yaml.

    batch_plan.py --config-json '<config as JSON>'

The shared values are config.yaml itself. The samples are either
  * `samples:` in config.yaml, a list whose items are a path or a mapping of `name`, `input` and any
    config key that differs for that sample (e.g. `nuclei`, `segment`, `organ`), or
  * `manifest:` a TSV file (columns order, name, path, platform, species,
    organ, dimension, condition, overrides, check); `path` is the sample's input, `overrides` a JSON
    object of config keys, and any further column named after a config key overrides that key. Rows
    are run in `order`; order `-` (or a blank order) leaves a row out.
Every sample is planned on its own with plan.py (the same inspection and validation as a single
run). A sample that does not validate is reported and skipped by the batch runner; the others run. When no sample validates there is nothing to run, and the
plan fails.

stdout: JSON {samples: [...], batch: {...}}; stderr: one line per sample (ok / FAIL) and each plan.
exit 2: no samples, no sample that validates, or a malformed sample list.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import plan as P  # noqa: E402

BATCH_KEYS = {"samples", "manifest", "batch_on_fail", "batch_retry_failed", "batch_prepare_only",
              "batch_min_free_gb"}
RUN_WIDE = {"threads", "gpu", "cores", "_cores"}   # one container runs every sample: these cannot differ
NOT_KEYS = {"name", "order", "check", "path", "overrides"}


class BatchError(Exception):
    pass


def unique_name(base: str, used: set) -> str:
    """A run folder name from base (letters, digits, _ . -), made unique among `used` with _2, _3 ..."""
    n = re.sub(r"[^A-Za-z0-9_.-]+", "_", base).strip("_")[:80] or "dataset"
    out, k = n, 2
    while out in used:
        out, k = f"{n}_{k}", k + 1
    used.add(out)
    return out


def rows_from_manifest(path: str) -> list[dict]:
    p = Path(P.input_path(path))
    if not p.exists():
        raise BatchError(f"manifest {str(p)!r} does not exist")
    with open(p, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    rows = [r for r in rows if (r.get("order") or "").strip() not in ("", "-")]
    rows.sort(key=lambda r: float(r["order"]))
    out = []
    for r in rows:
        s = {k: v for k, v in r.items() if k and v not in (None, "") and k not in ("order", "check", "overrides")}
        if r.get("path"):
            s["input"] = r["path"]
        s.pop("path", None)
        if (r.get("overrides") or "").strip():
            try:
                ov = json.loads(r["overrides"])
            except ValueError as e:
                raise BatchError(f"manifest row {r.get('name') or r.get('path')}: overrides is not JSON ({e})")
            if not isinstance(ov, dict):
                raise BatchError(f"manifest row {r.get('name')}: overrides must be a JSON object")
            s.update(ov)
        out.append(s)
    return out


def rows_from_config(samples) -> list[dict]:
    if not isinstance(samples, list):
        raise BatchError("samples must be a list (one item per sample: a path, or a mapping with input and the "
                         "keys that differ)")
    out = []
    for i, s in enumerate(samples, 1):
        if isinstance(s, str):
            s = {"input": s}
        if not isinstance(s, dict):
            raise BatchError(f"samples item {i} must be a path or a mapping, got {s!r}")
        if s.get("skip") in (True, "true", "yes"):
            continue
        out.append({k: v for k, v in s.items() if k != "skip"})
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config-json", required=True)
    a = ap.parse_args()
    C = json.loads(a.config_json)
    try:
        if C.get("samples") and not P.blank(C.get("manifest")):
            raise BatchError("give the samples either as samples: or as manifest:, not both")
        rows = rows_from_config(C["samples"]) if C.get("samples") else rows_from_manifest(C["manifest"])
        if not rows:
            raise BatchError("the batch has no samples")
    except BatchError as e:
        print(f"SSAF-T batch problem: {e}")
        sys.exit(2)

    base = {k: v for k, v in C.items() if k not in BATCH_KEYS}
    allowed = set(base) - RUN_WIDE
    used, samples = set(), []
    for i, r in enumerate(rows, 1):
        inp = r.get("input", base.get("input"))
        stem = Path(str(inp or f"sample{i}")).name
        stem = stem[:-5] if stem.endswith(".zarr") else stem
        name = unique_name(str(r.get("name") or (stem if stem not in ("outs", "output", "input", "") else f"sample{i}")), used)
        over = {k: v for k, v in r.items() if k not in NOT_KEYS}
        bad = sorted(set(over) - allowed)
        merged = dict(base, **over)
        rec = {"name": name, "order": i, "input": inp, "overrides": over, "config": merged}
        if bad:
            rec.update(ok=False, errors=[f"{', '.join(bad)}: not per-sample config keys "
                                         f"(threads and gpu are shared by the batch; valid keys: "
                                         f"{', '.join(sorted(allowed))})"])
        else:
            try:
                pl = P.plan(merged, ".ssaft")
                cfg = pl["cfg"]
                rec.update(ok=True, errors=[], platform=cfg["platform"], dimension=cfg["dimension"],
                           tissue=f"{cfg['species']}/{cfg['organ']}", stages=cfg["stages"],
                           segment=cfg["segment"].get("nuclei") if "segment" in cfg["stages"] else "off",
                           summary=pl["summary"])
            except P.PlanError as e:
                rec.update(ok=False, errors=[str(e)])
        samples.append(rec)

    try:
        batch = {"on_fail": "continue" if P.blank(C.get("batch_on_fail")) else str(C["batch_on_fail"]).strip(),
                 "retry_failed": False if P.blank(C.get("batch_retry_failed"))
                 else P.to_bool(C["batch_retry_failed"], "batch_retry_failed"),
                 "prepare_only": False if P.blank(C.get("batch_prepare_only"))
                 else P.to_bool(C["batch_prepare_only"], "batch_prepare_only"),
                 "min_free_gb": 200.0 if P.blank(C.get("batch_min_free_gb"))
                 else P.to_num(C["batch_min_free_gb"], "batch_min_free_gb"),
                 "threads": C.get("threads"), "gpu": C.get("gpu", "auto")}
    except P.PlanError as e:
        print(f"SSAF-T batch problem: {e}")
        sys.exit(2)
    if batch["on_fail"] not in ("continue", "stop"):
        print(f"SSAF-T batch problem: batch_on_fail must be continue or stop, got {batch['on_fail']!r}")
        sys.exit(2)

    ok = sum(s["ok"] for s in samples)
    print("=" * 30 + " SSAF-T batch plan " + "=" * 28, file=sys.stderr)
    print(f"{len(samples)} samples, {ok} ready; run one after another in this order; "
          f"on failure: {batch['on_fail']}; retry failed: {batch['retry_failed']}; "
          f"prepare only: {batch['prepare_only']}; stop below {batch['min_free_gb']:.0f} GB free", file=sys.stderr)
    for s in samples:
        if s["ok"]:
            print(f"ok   [{s['order']}] {s['name']}: {s['platform']} {s['dimension']}D {s['tissue']}, "
                  f"segment {s['segment']}, overrides {json.dumps(s['overrides'], default=str)}", file=sys.stderr)
        else:
            print(f"FAIL [{s['order']}] {s['name']}: " + " | ".join(s["errors"]), file=sys.stderr)
    for s in samples:
        if s["ok"]:
            print(f"\n----- {s['name']} -----\n{s['summary']}", file=sys.stderr)
    print("=" * 77, file=sys.stderr, flush=True)
    if ok == 0:
        print("SSAF-T batch problem: no sample of the batch can run, so there is nothing to do:\n  - "
              + "\n  - ".join(f"{s['name']}: " + " | ".join(s["errors"]) for s in samples))
        sys.exit(2)
    for s in samples:
        s.pop("summary", None)
    print(json.dumps({"samples": samples, "batch": batch}, default=str))


if __name__ == "__main__":
    main()
