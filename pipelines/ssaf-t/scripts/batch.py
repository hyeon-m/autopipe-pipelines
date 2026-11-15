#!/usr/bin/env python
"""Run a batch of samples one after another.

    batch.py --plan batch_plan.json --snakefile /pipeline/Snakefile --cores 16

Each sample runs this pipeline's single-sample workflow, from load to the report, in its own
directory runs/<name>/ (the same numbered stage folders as a single run), with its config
(config.yaml + the sample's own values) written to runs/<name>/config.yaml:

  * samples run strictly one at a time, in manifest order; the next starts when one ends;
  * a sample whose plan did not validate is 'skipped' (batch_plan.py reported errors);
  * a failed sample is recorded and the batch moves on (batch_on_fail: stop ends the batch there);
  * a sample that failed in an earlier run is skipped unless batch_retry_failed is true;
  * a finished sample is run again only as far as a changed value requires (each sample is its own
    Snakemake run; with no change nothing re-runs);
  * batch_prepare_only stops every sample before the annotation ('prepared'); a later run without
    it continues them from the annotation;
  * a STOP file in the output directory stops the batch before the next sample (it is removed on
    exit), and so does a disk with less than batch_min_free_gb free.

batch_state.json holds per sample: status (pending / running / done / prepared / failed / skipped),
the current step, what it is waiting for (`waiting`: GPU or host memory, set only while it waits),
start, end, hours, error. batch_report.md sums the batch up. The exit code is
non-zero when the batch stopped early (STOP, low disk, batch_on_fail: stop), so running it again
resumes where it stopped.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from layout import CHECK, STAGE_DIR, WORK  # noqa: E402

ROOT = Path.cwd()
STATE = ROOT / "batch_state.json"
STOPF = ROOT / "STOP"


def now():
    return datetime.now().isoformat(timespec="seconds")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


def read_state():
    return json.loads(STATE.read_text()) if STATE.exists() else {"samples": {}}


WAITS = re.compile(r"^\[\d\d:\d\d:\d\d\] .* (waits for \d+ GiB free on|needs ~\d+ GiB; .*: waiting)")
ROOM = re.compile(r"^\[\d\d:\d\d:\d\d\] .*(has room after|memory free after) \d+ min")


def set_state(name, **kw):
    st = read_state()
    st["samples"].setdefault(name, {}).update(kw)
    st["updated"] = now()
    tmp = STATE.with_suffix(".tmp")
    tmp.write_text(json.dumps(st, indent=1, default=str))
    tmp.replace(STATE)


def run_sample(s, snakefile, cores, prepare_only) -> tuple[str, str | None]:
    rd = ROOT / "runs" / s["name"]
    rd.mkdir(parents=True, exist_ok=True)
    cfgf = rd / "config.yaml"
    cfgf.write_text(json.dumps(s["config"], indent=1, default=str) + "\n")     # JSON is YAML
    # a killed run leaves the sample's lock, which would fail the sample; samples run one at a time
    # and this batch holds the output folder, so a lock found here is stale
    for lock in (rd / ".snakemake" / "locks").glob("*.lock"):
        lock.unlink()
    cmd = [shutil.which("snakemake") or "snakemake", "--snakefile", snakefile, "--directory", str(rd),
           "--configfile", str(cfgf), "--cores", str(cores), "--rerun-incomplete"]
    if prepare_only:
        # up to the clustering, with the AI CLI check the annotation will need
        cmd += [f"{WORK}/stages/{STAGE_DIR['cluster']}.done", f"{CHECK}/env.json"]
    env = dict(os.environ, SSAFT_BATCH_CHILD="1")
    tail = []
    with open(rd / "pipeline.log", "a") as fh:
        fh.write(f"\n===== {now()} {' '.join(cmd)}\n")
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env, bufsize=1)
        for line in p.stdout:
            fh.write(line)
            fh.flush()
            tail = (tail + [line.rstrip()])[-40:]
            m = re.search(r"wildcards: step=(\S+)", line) or re.search(r"rule (annotate|final|report|check_env):", line)
            if m:
                step = m.group(1)
                set_state(s["name"], step=step, waiting=None)
                print(f"[{s['name']}] {step}", flush=True)
            elif WAITS.search(line):                    # a stage waiting for room is not stuck: say so
                set_state(s["name"], waiting=line.strip())
                print(f"[{s['name']}] {line.rstrip()}", flush=True)
            elif ROOM.search(line):
                set_state(s["name"], waiting=None)
                print(f"[{s['name']}] {line.rstrip()}", flush=True)
            elif line.startswith(("[annotate]", "Error", "SSAF-T configuration problem")) or "Finished jobid" in line \
                    or re.match(r"\[\d\d:\d\d:\d\d\] == ", line):
                print(f"[{s['name']}] {line.rstrip()}", flush=True)
        rc = p.wait()
    if rc != 0:
        return "failed", "\n".join(tail[-25:])
    return ("prepared" if prepare_only else "done"), None


def report(plan):
    st = read_state()["samples"]
    lines = ["# SSAF-T batch report", "", f"- samples: {len(plan['samples'])} in this run, {len({s['name'] for s in plan['samples']} | set(st))} in the batch", f"- updated: {now()}", "",
             "| # | sample | status | platform | tissue | hours | cell types | kept voxels | cells | report |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    out = []
    # a sample left out of this run (skip: true) keeps its row: its results are still in runs/<name>/
    left_out = [{"name": n, "order": "-"} for n in st if n not in {s["name"] for s in plan["samples"]}]
    for s in plan["samples"] + left_out:
        r = st.get(s["name"], {})
        rep = {}
        rj = ROOT / "runs" / s["name"] / "report.json"
        if rj.exists():
            try:
                rep = json.loads(rj.read_text())
            except ValueError:
                rep = {}
        units = rep.get("units", {})
        u = next(iter(units.values()), {}) if units else {}
        row = {"name": s["name"], "status": r.get("status", "pending"),
               "platform": s.get("platform", rep.get("platform")),
               "dimension": s.get("dimension", rep.get("dimension")), "tissue": s.get("tissue", rep.get("tissue")),
               "hours": r.get("hours"),
               "celltypes": u.get("celltypes_on_map"), "kept_voxel_fraction": u.get("kept_voxel_fraction"),
               "cells": u.get("cells"), "error": r.get("error"),
               "report": f"runs/{s['name']}/report.md" if rj.exists() else None}
        out.append(row)
        lines.append(f"| {s['order']} | {s['name']} | {row['status']} | {row['platform'] or ''} "
                     f"{str(row['dimension']) + 'D' if row['dimension'] else ''} | {row['tissue'] or ''} | "
                     f"{row['hours'] if row['hours'] is not None else ''} | {row['celltypes'] or ''} | "
                     f"{row['kept_voxel_fraction'] if row['kept_voxel_fraction'] is not None else ''} | "
                     f"{row['cells'] if row['cells'] is not None else ''} | "
                     f"{'`' + row['report'] + '`' if row['report'] else ''} |")
    errs = [r for r in out if r["error"]]
    if errs:
        lines += ["", "## Problems", ""]
        for r in errs:
            lines += [f"### {r['name']} ({r['status']})", "", "```", str(r["error"])[-3000:], "```", ""]
    lines += ["", "Each sample's results are in `runs/<sample>/` (00_check ... 09_final, report.md)."]
    (ROOT / "batch_report.md").write_text("\n".join(lines) + "\n")
    (ROOT / "batch_report.json").write_text(json.dumps(out, indent=1, default=str))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", required=True, type=Path)
    ap.add_argument("--snakefile", required=True)
    ap.add_argument("--cores", type=int, default=16)
    a = ap.parse_args()
    plan = json.loads(a.plan.read_text())
    B = plan["batch"]
    samples = plan["samples"]
    log(f"batch: {len(samples)} samples" + (" (prepare only: stopping each before the annotation)"
                                              if B["prepare_only"] else ""))
    stopped = None
    try:
        for s in samples:
            name = s["name"]
            prev = read_state()["samples"].get(name, {}).get("status")
            if not s["ok"]:
                set_state(name, status="skipped", error="plan did not validate: " + " | ".join(s["errors"]),
                          step=None)
                log(f"[{s['order']}/{len(samples)}] {name}: skipped (plan did not validate)")
                continue
            if prev == "failed" and not B["retry_failed"]:
                log(f"[{s['order']}/{len(samples)}] {name}: failed earlier; skipped (batch_retry_failed: false)")
                continue
            if prev == "prepared" and B["prepare_only"]:
                log(f"[{s['order']}/{len(samples)}] {name}: already prepared")
                continue
            if STOPF.exists():
                stopped = f"STOP file found; stopping before {name}"
                break
            free = shutil.disk_usage(ROOT).free / 2**30
            if B["min_free_gb"] and free < B["min_free_gb"]:
                stopped = f"only {free:.0f} GB free (< batch_min_free_gb {B['min_free_gb']:.0f}); stopping before {name}"
                break
            log(f"[{s['order']}/{len(samples)}] {name} ({s['platform']} {s['dimension']}D, {s['tissue']}, "
                f"overrides {json.dumps(s['overrides'], default=str)})")
            t0 = time.time()
            before = dict(read_state()["samples"].get(name, {}))
            set_state(name, status="running", started=now(), finished=None, error=None, step="check_env",
                      run_dir=f"runs/{name}")
            result, err = run_sample(s, a.snakefile, a.cores, B["prepare_only"])
            if prev == "done" and result == "done" and time.time() - t0 < 60:
                # a resume found nothing to do for this sample: its times stay those of the real run
                set_state(name, status=result, started=before.get("started"), finished=before.get("finished"),
                          hours=before.get("hours"), step=None, error=None)
            else:
                set_state(name, status=result, finished=now(), hours=round((time.time() - t0) / 3600, 2),
                          step=None, error=err)
            log(f"{name}: {result} in {(time.time() - t0) / 60:.0f} min" + (f"\n{err}" if err else ""))
            report(plan)
            if result == "failed" and B["on_fail"] == "stop":
                stopped = f"{name} failed and batch_on_fail is stop"
                break
    finally:
        STOPF.unlink(missing_ok=True)
        report(plan)
    counts = {}
    for s in samples:
        k = read_state()["samples"].get(s["name"], {}).get("status", "pending")
        counts[k] = counts.get(k, 0) + 1
    if stopped:
        log(f"batch stopped: {stopped}. {counts}. Run again to resume.")
        sys.exit(3)
    log(f"batch finished: {counts}")


if __name__ == "__main__":
    main()
