#!/usr/bin/env python
"""The annotate stage: CeLLMan, run headless by Claude Code (claude -p) or Codex (codex exec).

    annotate.py --config .ssaft/run.json

The prompt (build_prompt) asks CeLLMan's annotate-mode to annotate every unit's clustered anchors with
one annotate-lead agent per unit, which merges `replicates` independent annotate-expert judgments.
A reference run (annotate-reference-mode) and `replicates: 1` use one annotate-expert per unit, because
CeLLMan's consensus does not take a reference. The prompt tells the session it runs headless and fully
autonomous, and that the run ends at the standard annotation: CeLLMan's interactive phase (the free
annotation after finalize) is not part of the pipeline, and outputs() refuses its files. Because a pipeline step has to finish on its own:

  * it waits for the session to end;
  * the session's working directory is the run directory (.ssaft), prepared with CeLLMan's skill and
    agent (.claude/ for Claude Code, .agents/ + .codex/ for Codex) and CeLLMan/ (the helper path the
    prompt names), so the agent finds them and Codex may write the annotation outputs;
  * the cellman-kg MCP server is registered for the session (--mcp-config / -c mcp_servers...);
  * Codex gets network access in its sandbox (CeLLMan reads the Cell-Ontology server);
  * a reference given as a list of cell types (config reference_celltypes) instead of an annotated
    h5ad: build_prompt knows only the h5ad form, so the prompt is switched to annotate-reference-mode
    here and given the `- reference : <label list>` line that skill reads (an inline list is written
    to .ssaft/reference_labels.json first).

Outputs go to 05_annotation/ (layout.py): cellman/ (prompt, session events, summary) and, once the
session has ended, annotation_out/ moved out of .ssaft/units/<unit>/ (a link is left there for the
map stage). The summary records the model asked for (llm_model, blank = the CLI's default) and the
model the session's events report (model_used; models_called also lists e.g. a subagent's model).
The stage fails, with the session's last message, when any unit lacks *_annotated.h5ad or
celltype_markers.yaml.

Re-runs: cellman/annotation_inputs.json records what a session was started for (the run.json values
the annotation reads and the clustering of every unit, layout.annotation_inputs). When the stage runs
again (e.g. the run was interrupted after the annotation, or during it once CeLLMan had written its
outputs) with the same inputs and every unit's outputs complete, those outputs are used and no new
session is started. With changed inputs, the previous annotation_out is moved aside
(annotation_out.prev-<time>) and CeLLMan runs again.
"""
from __future__ import annotations

import argparse
import fcntl
import json
import os
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

CELLMAN = Path(os.environ.get("CELLMAN_DIR", "/opt/CeLLMan"))
sys.path.insert(0, str(Path(__file__).parent))
from check_env import events, login, used_models  # noqa: E402
from layout import ANNOTATION_INPUTS as INPUTS, STAGE_DIR, annotation_inputs  # noqa: E402
from ssaft_stage import relocate, resolve_device  # noqa: E402

PY = sys.executable

AUTONOMOUS = """
## Where CeLLMan is

The helper is {helper}. Use that path, and give it to every agent you spawn. Do not
search for it: the data directories hold millions of files and a `find` over them does not finish.

## You are running headless and FULLY AUTONOMOUS

**Spawn every annotate-lead or annotate-expert agent with `run_in_background: false`.** Agents run in the background by
default, and a background agent does not hold your turn open: this session then ends with the agent still
working, the runtime kills it minutes later, and the dataset fails with nothing written. Your turn must not
end until the agent has returned its result — do not spawn and report "the agent is running", do not poll a
background task, and do not end your turn early. One turn: spawn (foreground), wait, report what it returned.
Tell every annotate-lead in its spawn prompt that the same holds for it: it spawns its judges with
`run_in_background: false`, has every judge's result before it goes on, and does not end its turn before
finalize has run.

Nobody reads this session and NOTHING you ask will ever be answered. Do not ask the user anything.
This is the default annotation, end to end: run it and finish.

Wherever annotate-mode would have you ask (module confirmation, a coarse or ambiguous tissue name, what
to do with the result after finalize, whether to explore something further), take the best-supported
option on the evidence you have, write one line in your final summary saying what you chose and why,
and carry on. Finalize the annotation as it stands: no optional extra exploration rounds, no recovery
passes, no relabelling beyond the method's own rules. STOP only if the input itself is unusable (a
missing h5ad, no module for the tissue) and report that as FAILED.

This run is the standard annotation only, finalized once into `cellman_annotation`. CeLLMan's interactive
phase is not part of it: never write `cellman_free_annotation` (no `--column`), do not annotate a second
time after finalize, and split clusters with `descend` only, not `subcluster`. Free-annotation outputs make
this pipeline stage fail.

When every dataset is finished, end with one line per dataset:
    DONE <input h5ad> -> <annotated h5ad>
or, for a dataset that stopped:
    FAILED <input h5ad>: <reason>
"""


def unit_dirs(run_dir: Path, cfg: dict) -> list[str]:
    """The run's units (state.json unit_list), 'sample' before load has written it."""
    st = run_dir / "state.json"
    units = json.loads(st.read_text()).get("unit_list") if st.exists() else None
    return units or ["sample"]


def build_prompt(cfg: dict, run_dir: Path, host: str, workdir: Path) -> str:
    """The caller prompt: annotate-mode (annotate-reference-mode with a reference h5ad) over every
    unit's clustered anchors, one agent per unit, and the AUTONOMOUS block. The agent is annotate-lead
    (the consensus of `replicates` judges), or annotate-expert for a reference or `replicates: 1`."""
    units = unit_dirs(run_dir, cfg)
    ref = cfg["annotate"].get("reference_h5ad")
    skill = "annotate-reference-mode" if ref else "annotate-mode"
    replicates = int(cfg["annotate"].get("replicates") or 3)
    listed = cfg["annotate"].get("reference_labels") or cfg["annotate"].get("reference_labels_file")
    consensus = replicates > 1 and not ref and not listed
    agent = "annotate-lead" if consensus else "annotate-expert"
    if host != "claude":
        agent = agent.replace("-", "_")
    if host == "claude":
        head = f"/{skill} "
    else:
        head = f"${skill} "
    tissue = f"{cfg['species']}, {cfg['organ']}" + (f", {cfg['condition']}" if cfg.get("condition") else "")
    python = cfg["python"]
    dev = cfg.get("gpu_device")
    note_bits = [f"input anchors are an SSAF-T spatial vector field ({cfg['platform']}, "
                 f"{cfg['dimension']}D), normalised and log1p-transformed"]
    if dev:
        note_bits.append(f"use GPU {dev}")
    if cfg["annotate"].get("note"):
        note_bits.append(cfg["annotate"]["note"])
    note = "; ".join(note_bits)

    lines = []
    if len(units) == 1:
        u = units[0]
        cls = run_dir / "units" / u / f"{cfg['species']}_{cfg['organ']}_cls.h5ad"
        lines += [f"{head}Annotate the h5ad file below using the {agent} agent.", "",
                  f"- input          : {cls}",
                  f"- tissue         : {tissue}",
                  f"- virtual env    : {Path(python).parent.parent}",
                  "- cluster column : leiden",
                  f"- out            : {run_dir / 'units' / u / 'annotation_out'}"]
    else:
        lines += [f"{head}Annotate the h5ad files below using the {agent} agent.",
                  "Run each dataset as an INDEPENDENT agent IN PARALLEL — one agent per dataset, each with its",
                  "OWN --out directory.",
                  "Do NOT annotate, read any out directory, or run helper.py yourself — only spawn the agents and",
                  "wait for their results; each out directory belongs to its agent alone.",
                  "Settle the KG module for EVERY dataset first, and only then start all datasets together.", "",
                  "Datasets:"]
        for i, u in enumerate(units, 1):
            cls = run_dir / "units" / u / f"{cfg['species']}_{cfg['organ']}_cls.h5ad"
            lines += [f"  [{i}] input  : {cls}",
                      f"      tissue : {tissue}",
                      f"      out    : {run_dir / 'units' / u / 'annotation_out'}"]
        lines += ["", "Common:",
                  f"- virtual env    : {Path(python).parent.parent}",
                  "- cluster column : leiden"]
    if consensus:
        lines.append(f"- replicates     : {replicates}")
    if cfg["annotate"].get("disease"):
        lines.append(f"- disease        : {cfg['annotate']['disease']}")
    if ref:
        lines += [f"- reference h5ad   : {ref}",
                  f"- reference column : {cfg['annotate']['reference_column']}",
                  f"- reference scope  : {cfg['annotate'].get('reference_scope', 'strict')}"]
    lines.append(f"- note           : {note}")
    helper = workdir / "CeLLMan" / "scripts" / "helper.py"
    return "\n".join(lines) + "\n" + AUTONOMOUS.format(helper=helper)



def now():
    return datetime.now().isoformat(timespec="seconds")


def read_json(p: Path):
    try:
        return json.loads(p.read_text())
    except (OSError, ValueError):
        return None


def mark(run_dir: Path, units, **rec):
    state_p = run_dir / "state.json"
    with open(run_dir / ".state.lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        st = json.loads(state_p.read_text()) if state_p.exists() else {"units": {}}
        for u in units:
            st["units"].setdefault(u, {}).setdefault("annotate", {}).update(rec)
        tmp = state_p.with_suffix(".tmp")
        tmp.write_text(json.dumps(st, indent=1))
        tmp.replace(state_p)
        fcntl.flock(lock, fcntl.LOCK_UN)


def place(src: Path, dst: Path):
    if dst.is_symlink() or dst.is_file():
        dst.unlink()
    elif dst.is_dir():
        shutil.rmtree(dst)
    dst.parent.mkdir(parents=True, exist_ok=True)
    if src.is_dir():
        shutil.copytree(src, dst)
    else:
        shutil.copy2(src, dst)


def prepare_workdir(wd: Path, host: str):
    """CeLLMan registered the way its README installs it, in the session's working directory."""
    link = wd / "CeLLMan"
    if link.is_symlink() or link.is_file():
        link.unlink()
    elif link.exists():
        shutil.rmtree(link)
    link.symlink_to(CELLMAN)
    if host == "claude":
        for s in ("annotate-mode", "annotate-reference-mode"):
            place(CELLMAN / "skills" / s, wd / ".claude" / "skills" / s)
        for a in ("annotate-expert", "annotate-lead"):
            place(CELLMAN / "agents" / f"{a}.md", wd / ".claude" / "agents" / f"{a}.md")
    else:
        for s in ("annotate-mode", "annotate-reference-mode"):
            place(CELLMAN / "skills" / s, wd / ".agents" / "skills" / s)
        for a in ("annotate-expert", "annotate-lead"):
            place(CELLMAN / "agents" / f"{a}.toml", wd / ".codex" / "agents" / f"{a}.toml")
        (wd / ".codex" / "config.toml").write_text("\n[sandbox_workspace_write]\nnetwork_access = true\n")


def reference_list(prompt: str, A: dict, run_dir: Path, host: str) -> str:
    """The prompt for a reference given as a list of cell types (annotate-reference-mode, `- reference`)."""
    labels, f = A.get("reference_labels"), A.get("reference_labels_file")
    if not labels and not f:
        return prompt
    if labels:
        f = run_dir / "reference_labels.json"
        f.write_text(json.dumps(labels, indent=1) + "\n")
    lead = "/" if host == "claude" else "$"
    prompt = prompt.replace(f"{lead}annotate-mode ", f"{lead}annotate-reference-mode ", 1)
    lines = (f"- reference        : {f}\n"
             f"- reference scope  : {A.get('reference_scope') or 'strict'}\n")
    i = prompt.find("- note")
    return prompt[:i] + lines + prompt[i:] if i >= 0 else prompt + "\n" + lines


def cli_version(exe):
    try:
        return subprocess.run([exe, "--version"], capture_output=True, text=True, timeout=60).stdout.strip()
    except Exception:  # noqa: BLE001
        return None


def outputs(out: Path):
    """(annotated h5ad, label counts) when the unit's annotation is complete, else None."""
    hits = sorted(out.glob("*_annotated.h5ad")) if out.exists() else []
    if not hits or not (out / "celltype_markers.yaml").exists():
        return None
    free = sorted(p.name for p in out.glob("*_free.*"))
    if free:
        print(f"[annotate] {out} holds free-annotation files, which this pipeline does not make: {free}", flush=True)
        return None
    import anndata as ad
    try:
        obs = ad.read_h5ad(hits[0], backed="r").obs
        if "cellman_free_annotation" in obs.columns:
            print(f"[annotate] {hits[0]} has a cellman_free_annotation column, which this pipeline does not make",
                  flush=True)
            return None
        labels = obs["cellman_annotation"].astype(str).value_counts().to_dict()
    except Exception as e:  # noqa: BLE001 -- e.g. a file cut short by an interruption
        print(f"[annotate] {hits[0]} is not readable ({type(e).__name__}: {e})", flush=True)
        return None
    return hits[0], labels


def collect(run_dir: Path, adir: Path, units, many: bool):
    """Every unit's annotation outputs, moved into 05_annotation/ (a link left in units/<unit>/)."""
    missing, labels, outs = [], {}, {}
    for u in units:
        out = run_dir / "units" / u / "annotation_out"
        got = outputs(out)
        if got is None:
            missing.append(u)
            continue
        h5, labels[u] = got
        dst = (adir / u if many else adir) / "annotation_out"
        relocate(out, dst)
        outs[u] = {"annotation_out": str(dst), "annotated_h5ad": str(dst / h5.name),
                   "markers": str(dst / "celltype_markers.yaml")}
    return missing, labels, outs


def print_labels(labels):
    for u, c in labels.items():
        print(f"[annotate] {u}: {len(c)} labels: " + ", ".join(f"{k} {v}" for k, v in list(c.items())[:30]),
              flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--config", required=True, type=Path)
    a = ap.parse_args()
    cfg = json.loads(a.config.read_text())
    run_dir = Path(cfg["run_dir"])
    adir = run_dir.parent / STAGE_DIR["annotate"]
    A = cfg["annotate"]
    host = A.get("host", "claude")
    model = A.get("model")
    cdir = adir / "cellman"
    cdir.mkdir(parents=True, exist_ok=True)

    units = unit_dirs(run_dir, cfg)
    many = units != ["sample"]
    state = json.loads((run_dir / "state.json").read_text()) if (run_dir / "state.json").exists() else {"units": {}}
    for u in units:
        cls = run_dir / "units" / u / f"{cfg['species']}_{cfg['organ']}_cls.h5ad"
        if not cls.exists():
            # species/organ changed after clustering (the annotation re-runs): the clustered
            # anchors are the same file under the earlier tissue name
            prev = (state["units"].get(u, {}).get("cluster", {}).get("summary") or {}).get("h5ad")
            if not prev or not Path(prev).exists():
                raise SystemExit(f"clustering output missing for unit {u}: {cls}")
            shutil.copy2(Path(prev).resolve(), cls)
            print(f"[annotate] {u}: clustered anchors {prev} -> {cls} (tissue renamed)", flush=True)

    # ---- same inputs as the session that wrote the outputs there: use them, start no new session
    inputs = annotation_inputs(cfg, state, units)
    before = read_json(cdir / INPUTS) or {}
    if before.get("inputs") == inputs and all(outputs(run_dir / "units" / u / "annotation_out") for u in units):
        missing, labels, outs = collect(run_dir, adir, units, many)
        info = dict(read_json(cdir / "annotation_summary.json") or {"host": host, "model": model or "(CLI default)"},
                    stage="annotate", units=units, outputs=outs, labels=labels, missing=missing,
                    reused=True, reused_at=now(),
                    reused_note="settings and clustering unchanged since this annotation was made: "
                                "its outputs are used and CeLLMan was not run again")
        (cdir / "annotation_summary.json").write_text(json.dumps(info, indent=1))
        (adir / "summary.json").write_text(json.dumps(info, indent=1))
        (cdir / INPUTS).write_text(json.dumps(dict(before, status="done", reused_at=now()), indent=1))
        mark(run_dir, units, status="done", finished=now(), reused=True, summary={"labels": labels})
        print(f"[annotate] settings and clustering unchanged since the annotation started "
              f"{before.get('started') or '(earlier run)'}: its outputs are used, no new {host} session", flush=True)
        print_labels(labels)
        return

    for u in units:
        out = run_dir / "units" / u / "annotation_out"
        target = out.resolve() if out.is_symlink() else out
        if out.is_symlink():
            out.unlink()
        if target.exists() and any(target.iterdir()):
            aside = target.with_name(f"annotation_out.prev-{datetime.now():%Y%m%d-%H%M%S}")
            target.rename(aside)
            print(f"[annotate] {u}: previous annotation kept at {aside}", flush=True)

    _, problem = login(host)
    if problem:
        raise SystemExit(problem)

    # the device the agent is told to use is the one attached to this run
    cfg_rt = dict(cfg, gpu_device=resolve_device(cfg))

    wd = run_dir
    prepare_workdir(wd, host)
    prompt = build_prompt(cfg_rt, run_dir, host, workdir=wd)
    prompt = reference_list(prompt, A, run_dir, host)
    (cdir / "prompt.md").write_text(prompt)
    kg = {"command": "python3", "args": [str(CELLMAN / "mcp" / "cellman_kg_mcp.py")],
          "env": {"CELLMAN_PYTHON": PY, "CELLMAN_HELPER": str(CELLMAN / "scripts" / "helper.py"),
                  "CELLONTO_SERVER_URL": cfg["cellman_server"]}}

    if host == "claude":
        exe = shutil.which("claude")
        mcp = cdir / "mcp.json"
        mcp.write_text(json.dumps({"mcpServers": {"cellman-kg": kg}}, indent=1))
        cmd = [exe or "claude", "-p", prompt, "--output-format", "stream-json", "--verbose",
               "--permission-mode", "acceptEdits",
               "--allowedTools", "Bash,Read,Write,Edit,Glob,Grep,Agent,SendMessage,Skill,mcp__cellman-kg",
               "--add-dir", str(run_dir), "--mcp-config", str(mcp)]
        if model:
            cmd += ["--model", model]
    else:
        exe = shutil.which("codex")
        toml_env = "{" + ", ".join(f'{k} = "{v}"' for k, v in kg["env"].items()) + "}"
        cmd = [exe or "codex", "exec"] + (["-m", model] if model else []) + \
              ["-C", str(wd), "-s", A.get("codex_sandbox", "danger-full-access"), "--json", "--skip-git-repo-check",
               "-c", "sandbox_workspace_write.network_access=true",
               "-c", f'mcp_servers.cellman-kg.command="{kg["command"]}"',
               "-c", f'mcp_servers.cellman-kg.args=["{kg["args"][0]}"]',
               "-c", f"mcp_servers.cellman-kg.env={toml_env}",
               prompt]
    if exe is None:
        raise SystemExit(f"{host} CLI not found on PATH")

    env = dict(os.environ)
    env["CLAUDE_CODE_PRINT_BG_WAIT_CEILING_MS"] = "0"
    env["CELLONTO_SERVER_URL"] = cfg["cellman_server"]   # the helper's server, the one the panel stage reads
    env["PATH"] = f"{Path(PY).parent}:{env.get('PATH', '')}"
    t0 = datetime.now()
    (cdir / INPUTS).write_text(json.dumps({"inputs": inputs, "status": "running",
                                           "started": t0.isoformat(timespec="seconds")}, indent=1))
    log = cdir / "caller.jsonl"
    start = log.stat().st_size if log.exists() else 0   # caller.jsonl is appended: this session starts here
    print(f"[annotate] {host} ({cli_version(exe)}) session started for {units} with model "
          f"{model or '(CLI default)'}; events -> {log}", flush=True)
    with open(log, "ab") as fh:
        proc = subprocess.Popen(cmd, cwd=wd, stdout=fh, stderr=subprocess.STDOUT, stdin=subprocess.DEVNULL,
                                env=env)
        (cdir / "caller.pid").write_text(str(proc.pid))
        mark(run_dir, units, status="running", started=t0.isoformat(timespec="seconds"), pid=proc.pid,
             log=str(log), host=host)
        rc = proc.wait()
    (cdir / "caller.pid").unlink(missing_ok=True)

    with open(log, "rb") as fh:
        fh.seek(start)
        evs = events(fh.read().decode(errors="replace"))
    # the session's last message (claude: the result event; codex: the last agent message)
    last = None
    for ev in evs:
        if ev.get("type") == "result":
            last = str(ev.get("result", ""))
        item = ev.get("item") or ev.get("msg") or {}
        if isinstance(item, dict) and item.get("type") in ("agent_message", "assistant_message") and item.get("text"):
            last = item["text"]
    (cdir / "caller_result.txt").write_text(last or "")
    used = used_models(evs, host, model)
    print(f"[annotate] model used: {used['model_used']}"
          + (f" (models called: {', '.join(used['models_called'])})" if len(used["models_called"]) > 1 else ""),
          flush=True)

    missing, labels, outs = collect(run_dir, adir, units, many)
    secs = round((datetime.now() - t0).total_seconds())
    info = {"stage": "annotate", "host": host, "cli_version": cli_version(exe), "model": model or "(CLI default)",
            "model_used": used["model_used"], "models_called": used["models_called"],
            "reported_by_cli": used.get("reported_by_cli"), "model_source": used.get("model_source"),
            "exit_code": rc, "seconds": secs, "units": units, "outputs": outs, "labels": labels,
            "missing": missing, "prompt": str(cdir / "prompt.md"),
            "events": str(log), "last_message": (last or "")[-4000:]}
    (cdir / "annotation_summary.json").write_text(json.dumps(info, indent=1))
    (adir / "summary.json").write_text(json.dumps(info, indent=1))
    if missing:
        mark(run_dir, units, status="failed", error=f"no annotation outputs for {missing}", finished=now(), seconds=secs)
        raise SystemExit(f"annotation ended (exit {rc}) without outputs for {missing}. "
                         f"Last message of the {host} session:\n{(last or '(none)')[-3000:]}")
    (cdir / INPUTS).write_text(json.dumps({"inputs": inputs, "status": "done",
                                           "started": t0.isoformat(timespec="seconds"), "finished": now()}, indent=1))
    mark(run_dir, units, status="done", finished=now(), seconds=secs, summary={"labels": labels},
         model_used=used["model_used"])
    print_labels(labels)


if __name__ == "__main__":
    main()
