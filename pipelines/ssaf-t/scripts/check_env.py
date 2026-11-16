#!/usr/bin/env python
"""First rule: is the chosen AI CLI installed, logged in and answering with the chosen model; which
GPU / CPU does this run get.

    check_env.py --host claude|codex [--model M] [--gpu auto|cpu|cuda] --out 00_check/env.json
                 --gpu-info 00_check/gpu_info.txt

Fails before any heavy stage when the container holds no login for the CLI (AutoPipe mounts the
server user's login from version 0.0.36 on), when the CLI cannot answer one short prompt with the
model asked for (llm_model; blank = the CLI's default model), or when gpu: cuda is requested and no
CUDA device is usable. env.json records the model that answered (model_used), and the server GPUs
the run was given (server_gpus: AutoPipe's run record .autopipe-run.json, `gpus`), since inside the
container the attached GPUs are numbered from 0.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

AUTOPIPE_MIN = "0.0.36"
PROMPT = "Reply with the single word OK."
RUN_RECORD = ".autopipe-run.json"


def login(host: str) -> tuple[Path, str | None]:
    """The CLI's login file in this container, and why the run cannot use it (None: it can)."""
    if host == "claude":
        name, var, key, how = "Claude Code", "CLAUDE_CONFIG_DIR", "ANTHROPIC_API_KEY", "claude"
        f = Path(os.environ.get(var, Path.home() / ".claude")) / ".credentials.json"
    else:
        name, var, key, how = "Codex", "CODEX_HOME", "OPENAI_API_KEY", "codex login"
        f = Path(os.environ.get(var, Path.home() / ".codex")) / "auth.json"
    if f.exists() or os.environ.get(key):
        return f, None
    return f, (f"No {name} login in the container ({f} does not exist). The annotation uses the server user's "
               f"{name} login, which AutoPipe mounts into the container from version {AUTOPIPE_MIN} on: run this "
               f"pipeline with AutoPipe {AUTOPIPE_MIN} or later. With such a version, check that the server user "
               f"is logged in (`{how}` on the server), then run again.")


def server_gpus(start: Path | None = None) -> dict | None:
    """The server GPUs AutoPipe attached to this run: `gpus` of its run record (.autopipe-run.json) in
    the output folder, or the batch's output folder above it. None when there is no record."""
    here = Path(start or Path.cwd()).resolve()
    for d in [here, *here.parents]:
        f = d / RUN_RECORD
        if f.is_file():
            try:
                j = json.loads(f.read_text())
            except (OSError, ValueError):
                return None
            return {"gpus": j.get("gpus"), "gpu_flags": j.get("gpu_flags"), "run_name": j.get("run_name"),
                    "started_at": j.get("started_at"), "record": str(f)}
    return None


def events(text: str) -> list[dict]:
    """The JSON objects of a CLI's JSON output (one object, or one per line)."""
    out = []
    for line in text.splitlines():
        line = line.strip()
        if line.startswith("{"):
            try:
                ev = json.loads(line)
            except ValueError:
                continue
            if isinstance(ev, dict):
                out.append(ev)
    return out


def _fields(obj, key):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == key and isinstance(v, str) and v:
                yield v
            yield from _fields(v, key)
    elif isinstance(obj, list):
        for v in obj:
            yield from _fields(v, key)


def codex_config_model() -> str | None:
    """`model = "..."` of the mounted Codex config.toml, if set."""
    f = Path(os.environ.get("CODEX_HOME", Path.home() / ".codex")) / "config.toml"
    try:
        m = re.search(r'^\s*model\s*=\s*"([^"]+)"', f.read_text(), re.M)
        return m.group(1) if m else None
    except OSError:
        return None


def codex_session_models(evs: list[dict]) -> list[str]:
    """Models of a Codex session and of every subagent it spawned, read from the session files Codex
    writes under $CODEX_HOME/sessions/YYYY/MM/DD/rollout-*-<thread id>.jsonl. `codex exec --json` does
    not report the model, but each `turn_context` line of a session file does, and the first line of a
    subagent's file (session_meta) names its parent thread. The session is found by the thread id of
    its `thread.started` event, so a concurrent run on the same login is never mixed in."""
    tid = next((str(e["thread_id"]) for e in evs if e.get("type") == "thread.started" and e.get("thread_id")), None)
    root = Path(os.environ.get("CODEX_HOME", Path.home() / ".codex")) / "sessions"
    if not tid or not root.is_dir():
        return []
    main = sorted(root.glob(f"*/*/*/rollout-*-{tid}.jsonl"))
    if not main:
        return []
    # subagents start after the main session: its day folder and any later ones (a run past midnight)
    day0 = main[0].parent
    days = sorted(d for d in root.glob("*/*/*") if d.is_dir() and d.relative_to(root).parts >= day0.relative_to(root).parts)
    parent = {}
    for d in days:
        for f in d.glob("rollout-*.jsonl"):
            try:
                with open(f, encoding="utf-8", errors="replace") as fh:
                    meta = json.loads(fh.readline()).get("payload") or {}
            except (OSError, ValueError):
                continue
            if meta.get("id"):
                parent[str(meta["id"])] = (str(meta.get("parent_thread_id") or ""), f)
    keep, grew = {tid}, True
    while grew:   # the main thread and its descendants (a lead's judges are depth 2)
        grew = False
        for i, (p, _) in parent.items():
            if p in keep and i not in keep:
                keep.add(i)
                grew = True
    files = [main[0]] + sorted(f for i, (_, f) in parent.items() if i in keep and f != main[0])
    seen = []
    for f in files:
        try:
            with open(f, encoding="utf-8", errors="replace") as fh:
                for line in fh:
                    if '"turn_context"' not in line:
                        continue
                    m = (json.loads(line).get("payload") or {}).get("model")
                    if m and m not in seen:
                        seen.append(str(m))
        except (OSError, ValueError):
            continue
    return seen


def used_models(evs: list[dict], host: str, asked: str | None) -> dict:
    """The model that ran a session, from its JSON events.

    claude: the init event's `model` (the session's model) and the result's `modelUsage` (every model
    the session called, e.g. a subagent's). codex: any `model` field of its events; when the CLI reports
    none, the models in Codex's own session files (codex_session_models: the session and its
    subagents); when those cannot be read either, the model asked for or the config.toml default.
    Neither is reported by the CLI (`reported_by_cli: false`); `model_source` says where it came from."""
    if host == "claude":
        main = next((e["model"] for e in evs if e.get("type") == "system" and e.get("model")), None)
        every = []
        for e in evs:
            every += list((e.get("modelUsage") or {}).keys())
        every = list(dict.fromkeys(([main] if main else []) + every))
        return {"model_used": main or (every[0] if len(every) == 1 else None), "models_called": every,
                "reported_by_cli": bool(every), "model_source": "cli events" if every else "unknown"}
    every = list(dict.fromkeys(m for e in evs for m in _fields(e, "model")))
    if every:
        return {"model_used": every[0], "models_called": every, "reported_by_cli": True, "model_source": "cli events"}
    files = codex_session_models(evs)
    if files:
        return {"model_used": files[0], "models_called": files, "reported_by_cli": False,
                "model_source": "codex session files"}
    guess = asked or codex_config_model()
    return {"model_used": guess or "(codex default; the CLI did not report its name)", "models_called": [],
            "reported_by_cli": False,
            "model_source": "llm_model" if asked else ("config.toml" if guess else "unknown")}


def reply_text(evs: list[dict], raw: str) -> str:
    """The answer of a one-prompt session (claude: result; codex: the agent message)."""
    got = [str(e["result"]) for e in evs if e.get("type") == "result" and "result" in e]
    for e in evs:
        item = e.get("item") or e.get("msg") or {}
        if isinstance(item, dict) and item.get("text"):
            got.append(str(item["text"]))
    return " ".join(got) if got else raw


def run(cmd, timeout=300, cwd=None):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, cwd=cwd,
                           stdin=subprocess.DEVNULL)
        return r.returncode, r.stdout, r.stderr
    except subprocess.TimeoutExpired:
        return 124, "", f"timed out after {timeout}s"
    except OSError as e:
        return 127, "", str(e)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--host", required=True, choices=["claude", "codex"])
    ap.add_argument("--model", default="")
    ap.add_argument("--gpu", default="auto")
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--gpu-info", required=True, type=Path)
    a = ap.parse_args()
    rep = {"host": a.host, "model": a.model or "(CLI default)"}
    a.out.parent.mkdir(parents=True, exist_ok=True)

    # ---- GPU / CPU (written first, so it is there even when the CLI check fails)
    rc, full, _ = run(["nvidia-smi"], 60)
    a.gpu_info.write_text(full if rc == 0 else "no GPU attached to this run (nvidia-smi unavailable)\n")
    rc, out, _ = run(["nvidia-smi", "-L"], 60)
    rep["nvidia_smi"] = out.strip().splitlines() if rc == 0 else []
    # the container numbers its GPUs from 0; the server's numbers are those the run was given
    rep["server_gpus"] = server_gpus()

    # ---- AI CLI: logged in, and answering with the model asked for
    exe = shutil.which(a.host)
    if exe is None:
        sys.exit(f"{a.host} CLI is not installed in the image")
    rep["cli_version"] = run([exe, "--version"], 60)[1].strip()
    f, problem = login(a.host)
    rep["login_file"] = str(f)
    if problem:
        sys.exit(problem)
    if a.host == "claude":
        cmd = [exe, "-p", PROMPT, "--max-turns", "1", "--output-format", "stream-json", "--verbose"]
        cmd += ["--model", a.model] if a.model else []
    else:
        cmd = [exe, "exec"] + (["-m", a.model] if a.model else []) + \
              ["--skip-git-repo-check", "-s", "read-only", "--json", PROMPT]
    with tempfile.TemporaryDirectory() as td:
        rc, out, err = run(cmd, 300, cwd=td)
    evs = events(out)
    reply = reply_text(evs, out)
    errored = any(e.get("type") == "result" and e.get("is_error") for e in evs)
    rep.update(used_models(evs, a.host, a.model or None))
    rep["cli_check"] = {"exit_code": rc, "reply": reply.strip()[-300:]}
    if rc != 0 or errored or "OK" not in reply.upper():
        why = (reply.strip() if evs else "") or err.strip() or out.strip()
        why = why[-1200:] or "(no output)"
        if a.model:
            sys.exit(f"{a.host} could not answer with the model {a.model!r} on the server user's account "
                     f"(exit {rc}): {why}\nCheck the model name and that this account may use it, or leave "
                     f"llm_model blank to use the CLI's default model.")
        sys.exit(f"{a.host} did not answer a one-line prompt with its default model (exit {rc}): {why}")
    if a.model and rep.get("reported_by_cli") and rep.get("model_used") and \
            a.model not in rep["model_used"] and rep["model_used"] not in a.model:
        rep["model_note"] = (f"llm_model {a.model!r} is answered by {rep['model_used']!r} "
                             f"(the CLI's name for it)")

    try:
        import torch
        ok = torch.cuda.is_available()
        rep["torch"] = {"version": torch.__version__, "cuda": torch.version.cuda, "cuda_available": ok,
                        "devices": [torch.cuda.get_device_name(i) for i in range(torch.cuda.device_count())] if ok else []}
    except Exception as e:  # noqa: BLE001
        ok = False
        rep["torch"] = {"error": f"{type(e).__name__}: {e}"}
    rep["gpu_request"] = a.gpu
    rep["device"] = "cuda:0" if (ok and a.gpu != "cpu") else "cpu"
    if a.gpu == "cuda" and not ok:
        sys.exit("gpu: cuda was requested but no CUDA device is usable (no GPU attached to the run, or the "
                 "server driver does not support CUDA 13). Attach a GPU or set gpu: auto / cpu.")
    if rep["nvidia_smi"] and not ok:
        rep["warning"] = ("a GPU is attached but torch cannot use it (the CUDA 13 stack needs a driver >= 580); "
                          "the run continues on the CPU")
    rep["cpu_cores"] = os.cpu_count()
    try:
        rep["avx512f"] = "avx512f" in Path("/proc/cpuinfo").read_text()
    except OSError:
        rep["avx512f"] = None
    try:
        import ssaft
        rep["ssaft"] = ssaft.__version__
    except Exception as e:  # noqa: BLE001
        sys.exit(f"ssaft does not import in the pipeline environment: {e}")
    a.out.write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
