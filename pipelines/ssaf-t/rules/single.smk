# One sample. Stages: load -> zstack -> tile -> panel_vf -> downsample -> cluster ->
# annotate -> map -> filter -> segment; which of them run depends on the platform, the map's dimension,
# the image downsampling uses (a z-stack adds the zstack stage), tiling and segmentation. Each stage
# writes into its numbered folder (scripts/layout.py: 00_check ... 09_final).
# Each stage gets the run.json sections it reads as params, so a changed value re-runs that stage and
# everything after it; earlier stages are kept.
# Params are functions throughout: run.json values carry braces (e.g. "mosaic_{stain}_z*.tif") that
# Snakemake would otherwise read as wildcards. The target (report.md) is rule all in the Snakefile.
#
# No stage depends on a job that can run again: run.json is written here when a run starts, not by a
# rule, and check_env runs first but is not an input of any stage. (When a job whose output others
# take as ancient() input runs again, Snakemake does not look at those others in that run: changed
# params and even missing outputs below it would go unnoticed.)

# ---- plan (printed at dry_run): inspect the input and build run.json (scripts/plan.py).
# _cores: the cores this run was given, used when config threads is blank
_plan = subprocess.run([PY, f"{PIPE}/scripts/plan.py", "--config-json", json.dumps(dict(config, _cores=workflow.cores)),
                        "--run-dir", WORK], capture_output=True, text=True)
sys.stderr.write(_plan.stderr)
if _plan.returncode != 0:
    raise WorkflowError(_plan.stdout.strip() or _plan.stderr.strip()[-4000:])
PLAN = json.loads(_plan.stdout)
CFG = PLAN["cfg"]
STAGES = CFG["stages"]
STEPS = [STAGE_DIR[s] for s in STAGES]          # e.g. 01_load, 03_downsample, ...
HOST = CFG["annotate"]["host"]
NCORES = int(CFG["ncores"])                     # config threads, or the cores the run was given


# ---- run.json for the stage scripts (in .ssaft/) and a copy of the plan in 00_check/
def _read(path):
    try:
        with open(path) as fh:
            return fh.read()
    except OSError:
        return None


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path + ".tmp", "w") as fh:
        fh.write(text)
    os.replace(path + ".tmp", path)


_DRY = any(a in ("-n", "--dry-run", "--dryrun") or re.fullmatch(r"-[A-Za-z]*n[A-Za-z]*", a) for a in sys.argv[1:])
_LOCKED = os.path.isdir(".snakemake/locks") and bool(os.listdir(".snakemake/locks"))
_FILES = {f"{WORK}/run.json": json.dumps(CFG, indent=1, default=str),
          f"{CHECK}/run.json": json.dumps(CFG, indent=1, default=str),
          f"{CHECK}/inspect.json": json.dumps(PLAN["inspect"], indent=1, default=str),
          f"{CHECK}/answers.json": json.dumps(PLAN["answers"], indent=1, default=str),
          f"{CHECK}/plan.txt": PLAN["summary"] + "\n"}
if PLAN.get("zplanes"):
    # molecules per z plane, the recommended range and the planes this run keeps (a z-stack)
    _FILES[f"{CHECK}/zstack_planes.json"] = json.dumps(PLAN["zplanes"], indent=1, default=str)
if PLAN.get("images"):
    # the images downsampling can use, the recommended one and the one this run uses
    _FILES[f"{CHECK}/images.json"] = json.dumps(PLAN["images"], indent=1, default=str)
# a dry run (or a run that cannot start because another one holds this folder) only creates what is
# missing; a run rewrites what changed
for _path, _text in _FILES.items():
    _old = _read(_path)
    if _old != _text and (_old is None or not (_DRY or _LOCKED)):
        _write(_path, _text)


def stage_params(stage):
    """The run.json values a stage reads (a change re-runs that stage and what follows)."""
    c = CFG
    if stage == "load":
        p = {"dataset": c["dataset_path"], "platform": c["platform"], "load": c["load"],
             "image": {k: c["image"].get(k) for k in ("key", "channel")}, "input_kind": c["input_kind"]}
        if c["platform"] in ("merscope", "cosmx") or c["input_kind"] == "zarr":
            p["dimension"] = c["dimension"]
        if c.get("merscope_images"):
            p["merscope_images"] = c["merscope_images"]         # the 2D mosaic and the z-stack written
        if c["input_kind"] == "zarr":
            p["zstack_element"] = c["zarr"]["zstack_element"]
            if c["zarr"]["zstack_element"] and c["zstack"].get("planes"):
                p["zstack_planes"] = c["zstack"]["planes"]      # the element is cut down at load
        return p
    if stage == "zstack":
        return {"zstack": c["zstack"], "dimension": c["dimension"], "tiled": "tile" in STAGES}
    if stage == "tile":
        return {"tile": c["tile"], "image": c["image"], "dimension": c["dimension"]}
    if stage == "panel_vf":
        return {"species": c["species"], "organ": c["organ"], "panel": c["panel"], "vf": c["vf"],
                "server": c["cellman_server"]}
    if stage == "downsample":
        p = {"downsample": c["downsample"], "image": c["image"], "dimension": c["dimension"],
             "tiled": "tile" in STAGES}
        if c.get("downsample_image"):
            p["downsample_image"] = c["downsample_image"]      # the 2D image or the z-stack sampled
        return p
    if stage == "cluster":
        return {"cluster": c["cluster"]}
    if stage == "annotate":
        return {"species": c["species"], "organ": c["organ"], "condition": c["condition"],
                "annotate": c["annotate"]}
    return {stage: c[stage]}


def stage_values(stage):
    return json.dumps(stage_params(stage), sort_keys=True, default=str)


def done(step):
    return f"{WORK}/stages/{step}.done"


def prev_done(step):
    i = STEPS.index(step)
    return done(STEPS[i - 1]) if i else []


wildcard_constraints:
    step="|".join(re.escape(s) for s in STEPS if s != STAGE_DIR["annotate"])


rule check_env:
    """The chosen AI CLI answers one prompt; GPU / CPU of this run recorded. Runs before the stages
    (highest priority), and a failure stops the run before them."""
    output:
        env=f"{CHECK}/env.json",
        gpu=f"{CHECK}/gpu_info.txt",
    params:
        host=lambda w: HOST,
        model=lambda w: CFG["annotate"].get("model") or "",
        gpu=lambda w: CFG["gpu_request"],
    priority: 100
    log:
        f"{CHECK}/check_env.log",
    shell:
        "{PY} {PIPE}/scripts/check_env.py --host {params.host} --model '{params.model}' --gpu {params.gpu}"
        " --out {output.env} --gpu-info {output.gpu} 2>&1 | tee {log}"


rule stage:
    """One SSAF-T stage for every unit (ssaft_stage.py); outputs into its folder."""
    input:
        cfg=ancient(f"{WORK}/run.json"),
        prev=lambda w: prev_done(w.step),
    output:
        f"{WORK}/stages/{{step}}.done",
    params:
        stage=lambda w: DIR_STAGE[w.step],
        values=lambda w: stage_values(DIR_STAGE[w.step]),
    threads: NCORES
    log:
        "{step}/stage.log",
    shell:
        "{PY} {PIPE}/scripts/ssaft_stage.py --config {input.cfg} --stage {params.stage} --cores {threads}"
        " 2>&1 | tee {log}; touch {output}"


rule annotate:
    """CeLLMan annotation by Claude Code / Codex, headless (annotate.py)."""
    input:
        cfg=ancient(f"{WORK}/run.json"),
        prev=prev_done(STAGE_DIR["annotate"]),
    output:
        done(STAGE_DIR["annotate"]),
    params:
        values=lambda w: stage_values("annotate"),
    threads: NCORES
    log:
        f"{STAGE_DIR['annotate']}/annotate.log",
    shell:
        "{PY} {PIPE}/scripts/annotate.py --config {input.cfg} 2>&1 | tee {log}; touch {output}"


rule final:
    """09_final: the SpatialData with every input and SSAF-T element."""
    input:
        cfg=ancient(f"{WORK}/run.json"),
        last=done(STEPS[-1]),
    output:
        done(FINAL),
    log:
        f"{FINAL}/final.log",
    shell:
        "{PY} {PIPE}/scripts/ssaft_stage.py --config {input.cfg} --stage final 2>&1 | tee {log}; touch {output}"


rule report:
    """Per unit: cell types, kept voxel fraction, cells; output locations; map preview."""
    input:
        cfg=ancient(f"{WORK}/run.json"),
        final=done(FINAL),
    output:
        "report.md",
        "report.json",
    log:
        f"{FINAL}/report.log",
    shell:
        "{PY} {PIPE}/scripts/report.py --config {input.cfg} 2>&1 | tee {log}"
