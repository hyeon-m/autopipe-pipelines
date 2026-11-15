# A batch. config.yaml holds the shared values, `samples:` or `manifest:` the samples and what
# differs per sample.
# Every sample is planned at parse time (batch_plan.py): its plan is printed, a sample that does
# not validate is reported and skipped. scripts/batch.py then runs the samples one after another, each
# with the single-sample workflow in runs/<name>/. Re-running the batch re-runs each
# sample only as far as a changed value requires. The target (.batch/batch.done) is rule all in the
# Snakefile.

# _cores: the cores this run was given, used when config threads is blank
_bplan = subprocess.run([PY, f"{PIPE}/scripts/batch_plan.py", "--config-json",
                         json.dumps(dict(config, _cores=workflow.cores))], capture_output=True, text=True)
sys.stderr.write(_bplan.stderr)
if _bplan.returncode != 0:
    raise WorkflowError(_bplan.stdout.strip() or _bplan.stderr.strip()[-4000:])
BPLAN = json.loads(_bplan.stdout)


rule batch_configure:
    """The batch plan: every sample's config and validation (batch_plan.py)."""
    output:
        "batch_plan.json",
    params:
        plan=lambda w: json.dumps(BPLAN, sort_keys=True, default=str),
    run:
        with open(output[0], "w") as fh:
            json.dump(BPLAN, fh, indent=1, default=str)
        with open("batch_manifest.tsv", "w") as fh:
            fh.write("order\tname\tinput\tstatus\tplatform\tdimension\ttissue\toverrides\terrors\n")
            for s in BPLAN["samples"]:
                fh.write("\t".join(str(x) for x in (
                    s["order"], s["name"], s["input"], "ok" if s["ok"] else "FAIL", s.get("platform", ""),
                    s.get("dimension", ""), s.get("tissue", ""), json.dumps(s["overrides"], default=str),
                    " | ".join(s["errors"]).replace("\n", " "))) + "\n")


rule batch:
    """Every sample, one after another, from load to its report (batch.py)."""
    input:
        "batch_plan.json",
    output:
        touch(".batch/batch.done"),
    threads: workflow.cores
    log:
        ".batch/last_run.log",
    shell:
        # batch.log is not the rule's log: Snakemake deletes a rule's log when the rule starts, and a
        # resumed batch must keep what the earlier attempt wrote
        "{PY} {PIPE}/scripts/batch.py --plan {input} --snakefile {SNAKEFILE} --cores {threads} 2>&1 | tee -a batch.log {log}"
