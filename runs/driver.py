#!/usr/bin/env python3
"""
Dakota analysis driver for the Bubble Model fits (replaces prepro.py + params_bm.txt.stencil +
output_for_dakota_err*).

Called by Dakota (fork interface) inside the run's work directory as
    driver.sh <parameters_file> <results_file>
It
  1. reads the free parameters from Dakota's parameters file and combines them with the fixed and
     tied parameters of the run (config.json, written by launch.py) -> config/params_bm.txt;
  2. runs one model evaluation, `-N <n_sims> -seed <seed>` with seed = hash of (search name,
     eval_id), so a re-evaluated point gets the same result: through the model server started by
     launch.py, or `bubble-model -noplot -d<domain> -N <n_sims> -seed <seed>` (through the site
     launcher, e.g. srun) if none is running;
  3. reads results/error.txt and writes the chosen objective (E_A or E_B) to the results file.
     If the model aborted (MAX_LESIONS / MAX_* limits / MAX_SIMU_TIME: error.txt holds the
     integer line '100, 100, ...') it returns the run's penalty instead (config fit.penalty), so
     that DIRECT treats the point as very poor but stays finite.
     If the run did not complete at all (crash, or the launcher failed, e.g. srun at the job's
     time limit) nothing is written: the driver exits with an error, Dakota stops without
     recording the evaluation, and a restart evaluates it again.
  4. appends one line per evaluation to evals.csv (eval id, parameters, all errors, status,
     seconds, seed). evals.csv is the record used to pick the optimum (minimum objective).
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bm  # noqa: E402


def main():
    params_file, results_file = sys.argv[1], sys.argv[2]
    work = Path.cwd()
    cfg = json.loads((work / "config.json").read_text())
    fit = cfg["fit"]

    free, eval_id = bm.read_dakota_params(params_file)
    if set(free) != set(fit["free"]):
        sys.exit(f"driver: Dakota variables {sorted(free)} != free parameters {sorted(fit['free'])}")
    p = {**cfg["resolved_params"], **free}
    bm.apply_tied(p, cfg.get("tied", {}), explicit=set(cfg.get("params", {})) | set(free))

    out = bm.run_model(work, work / "bin" / "bubble-model", p, cfg["n_sims"], cfg.get("domain", "lobe"),
                       cfg["launcher"], work / "model.log", bm.seed_for(cfg["name"], eval_id))
    if out["status"] == bm.CRASH:
        sys.exit(f"driver: model run of evaluation {eval_id} did not complete (return code "
                 f"{out['returncode']}, see model.log); not recorded")
    if out["status"] == bm.OK:
        value = out["errors"][fit["objective"]]
    else:
        value = float(fit["penalty"])
    Path(results_file).write_text(f"{value:.10e} {fit['objective']}\n")

    errs = out["errors"] or {}
    bm.append_row(work / "evals.csv",
                  ["eval_id", *bm.PARAM_NAMES, *bm.ERROR_NAMES, "objective", "status", "seconds", "seed"],
                  [eval_id, *(repr(p[k]) for k in bm.PARAM_NAMES), *(errs.get(k, "") for k in bm.ERROR_NAMES),
                   repr(value), out["status"], out["seconds"], out["seed"]])


if __name__ == "__main__":
    main()
