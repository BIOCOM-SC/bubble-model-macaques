# runs/ - the simulation runs of the paper

Every result of the paper comes from one run, described by a config file `configs/<name>.toml`
(merged over `configs/_base.toml`) and executed by `launch.py` in its own work directory
`work/<name>/` (not included in the repository). Python >= 3.11, standard library only. Fits need
[Dakota](https://dakota.sandia.gov) (DIRECT, `coliny_direct`).

## Runs

Fits (Dakota DIRECT; <= 5000 global + 5 x 500 local + 500 refinement evaluations, each of 10^4
simulated animals):

- `fitA`: 7 free parameters (bounds of Table 1), objective E_A.
- `fitB`: the same with objective E_B.
- `noimm`: alternative model without adaptive immunity (f_imm = 1), 5 free parameters.
- `timm21`: alternative model with immunity starting at t_imm = 21 days, 6 free parameters.

Case sets at the fit-A optimum (number of cases, each of 10^4 simulated animals):

- `thresholds` (9): detection thresholds, micronodule {0.5, 1, 1.5} x consolidation
  {3.5, 4.5, 5.5} mm (compile-time).
- `experiments` (23): control, immunosuppressed (f_imm = 1), lesions limited by segment walls,
  f_imm sweep (20 values).
- `immsplit` (3): immunity on growth rate only / on maximum radius only / both.
- `sweeps` (160): single-parameter sweeps, +/- 25 %, 20 values, 8 parameters.
- `immgrid` (874): t_imm x f_imm grid (t_imm 0-45 days, f_imm 0.10-1.00): Fig. 4 and the
  vaccinated macaques.
- `prcc` (1000): PRCC, Latin-hypercube samples, +/- 25 % around fit A, 8 parameters.

Group: `chainA` runs `fitA` and then every case set above, in one job.

All case sets read the fit-A optimum (`base = "fitA"` -> `work/fitA/best.json`), so `fitA` must
finish first; `chainA` does this in one job.

## Usage

Cluster settings (account, queues, modules, Dakota path, time limit of each job) go in
`hpc_site.toml` at the repository root; copy `hpc_site.example.toml` and fill it in.

```sh
python3 runs/launch.py submit fitA --test        # check the setup: small fit on the test queue
python3 runs/launch.py submit all                # the paper's runs: chainA, fitB, noimm, timm21
python3 runs/launch.py list                      # state of every run and of its Slurm job
python3 runs/launch.py run experiments --local   # run now on this machine, without Slurm
python3 runs/launch.py --help                    # all commands and options
```

`--test` runs a small version (500 simulations per evaluation, a few evaluations, 3 cases) in
`work/test/`, to check a setup in minutes; its numbers are meaningless. Case sets in test mode need
the test fit A (`run fitA --test`).

On the cluster used for the paper (2 nodes x 112 cores per job, 2 MPI processes x 56 threads per
node) a fit evaluation took 3-5 s: each fit 6-11 h, the case sets about 3.5 h in total. Re-running
or re-submitting a run continues it: fits restart the unfinished DIRECT search from Dakota's
`dakota.rst`, case sets skip finished cases.

## How a fit works

Each evaluation averages 10^4 simulations, so the objective is noisy (SD ~ 0.002 at the fit-A
optimum, E_A ~ 0.22). A fit has three steps:

1. **Global DIRECT** (Dakota `coliny_direct`) over the whole search box, in `global/`, at most
   5000 evaluations. `local_balance_parameter = 0.035` stops dividing boxes three trisection
   levels below the largest box: promising regions are resolved to about 1/27 of each parameter
   range, finer resolution is left to step 2. `division = all_dimensions` divides every longest
   side of a box, as in Jones et al. (1993) (Dakota's default divides the first parameter far
   more often than the others).
2. **Local DIRECT** in a box of +/- 0.1 x range around each of the 5 best distinct points of the
   global search (at least 0.2 x range apart in some parameter), in `local<i>/`, 500 evaluations
   each, same options.
3. **Refinement**: the 10 best distinct points of each local search are re-evaluated 10 times
   (`refine.csv`); the optimum is the point with the lowest mean (the minimum of single noisy
   evaluations is biased low).

If the model aborts an evaluation (lesion limits or 600 s, see `model/README.md`), Dakota receives
a penalty of 10^6. `best.json` holds the optimum (mean, SD and per-output errors), the noise of the
objective, all candidates and a summary of every search; `best_results.txt` the daily curves of the
optimum.

## Files

- `launch.py` - the commands above; `bm.py` - shared code (configs, parameters, building and
  running the model, case generation); `driver.py` - called by Dakota for each evaluation.
- `templates/` - the Dakota input and the Slurm job script, filled in by `launch.py`.
- `configs/` - one file per run; `_base.toml` holds the fixed parameters (Table 3 of the paper).

Every evaluation has its own seed (a hash of the run name and the evaluation), recorded in the
outputs, so any evaluation can be repeated exactly. A fit, and a case set with a single model build,
send all evaluations to one model process (`bubble-model -server`) instead of starting the model
each time; the results are the same.

Work directory of a run (`work/<name>/`): `config.json` (resolved config), `manifest.json` (git
commit, build commands, host), `model.log`; fits: `global/`, `local<i>/` (Dakota files and
`evals.csv`, one line per evaluation), `starts.json`, `refine.csv`, `best.json`, `best_params.txt`;
case sets: `cases.json`, `cases_out.csv` (one line per case: parameters, errors, status, time,
seed) and `results/cases/<id>_<name>.txt` (the model's `results.txt` for each case).
