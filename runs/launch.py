#!/usr/bin/env python3
"""
Runs of the Bubble Model paper: fits, in-silico experiments, PRCC, sweeps.

Each run is one file runs/configs/<name>.toml (merged over configs/_base.toml), executed in its own
work directory runs/work/<name>/ (runs/work/test/<name>/ with --test).

Commands
  launch.py run <name> [--local] [--test]    run now, in this shell (inside a Slurm job, or locally)
  launch.py submit <name|all> [--test] [--time HH:MM:SS] [--dry-run]
                                             write runs/work/<name>/job.sh (Slurm) and submit it;
                                             `all` = every config not part of a group
  launch.py list [--long] [--test]           state of every run and of its Slurm job

  --local   build with `make local` and run without a launcher (one MPI process); no Slurm.
  --test    small version of the run, to check a setup in minutes (numbers are meaningless):
            500 simulations per evaluation; fits: 3 global + 2 x 3 local evaluations + 2 x 2
            refinement re-evaluations; case sets: 3 cases. `submit --test` uses the test QOS.
            Case sets based on fit A need the test fit A first (`run fitA --test`).

Configs: type = "fit" (Dakota DIRECT), "cases" (a list of parameter sets) or "group" (several
runs in order in one job, e.g. chainA = fitA and every case set that starts from its optimum).
Cluster settings (account, QOS, modules, Dakota path, job time limits): hpc_site.toml, see
hpc_site.example.toml.

Work directory of a run:
  config.json         resolved config (what was actually run)
  manifest.json       git commit, model-source and targets checksums, build commands, host
  config/             params_bm.txt, experimental_results.txt, macaque -> lung/macaque
  bin/                model binaries (built at run time, one per set of compile-time defs)
  results/            model outputs of the last evaluation
  model.log           model output of every evaluation
  fits:       global/, local<i>/   one DIRECT search each (dakota.in, dakota.*, evals.csv = one
                                   line per evaluation); starts.json = starts of the local searches
              refine.csv           re-evaluations of the best points of each local search
              best.json, best_params.txt   optimum = refined point with the lowest mean objective
  case sets:  cases.json           the cases (parameters, domain, compile-time defs)
              cases_out.csv        one line per finished case (errors, status, seconds, seed)
              results/cases/<id>_<name>.txt   daily mean and SD of every output, per case

Every evaluation has an explicit seed (bm.seed_for). A fit, and a case set with one build, run all
evaluations through one persistent model process (`bubble-model -server`). Re-running or
re-submitting a run continues it: fits restart the unfinished DIRECT search from dakota.rst, case
sets skip the cases already in cases_out.csv.
"""
import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import bm  # noqa: E402

TEMPLATES = bm.RUNS / "templates"


def root_for(test):
    return bm.WORK / "test" if test else bm.WORK


def test_overrides(cfg):
    t = cfg.get("test", {})
    cfg = dict(cfg)
    cfg["n_sims"] = t.get("n_sims", 500)
    if cfg["type"] == "fit":
        cfg["fit"] = {**cfg["fit"], "max_evals": t.get("max_evals", 3), "starts": t.get("starts", 2),
                      "global": {**cfg["fit"]["global"], "max_evals": t.get("global_evals", 4)},
                      "refine_per_search": t.get("refine_per_search", 2), "refine_reps": t.get("refine_reps", 2)}
    cfg["max_cases"] = t.get("max_cases", 3)
    return cfg


def resolve(name, test, local):
    cfg = bm.load_config(name)
    if test:
        cfg = test_overrides(cfg)
    site = None if local else bm.load_site()
    cfg["mode"] = "local" if local else "hpc"
    cfg["test"] = bool(test)
    cfg["launcher"] = "" if local else site.get("launcher", "srun")
    return cfg, site


class ModelServer:
    """Persistent model process for one binary (`with ModelServer(...)`): started in the background
    (through the launcher) in <work>/server/, used by bm.run_model via BM_SERVER, stopped at the end
    (stop request; terminated if it does not stop). If the process dies, the evaluation in flight is
    a crash (not recorded) and the run stops, as when a one-shot model run fails."""

    def __init__(self, work, exe, cfg):
        self.work, self.exe, self.cfg, self.proc = Path(work), Path(exe), cfg, None
        self.dir = self.work / "server"

    def __enter__(self):
        if not self.cfg.get("server", True):
            return self
        if self.dir.exists():
            shutil.rmtree(self.dir)
        self.dir.mkdir()
        (self.dir / "exe").write_text(str(self.exe.resolve()) + "\n")
        cmd = shlex.split(self.cfg["launcher"]) + [str(self.exe), "-noplot", "-server", str(self.dir)]
        self.log = open(self.work / "model.log", "a")
        self.log.write(f"\n### {time.strftime('%Y-%m-%d %H:%M:%S')}  server: {shlex.join(cmd)}\n")
        self.log.flush()
        self.proc = subprocess.Popen(cmd, cwd=self.work, stdout=self.log, stderr=subprocess.STDOUT)
        (self.dir / "pid").write_text(f"{self.proc.pid}\n")

        def reap(proc=self.proc, d=self.dir):
            rc = proc.wait()
            if d.exists():
                (d / "exited").write_text(f"{rc}\n")
        threading.Thread(target=reap, daemon=True).start()
        os.environ[bm.SERVER_ENV] = str(self.dir)
        print(f"model server started (pid {self.proc.pid}): {shlex.join(cmd)}", flush=True)
        return self

    def __exit__(self, *exc):
        if self.proc is None:
            return False
        os.environ.pop(bm.SERVER_ENV, None)
        if self.proc.poll() is None:
            if bm.server_alive(self.dir):
                bm.server_request(self.dir, "stop\n", timeout=60)
            try:
                self.proc.wait(timeout=60)
            except subprocess.TimeoutExpired:
                self.proc.terminate()
                try:
                    self.proc.wait(timeout=20)
                except subprocess.TimeoutExpired:
                    self.proc.kill()
                    self.proc.wait()
        print(f"model server stopped (return code {self.proc.returncode})", flush=True)
        self.log.close()
        return False


def targets_path(cfg):
    t = bm.ROOT / cfg.get("targets", "data/processed/experimental_results.txt")
    if not t.exists():
        sys.exit(f"targets file not found: {t}")
    return t


# ------------------------------------------------------------------------------------------------
# fits
# ------------------------------------------------------------------------------------------------
def write_dakota_input(work, cfg, restart):
    free = cfg["fit"]["free"]
    names = list(free)
    txt = (TEMPLATES / "dakota.in").read_text().format(
        name=cfg["name"],
        restart="        read_restart = 'dakota.rst'\n" if restart else "",
        method=cfg["fit"].get("method", "coliny_direct"),
        max_evals=int(cfg["fit"]["max_evals"]),
        n_free=len(names),
        descriptors=" ".join(f"'{k}'" for k in names),
        lower=" ".join(repr(float(free[k][0])) for k in names),
        upper=" ".join(repr(float(free[k][1])) for k in names),
        driver=str(work / "driver.sh"),
        objective=cfg["fit"]["objective"],
        # numbers as `key = value`, keyword choices (strings) as `key choice`, e.g. division all_dimensions
        method_extra="".join(f"    {k} {v}\n" if isinstance(v, str) else f"    {k} = {v}\n"
                             for k, v in cfg["fit"].get("dakota_options", {}).items()),
    )
    (work / "dakota.in").write_text(txt)


def write_driver_sh(work, python):
    sh = work / "driver.sh"
    sh.write_text(f"#!/bin/bash\n# Dakota analysis driver (see runs/driver.py)\n"
                  f"exec {python} {bm.RUNS / 'driver.py'} \"$@\"\n")
    sh.chmod(0o755)


def archive_previous_attempt(work):
    """Before a restart: keep the Dakota outputs of the previous attempt (dakota.rst is reused)."""
    if not (work / "dakota.rst").exists():
        return False
    n = 1
    while (work / f"attempt{n}").exists():
        n += 1
    keep = work / f"attempt{n}"
    keep.mkdir()
    for f in ("dakota.in", "dakota.out", "dakota.err", "dakota.tab", "dakota_results.txt", "dakota.stdout"):
        if (work / f).exists():
            shutil.move(str(work / f), keep / f)
    shutil.copy2(work / "dakota.rst", keep / "dakota.rst")
    print(f"restart: previous Dakota outputs kept in {keep}", flush=True)
    return True


def n_rows(path):
    return sum(1 for _ in open(path)) - 1 if path.exists() else 0


# --- fits: (1) global DIRECT -> (2) local DIRECT around its best distinct points -> refinement ---
# A "point" is one model evaluation with a label ("global:412" = evaluation 412 of the global
# search, "local2:57" = evaluation 57 of local search 2), all 19 parameters and all errors.
def _point(label, r):
    return {"point": label, "params": {k: float(r[k]) for k in bm.PARAM_NAMES},
            "errors": {k: float(r[k]) for k in bm.ERROR_NAMES}}


def search_points(d):
    """Successful evaluations of one DIRECT search directory (global/ or local<i>/)."""
    if not (d / "evals.csv").exists():
        return []
    with open(d / "evals.csv") as f:
        return [_point(f"{d.name}:{r['eval_id']}", r) for r in csv.DictReader(f) if r["status"] == bm.OK]


def choose_starts(cfg, work):
    """Starts of the local searches: the `starts` best points of the global search that are
    distinct, i.e. at least 2 x radius apart (in units of each parameter's range) in some
    parameter, so that their local boxes do not overlap. Written once to starts.json."""
    f = work / "starts.json"
    if f.exists():
        return json.loads(f.read_text())
    obj, free = cfg["fit"]["objective"], cfg["fit"]["free"]
    sep = 2 * float(cfg["fit"]["radius"])
    pts = search_points(work / "global")
    best = []
    for p in sorted(pts, key=lambda p: p["errors"][obj]):
        if all(max(abs(p["params"][k] - q["params"][k]) / (hi - lo) for k, (lo, hi) in free.items()) >= sep
               for q in best):
            best.append(p)
        if len(best) == int(cfg["fit"]["starts"]):
            break
    out = {"objective": obj, "from": "global", "separation": sep, "n_global_ok": len(pts), "starts": best}
    f.write_text(json.dumps(out, indent=2) + "\n")
    return out


def local_box(cfg, start):
    """Local search box: start +- radius x (range of the search box), clipped to the search box."""
    r = float(cfg["fit"]["radius"])
    box = {}
    for k, (lo, hi) in cfg["fit"]["free"].items():
        x, h = start["params"][k], r * (hi - lo)
        box[k] = [max(lo, x - h), min(hi, x + h)]
    return box


def run_dakota(work, cfg, site):
    dakota = (site or {}).get("dakota") or shutil.which("dakota")
    if not dakota or not Path(dakota).exists():
        sys.exit("Dakota not found (hpc_site.toml: dakota = ..., or dakota on the PATH)")
    print(f"dakota: {dakota} -i dakota.in  (in {work}; output in dakota.stdout)", flush=True)
    with open(work / "dakota.stdout", "w") as out:
        subprocess.run([dakota, "-i", "dakota.in"], cwd=work, check=True, stdout=out,
                       stderr=subprocess.STDOUT)
    print(f"  {stop_reason(work) or 'no termination line in dakota.stdout'}", flush=True)
    (work / "dakota.done").write_text(f"{int(cfg['fit']['max_evals'])}\n")


def stop_reason(d):
    """Dakota/SCOLIB termination line, e.g. '[ Termination: Box-Size-Limit - ... ]' (printed to
    Dakota's stdout, not to dakota.out)."""
    f = d / "dakota.stdout"
    if not f.exists():
        return ""
    return next((l.strip() for l in f.read_text(errors="replace").splitlines() if "Termination" in l), "")


def local_done(lw, max_evals):
    m = lw / "dakota.done"
    return (m.exists() and int(m.read_text().split()[0]) == int(max_evals)) or \
        n_rows(lw / "evals.csv") >= int(max_evals)


def run_search(cfg, site, work, sub, box, max_evals, options, start, targets):
    """One DIRECT search (Dakota) over `box`, in work/<sub>/ (same driver and model binary as the
    fit). Stops at max_evals."""
    lw = work / sub
    lw.mkdir(exist_ok=True)
    lcfg = json.loads(json.dumps(cfg))
    lcfg["name"] = f"{cfg['name']}/{sub}"
    lcfg["fit"]["free"] = box
    lcfg["fit"]["max_evals"] = max_evals
    lcfg["fit"]["dakota_options"] = options
    lcfg["local_start"] = start
    if local_done(lw, max_evals):
        return
    restart = archive_previous_attempt(lw)
    bm.make_workdir(lw, lcfg, targets)
    write_driver_sh(lw, (site or {}).get("python", sys.executable))
    write_dakota_input(lw, lcfg, restart)
    if not (lw / "bin").exists():
        (lw / "bin").symlink_to("../bin")
    obj = cfg["fit"]["objective"]
    what = f"start {start['point']} ({obj} = {start['errors'][obj]:.4f})" if start else "whole search box"
    print(f"{cfg['name']}: DIRECT search {sub}, {what}, <= {max_evals} evaluations, {options}", flush=True)
    run_dakota(lw, lcfg, site)


def run_global(cfg, site, work, targets):
    """(1) Global DIRECT over the whole search box, in work/global/. Dakota options from
    [fit.global] (e.g. local_balance_parameter = SCOLIB max_boxsize_ratio: boxes smaller than
    this fraction of the largest box are not divided)."""
    n = int(cfg["fit"]["global"]["max_evals"])
    run_search(cfg, site, work, "global", cfg["fit"]["free"], n, direct_options(cfg), None, targets)


def direct_options(cfg):
    """Dakota DIRECT options, the same for the global and the local searches: `division` from [fit]
    and every [fit.global] entry except max_evals (e.g. local_balance_parameter)."""
    return {"division": cfg["fit"]["division"],
            **{k: v for k, v in cfg["fit"]["global"].items() if k != "max_evals"}}


def run_local(cfg, site, work, i, start, targets):
    """(2) Local search i: DIRECT in a box around start i, in work/local<i>/, with the same DIRECT
    options as the global search; stops at max_evals."""
    run_search(cfg, site, work, f"local{i}", local_box(cfg, start), int(cfg["fit"]["max_evals"]),
               direct_options(cfg), start, targets)


def local_dirs(work):
    return sorted((d for d in work.glob("local*") if d.is_dir()), key=lambda d: int(d.name[5:] or 0))


def local_points(work):
    return [p for d in local_dirs(work) for p in search_points(d)]


def candidates(cfg, work, per_search):
    """The `per_search` best distinct points (by their single evaluation) of EACH local search,
    so that the replicates compare the basins with each other as well as points within a basin."""
    if per_search <= 0:
        return []
    obj = cfg["fit"]["objective"]
    out, seen = [], set()
    for d in local_dirs(work):
        n = 0
        for p in sorted(search_points(d), key=lambda p: p["errors"][obj]):
            key = tuple(p["params"][k] for k in bm.PARAM_NAMES)
            if key not in seen:
                seen.add(key)
                out.append(p)
                n += 1
            if n == per_search:
                break
    return out


def _safe(label):
    return label.replace("/", "_").replace(":", "_").replace("#", "")


def refine(cfg, work, cands):
    """Re-evaluate the candidates `refine_reps` times each (independent seeds), so that the
    optimum is chosen by the mean error, not by one noisy evaluation (the minimum of many noisy
    evaluations is biased low). Restartable: refine.csv records finished replicates."""
    reps = int(cfg["fit"]["refine_reps"])
    if not cands or reps <= 0:
        return
    out_csv = work / "refine.csv"
    done = set()
    if out_csv.exists():
        with open(out_csv) as f:
            done = {(r["point"], r["rep"]) for r in csv.DictReader(f)}
    (work / "results" / "refine").mkdir(exist_ok=True)
    print(f"{cfg['name']}: refining the {len(cands)} best points x {reps} replicates", flush=True)
    obj = cfg["fit"]["objective"]
    # Replicates are independent: each has its own seed (hash of point and replicate number)
    for rank, c in enumerate(cands, 1):
        for rep in range(1, reps + 1):
            if (c["point"], str(rep)) in done:
                continue
            o = bm.run_model(work, work / "bin" / "bubble-model", c["params"], cfg["n_sims"],
                             cfg.get("domain", "lobe"), cfg["launcher"], work / "model.log",
                             bm.seed_for(cfg["name"], "refine", c["point"], rep))
            stop_if_crashed(o, f"{c['point']} rep {rep}")
            res = work / "results" / "results.txt"
            if o["status"] == bm.OK and res.exists():
                shutil.copyfile(res, work / "results" / "refine" / f"{_safe(c['point'])}_rep{rep}.txt")
            e = o["errors"] or {}
            bm.append_row(out_csv, ["rank", "point", "rep", *bm.PARAM_NAMES, *bm.ERROR_NAMES, "status",
                                    "seconds", "single", "seed"],
                          [rank, c["point"], rep, *(repr(c["params"][k]) for k in bm.PARAM_NAMES),
                           *(e.get(k, "") for k in bm.ERROR_NAMES), o["status"], o["seconds"],
                           repr(c["errors"][obj]), o["seed"]])
            val = f"{obj} = {e[obj]:.4f}" if e else o["status"]
            print(f"  point {rank} ({c['point']}, single {c['errors'][obj]:.4f}) rep {rep}: {val}", flush=True)


def stop_if_crashed(o, what):
    """A model run that did not complete (crash, or the launcher failed, e.g. srun at the job's
    time limit) is not recorded: stop here; re-running continues with this run."""
    if o["status"] == bm.CRASH:
        sys.exit(f"{what}: model run did not complete (return code {o['returncode']}, see model.log); "
                 f"not recorded, re-run to continue")


def _mean_sd(x):
    m = sum(x) / len(x)
    sd = (sum((v - m) ** 2 for v in x) / (len(x) - 1)) ** 0.5 if len(x) > 1 else float("nan")
    return m, sd


def write_best(work):
    """Optimum of a fit: the refined point with the minimum mean objective over its replicates
    (all replicates completed). Also records the starts, the local searches and the best single
    evaluation."""
    cfg = json.loads((work / "config.json").read_text())
    obj = cfg["fit"]["objective"]
    pts = search_points(work / "global") + local_points(work)
    if not pts:
        sys.exit(f"{work}: no successful evaluation in the DIRECT searches")
    single = min(pts, key=lambda p: p["errors"][obj])
    out = {"run": cfg["name"], "objective": obj,
           "best_single": {"point": single["point"], "value": single["errors"][obj]}}
    params = single["params"]
    if (work / "refine.csv").exists():
        with open(work / "refine.csv") as f:
            ref = list(csv.DictReader(f))
        by = {}
        for r in ref:
            by.setdefault(r["point"], []).append(r)
        cands = []
        for key, rr in by.items():
            good = [r for r in rr if r["status"] == bm.OK]
            vals = [float(r[obj]) for r in good]
            m, sd = _mean_sd(vals) if vals else (float("inf"), float("nan"))
            cands.append({"point": key, "rank": int(rr[0]["rank"]), "single": float(rr[0]["single"]),
                          "mean": m, "sd": sd, "n": len(vals), "n_runs": len(rr),
                          "errors_mean": {k: _mean_sd([float(r[k]) for r in good])[0] for k in bm.ERROR_NAMES}
                          if good else None,
                          "params": {k: float(rr[0][k]) for k in bm.PARAM_NAMES}})
        cands.sort(key=lambda c: c["mean"])
        full = [c for c in cands if c["n"] == c["n_runs"] and c["n"] > 0]
        if full:
            b = full[0]
            sds = [c["sd"] for c in full if c["n"] > 1]
            out["refine"] = {"per_search": cfg["fit"]["refine_per_search"], "reps": cfg["fit"]["refine_reps"],
                             "sigma_pooled": (sum(v * v for v in sds) / len(sds)) ** 0.5 if sds else None,
                             "candidates": cands}
            out.update({"point": b["point"], "value": b["mean"], "value_sd": b["sd"], "value_n": b["n"],
                        "errors": b["errors_mean"]})
            params = b["params"]
            rep1 = work / "results" / "refine" / f"{_safe(b['point'])}_rep1.txt"
            if rep1.exists():
                shutil.copyfile(rep1, work / "best_results.txt")
    if "refine" not in out:
        out.update({"point": single["point"], "value": single["errors"][obj], "errors": single["errors"]})
    def stage(d):
        lc = json.loads((d / "config.json").read_text())
        lp = search_points(d)
        st = lc.get("local_start")
        return {"dir": d.name, "start": st["point"] if st else None,
                "start_value": st["errors"][obj] if st else None, "box": lc["fit"]["free"],
                "dakota_options": lc["fit"].get("dakota_options", {}),
                "n_evals": n_rows(d / "evals.csv"), "n_ok": len(lp),
                "best_single": min((p["errors"][obj] for p in lp), default=None), "stop": stop_reason(d)}
    out["global"] = stage(work / "global")
    out["local"] = [stage(d) for d in local_dirs(work)]
    stages = [out["global"], *out["local"]]
    out["n_evals"] = sum(s["n_evals"] for s in stages)
    out["n_ok"] = sum(s["n_ok"] for s in stages)
    out["params"] = params
    (work / "best.json").write_text(json.dumps(out, indent=2) + "\n")
    bm.write_params(work / "best_params.txt", params)
    shutil.copyfile(work / "best_params.txt", work.parent / f"final_{cfg['name']}_params.txt")
    how = (f"mean of {out['value_n']} replicates, sd {out['value_sd']:.4f}; sigma_pooled "
           f"{out['refine']['sigma_pooled']:.4f}" if "refine" in out and out['refine']['sigma_pooled']
           else "single evaluation")
    print(f"{cfg['name']}: best {obj} = {out['value']:.4f} at {out['point']} ({how}); best single "
          f"evaluation {single['errors'][obj]:.4f} at {single['point']} ({out['n_ok']}/{out['n_evals']} ok)",
          flush=True)
    return out


def run_fit(cfg, site, work):
    """Fit: (1) global DIRECT over the search box; (2) a local DIRECT search around each of its
    `starts` best distinct points; (3) refinement of the best distinct points; best.json."""
    p = {**bm.base_params(cfg, work.parent), **cfg.get("params", {})}
    for k, (lo, hi) in cfg["fit"]["free"].items():
        p[k] = 0.5 * (lo + hi)          # placeholder; Dakota sets the free parameters
    bm.apply_tied(p, cfg.get("tied", {}), explicit=set(cfg.get("params", {})) | set(cfg["fit"]["free"]))
    cfg["resolved_params"] = p
    targets = targets_path(cfg)
    bm.make_workdir(work, cfg, targets)
    builds = [bm.build(work, cfg.get("defs", {}), cfg["mode"])]
    bm.write_manifest(work, cfg, targets, builds, cfg["mode"])
    with ModelServer(work, builds[0]["exe"], cfg):
        return fit_stages(cfg, site, work, targets)


def fit_stages(cfg, site, work, targets):
    run_global(cfg, site, work, targets)
    st = choose_starts(cfg, work)
    for i, s in enumerate(st["starts"], 1):
        run_local(cfg, site, work, i, s, targets)
    refine(cfg, work, candidates(cfg, work, int(cfg["fit"]["refine_per_search"])))
    return write_best(work)


# ------------------------------------------------------------------------------------------------
# case sets
# ------------------------------------------------------------------------------------------------
def subset(cases, k):
    """Test mode: k evenly spaced cases (always the first and the last)."""
    if len(cases) <= k:
        return cases
    idx = sorted({round(i * (len(cases) - 1) / (k - 1)) for i in range(k)}) if k > 1 else [0]
    return [cases[i] for i in idx]


def run_cases(cfg, site, work):
    cases = bm.generate_cases(cfg, work.parent)
    if cfg.get("test"):
        cases = subset(cases, cfg["max_cases"])
    targets = targets_path(cfg)
    bm.make_workdir(work, cfg, targets)
    cj = work / "cases.json"
    if cj.exists() and json.loads(cj.read_text()) != cases:
        sys.exit(f"{cj} exists and differs from the cases generated now (base changed?). "
                 f"Move the old work directory away to start again.")
    cj.write_text(json.dumps(cases, indent=1) + "\n")

    builds = {}
    for c in cases:
        key = bm.defs_string(c["defs"])
        if key not in builds:
            builds[key] = bm.build(work, c["defs"], cfg["mode"])
    bm.write_manifest(work, cfg, targets, list(builds.values()), cfg["mode"], {"n_cases": len(cases)})

    out_csv = work / "cases_out.csv"
    done = set()
    if out_csv.exists():
        with open(out_csv) as f:
            done = {int(r["id"]) for r in csv.DictReader(f)}
    (work / "results" / "cases").mkdir(exist_ok=True)
    todo = [c for c in cases if c["id"] not in done]
    print(f"{cfg['name']}: {len(cases)} cases, {len(done)} done, {len(todo)} to run", flush=True)
    # one build: one model server for all cases; several builds (thresholds): model started per case
    one = len(builds) == 1 and todo
    with ModelServer(work, next(iter(builds.values()))["exe"], cfg if one else {**cfg, "server": False}):
        run_case_list(cfg, work, builds, todo, out_csv)


def run_case_list(cfg, work, builds, todo, out_csv):
    for i, c in enumerate(todo, 1):
        exe = builds[bm.defs_string(c["defs"])]["exe"]
        o = bm.run_model(work, exe, c["params"], cfg["n_sims"], c["domain"], cfg["launcher"],
                         work / "model.log", bm.seed_for(cfg["name"], "case", c["id"]))
        stop_if_crashed(o, f"case {c['id']}")
        res = work / "results" / "results.txt"
        if o["status"] == bm.OK and res.exists():
            shutil.copyfile(res, work / "results" / "cases" / f"{c['id']:04d}_{c['name']}.txt")
        e = o["errors"] or {}
        bm.append_row(out_csv, ["id", "name", "domain", "defs", *bm.PARAM_NAMES, *bm.ERROR_NAMES,
                                "status", "seconds", "seed"],
                      [c["id"], c["name"], c["domain"], bm.defs_string(c["defs"]),
                       *(repr(c["params"][k]) for k in bm.PARAM_NAMES),
                       *(e.get(k, "") for k in bm.ERROR_NAMES), o["status"], o["seconds"], o["seed"]])
        ea = f"E_A = {e['E_A']:.3f}, E_B = {e['E_B']:.3f}" if e else ""
        print(f"  [{i}/{len(todo)}] #{c['id']} {c['name']}: {o['status']} {ea} ({o['seconds']} s)", flush=True)


# ------------------------------------------------------------------------------------------------
# commands
# ------------------------------------------------------------------------------------------------
def expand(name):
    """A group config (type = "group") stands for its runs, executed in order in one job."""
    c = bm.load_config(name)
    return [m for r in c["runs"] for m in expand(r)] if c["type"] == "group" else [name]


def all_runs(cfgs):
    """`submit all`: every config that is not part of a group (groups count as one run)."""
    in_groups = {m for n, c in cfgs.items() if c["type"] == "group" for m in expand(n)}
    return [n for n in cfgs if n not in in_groups]


def cmd_run(args):
    for name in expand(args.name):
        cfg, site = resolve(name, args.test, args.local)
        work = root_for(args.test) / name
        work.mkdir(parents=True, exist_ok=True)
        print(f"=== {name} ({cfg['type']}, {'test' if args.test else 'full'}, {cfg['mode']}) -> {work}", flush=True)
        if cfg["type"] == "fit":
            run_fit(cfg, site, work)
        else:
            run_cases(cfg, site, work)


def write_job(name, site, test, time):
    log_dir = root_for(test) / name
    log_dir.mkdir(parents=True, exist_ok=True)
    env = "\n".join(f"export {k}={shlex.quote(str(v))}" for k, v in site.get("env", {}).items())
    extra = "\n".join(f"#SBATCH {x}" for x in site.get("sbatch_extra", []))
    txt = (TEMPLATES / "job.sh").read_text().format(
        name=name, job_name=("t_" if test else "bm_") + name, account=site["account"],
        qos=site["qos_test"] if test else site["qos"], time=time,
        nodes=site["nodes"], ntasks_per_node=site["ntasks_per_node"], cpus_per_task=site["cpus_per_task"],
        log_dir=log_dir, extra_sbatch=extra,
        modules="\n".join(f"module load {m}" for m in site.get("modules", [])), env=env,
        root=bm.ROOT, python=site.get("python", "python3"), run_args=name + (" --test" if test else ""))
    job = log_dir / "job.sh"
    job.write_text(txt)
    return job


def cmd_submit(args):
    site = bm.load_site()
    if not args.dry_run and "dakota" in site and not Path(site["dakota"]).exists():
        sys.exit(f"Dakota not found at {site['dakota']} (hpc_site.toml)")
    cfgs = {n: bm.load_config(n) for n in bm.config_names()}
    for n in all_runs(cfgs) if args.name == "all" else [args.name]:
        times = site.get("job_time", {})
        time = args.time or (site.get("test_time", "01:00:00") if args.test
                             else times.get(n, times.get("default", "12:00:00")))
        job = write_job(n, site, args.test, time)
        cmd = ["sbatch", "--parsable", str(job)]
        if args.dry_run:
            print("dry-run:", shlex.join(cmd))
            continue
        jid = subprocess.run(cmd, check=True, capture_output=True, text=True).stdout.strip().split(";")[0]
        print(f"submitted {n}: job {jid}")


class Style:
    """ANSI styles, only when printing to a terminal (and NO_COLOR is not set)."""
    def __init__(self, on):
        c = (lambda code: f"\033[{code}m") if on else (lambda code: "")
        self.bold, self.dim, self.reset = c("1"), c("2"), c("0")
        self.green, self.yellow, self.red, self.cyan = c("32"), c("33"), c("31"), c("36")


def run_state(name, cfg, root, st):
    """(symbol, coloured state text) of one run, from its work directory."""
    w = root / name
    if cfg["type"] == "fit":
        obj = cfg["fit"]["objective"]
        if (w / "best.json").exists():
            b = json.loads((w / "best.json").read_text())
            if "refine" in b:
                val = f"{obj} = {b['value']:.3f} ± {b['value_sd']:.3f} (n = {b['value_n']})"
            else:
                val = f"{obj} = {b['value']:.3f}"
            bad = b["n_evals"] - b["n_ok"]
            extra = f"  {st.red}{bad} aborted{st.reset}" if bad else ""
            return f"{st.green}✔", f"{st.green}done{st.reset}  {st.bold}{val}{st.reset}{st.dim}  " \
                   f"({len(b['local'])} local searches, {b['n_evals']} evals){st.reset}{extra}"
        if (w / "starts.json").exists():
            pts = local_points(w)
            best = f"  best so far {min(p['errors'][obj] for p in pts):.3f}" if pts else ""
            nref = n_rows(w / "refine.csv")
            if nref:
                tot = len(local_dirs(w)) * int(cfg["fit"]["refine_per_search"]) * int(cfg["fit"]["refine_reps"])
                return f"{st.yellow}●", f"{st.yellow}refining{st.reset} {nref}/{tot} replicates{best}"
            n_st = len(json.loads((w / "starts.json").read_text())["starts"])
            locs = local_dirs(w)
            cur = f"{len(locs)}/{n_st} ({n_rows(locs[-1] / 'evals.csv')} evals)" if locs else f"0/{n_st}"
            return f"{st.yellow}●", f"{st.yellow}local search{st.reset} {cur}{best}"
        if (w / "global").exists():
            pts = search_points(w / "global")
            best = f"  best so far {min(p['errors'][obj] for p in pts):.3f}" if pts else ""
            n, tot = n_rows(w / "global" / "evals.csv"), int(cfg["fit"]["global"]["max_evals"])
            return f"{st.yellow}●", f"{st.yellow}global search{st.reset} {n}/{tot} evals{best}"
    elif (w / "cases.json").exists():
        tot = len(json.loads((w / "cases.json").read_text()))
        stat = []
        if (w / "cases_out.csv").exists():
            with open(w / "cases_out.csv") as f:
                stat = [r["status"] for r in csv.DictReader(f)]
        bad = sum(1 for x in stat if x != bm.OK)
        extra = f"  {st.red}{bad} aborted{st.reset}" if bad else ""
        if len(stat) >= tot:
            return f"{st.green}✔", f"{st.green}done{st.reset}  {tot}/{tot} cases{extra}"
        return f"{st.yellow}●", f"{st.yellow}running{st.reset}  {len(stat)}/{tot} cases{extra}"
    return f"{st.dim}○", f"{st.dim}not started{st.reset}"


def slurm_jobs():
    """{job name: 'RUNNING 3:12:05' | 'PENDING (Dependency)'} for the user's jobs, if squeue exists."""
    try:
        out = subprocess.run(["squeue", "-h", "-u", os.environ.get("USER", ""), "-o", "%j|%T|%M|%r"],
                             capture_output=True, text=True, check=True, timeout=20).stdout
    except (OSError, subprocess.SubprocessError):
        return None
    jobs = {}
    for line in out.splitlines():
        j, state, t, reason = (line.split("|") + ["", "", ""])[:4]
        jobs[j] = f"{state} {t}" if state == "RUNNING" else f"{state} ({reason})"
    return jobs


def cmd_list(args):
    st = Style(sys.stdout.isatty() and not os.environ.get("NO_COLOR"))
    root = root_for(args.test)
    cfgs = {n: bm.load_config(n) for n in bm.config_names()}
    if args.test:
        cfgs = {n: c if c["type"] == "group" else test_overrides(c) for n, c in cfgs.items()}
    jobs = slurm_jobs() if not args.test else None

    def job_text(name):
        if jobs is None:
            return ""
        j = jobs.get(("t_" if args.test else "bm_") + name)
        if not j:
            return f"{st.dim}no job{st.reset}"
        col = st.cyan if j.startswith("RUNNING") else st.dim
        return f"{col}{j}{st.reset}"

    def row(name, indent="", job=""):
        sym, text = run_state(name, cfgs[name], root, st)
        job = f"   {st.dim}│{st.reset} {job}" if job else ""
        print(f"  {indent}{sym}{st.reset} {name:<{14 - len(indent)}} {st.dim}{cfgs[name]['type']:<5}{st.reset}  {text}{job}")

    print(f"{st.bold}Runs in {root}{st.reset}")
    print(f"{st.dim}  ✔ done   ● in progress   ○ not started{st.reset}\n")
    print(f"{st.bold}submit all{st.reset}")
    for n in all_runs(cfgs):
        if cfgs[n]["type"] == "group":
            job = job_text(n)
            print(f"  {st.bold}{n}{st.reset} {st.dim}(one job: these runs in order){st.reset}"
                  + (f"   {st.dim}│{st.reset} {job}" if job else ""))
            for m in expand(n):   # a member submitted on its own (e.g. fitA) shows its own job
                own = jobs and jobs.get(("t_" if args.test else "bm_") + m)
                row(m, "  ", job=job_text(m) if own else "")
    for n in all_runs(cfgs):
        if cfgs[n]["type"] != "group":
            row(n, job=job_text(n))
    if args.long:
        print(f"\n{st.bold}descriptions{st.reset}")
        for n, c in cfgs.items():
            print(f"  {n:<15} {st.dim}{c.get('description', '')}{st.reset}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run", help="run a config now, in this shell (inside a Slurm job, or locally)")
    r.add_argument("name", help="config name (runs/configs/<name>.toml)")
    r.add_argument("--local", action="store_true", help="build with `make local`, no launcher, no Slurm")
    r.add_argument("--test", action="store_true", help="small version of the run, in runs/work/test/")
    s = sub.add_parser("submit", help="write the Slurm job of a config (or of `all`) and submit it")
    s.add_argument("name", help="config name, or `all` (every config not part of a group)")
    s.add_argument("--test", action="store_true",
                   help="small version of the run, in runs/work/test/, on the test QOS (qos_test, test_time)")
    s.add_argument("--time", help="Slurm time limit (HH:MM:SS) instead of [job_time] in hpc_site.toml")
    s.add_argument("--dry-run", action="store_true", help="write job.sh but do not submit it")
    ls = sub.add_parser("list", help="state of every run (and its Slurm job)")
    ls.add_argument("--test", action="store_true", help="show runs/work/test/")
    ls.add_argument("--long", action="store_true", help="also print each run's description")
    args = ap.parse_args()
    {"run": cmd_run, "submit": cmd_submit, "list": cmd_list}[args.cmd](args)


if __name__ == "__main__":
    main()
