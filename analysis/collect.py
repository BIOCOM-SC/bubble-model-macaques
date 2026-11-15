"""
Collect the run outputs (runs/work/<run>/) into the plain-text tables that the gnuplot scripts in
analysis/plots/ read (analysis/data/) and the summary tables quoted in the paper (analysis/tables/).

    python analysis/collect.py            # everything
    python analysis/collect.py fits prcc  # selected parts (vaccinated needs immgrid first)

Parts: exp animals fits minerr experiments immgrid vaccinated small sweeps prcc tradeoff
(see analysis/README.md for the figure or table each one feeds).

Conventions
- Model curves: results.txt row k holds the state at the end of simulated day k+1 (the fit error
  uses rows 20/52/76 = days 21/53/77), so the day column written here is k + 1.
- Every curve file has the same columns as data/exp.dat:
    day  NC sd  NM sd  NMD sd  NMI sd  VFC sd  MAC sd
  (mean and SD over the simulated animals).
- Missing values (aborted cases) are written as NaN.
"""
import argparse
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
WORK = ROOT / "runs" / "work"
DATA = HERE / "data"
TABLES = HERE / "tables"
TARGETS = ROOT / "data" / "processed" / "experimental_results.txt"
PER_ANIMAL = ROOT / "data" / "processed" / "per_animal.csv"
CONSOL = ROOT / "data" / "processed" / "consolidations.csv"
TARGETS_RIDV = ROOT / "data" / "processed" / "ridv" / "experimental_results.txt"   # process.py --cohort RIDV
DAYS_RIDV = [21, 49, 77]

DAYS = [21, 53, 77]
METRICS = ["nc", "nm", "nmd", "nmi", "vfc", "mac"]
# 0-based columns of results.txt (mean; SD is the next column)
RES_COL = {"nc": 3, "nm": 5, "vfc": 7, "nmd": 9, "nmi": 11, "mac": 15}
# row order of experimental_results.txt (mean row, SD row per metric)
EXP_ROW = {"nc": 0, "nm": 2, "nmd": 4, "nmi": 6, "vfc": 8, "mac": 10}
FITS = ["fitA", "fitB", "noimm", "timm21"]
FREE_ORDER = ["v0", "fvimm", "r0", "b", "nata", "tc", "n0"]
PARAM_TEX = {"v0": "v_0", "fvimm": "f_imm", "r0": "R_max0", "b": "beta", "nata": "rho",
             "tc": "t_imm", "n0": "N_0", "rmin": "r_min"}
ERRS = ["nc", "nm", "nmi", "nmd", "vfc", "E_A", "E_B"]


def fmt(x, nd=6):
    return "NaN" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}g}"


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def read_results(path):
    """results.txt -> array (tmax, 13): day, mean/sd for the six metrics in METRICS order."""
    raw = np.loadtxt(path, delimiter=",", usecols=range(17))
    out = np.empty((raw.shape[0], 13))
    out[:, 0] = raw[:, 0] + 1
    for i, m in enumerate(METRICS):
        out[:, 1 + 2 * i] = raw[:, RES_COL[m]]
        out[:, 2 + 2 * i] = raw[:, RES_COL[m] + 1]
    return out


HEADER = "# day  NC sd  NM sd  NMD sd  NMI sd  VFC sd  MAC sd"


def write_curve(arr, path, comment=""):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        if comment:
            f.write(f"# {comment}\n")
        f.write(HEADER + "\n")
        for row in arr:
            f.write(" ".join(fmt(float(v), 7) for v in row) + "\n")


def nan_curve(tmax=80):
    a = np.full((tmax, 13), np.nan)
    a[:, 0] = np.arange(1, tmax + 1)
    return a


# ------------------------------------------------------------------------------------------------
def read_targets(path, days):
    """experimental_results.txt -> array (3, 13) in the curve-file layout (day, mean/sd x 6)."""
    rows = [[float(x) for x in line.split(",")] for line in path.read_text().split("\n") if line.strip()]
    arr = np.empty((3, 13))
    arr[:, 0] = days
    for i, m in enumerate(METRICS):
        arr[:, 1 + 2 * i] = rows[EXP_ROW[m]]
        arr[:, 2 + 2 * i] = rows[EXP_ROW[m] + 1]
    return arr


def experimental():
    arr = read_targets(TARGETS, DAYS)
    write_curve(arr, DATA / "exp.dat", "experimental targets, mean and sample SD over n = 10 animals")
    return arr


def per_animal():
    """Per-animal trajectories (paired data): one block per animal, gnuplot `index`."""
    rows = read_csv(PER_ANIMAL)
    cons = read_csv(CONSOL)
    mac = {}
    for r in cons:
        if r["is_consolidation"] == "True" and r["present"] == "True":
            mac.setdefault((r["animal"], int(r["day"])), []).append(float(r["main_axis"]))
    animals = sorted({r["animal"] for r in rows})
    with open(DATA / "per_animal.dat", "w") as f:
        f.write("# one block per animal (gnuplot index i = animal i), columns: " + HEADER[2:] +
                " (sd columns = 0; MAC = mean main axis of that animal's consolidations, NaN if none)\n")
        for a in animals:
            f.write(f"# {a}\n")
            for r in sorted((r for r in rows if r["animal"] == a), key=lambda r: int(r["day"])):
                d = int(r["day"])
                m = mac.get((a, d))
                vals = [r["NC"], r["NM"], r["NMD"], r["NMI"], r["VFC"], np.mean(m) if m else float("nan")]
                f.write(f"{d} " + " ".join(f"{float(v):.6g} 0" for v in vals) + "\n")
            f.write("\n\n")
    (DATA / "per_animal_names.txt").write_text(" ".join(animals) + "\n")
    print(f"per_animal: {len(animals)} animals")


# ------------------------------------------------------------------------------------------------
def at_bound(v, lo, hi, tol=0.02):
    r = hi - lo
    if v <= lo + tol * r:
        return "lower"
    if v >= hi - tol * r:
        return "upper"
    return ""


def fits():
    """Fit table (parameters, objective mean ± sd, per-metric errors, sigma_E, bounds hit) and the
    model curves at each optimum (best_results.txt = replicate 1 of the optimum, 10^4 animals)."""
    table = []
    for fit in FITS:
        w = WORK / fit
        best = json.loads((w / "best.json").read_text())
        cfg = json.loads((w / "config.json").read_text())
        free = cfg["fit"]["free"]
        p = best["params"]
        bounds = {k: at_bound(p[k], *free[k]) for k in free}
        table.append({"fit": fit, "objective": best["objective"], "value": best["value"],
                      "sd": best["value_sd"], "n": best["value_n"],
                      "sigma_pooled": best["refine"]["sigma_pooled"],
                      "best_single": best["best_single"]["value"], "point": best["point"],
                      "n_evals": best["n_evals"], "n_ok": best["n_ok"],
                      "errors": best["errors"], "params": p, "free": free, "bounds": bounds})
        write_curve(read_results(w / "best_results.txt"), DATA / "curves" / f"{fit}.dat",
                    f"{fit} optimum ({best['point']}), replicate 1, 10^4 simulated animals")
    # values used in figure labels/markers (loaded by the gnuplot scripts; nothing typed by hand)
    by = {t["fit"]: t for t in table}
    pa = by["fitA"]["params"]
    with open(DATA / "fit_values.gp", "w") as f:
        f.write("# written by collect.py (fits): optima used in figure labels and markers\n")
        f.write(f"E0 = {by['fitA']['value']:.6g}\nSIGE = {by['fitA']['sigma_pooled']:.6g}\n")
        f.write(f"TFIT = {pa['tc']:.6g}; FFIT = {pa['fvimm']:.6g}\n")
        f.write(f"N0F = {pa['n0']:.6g}; V0F = {pa['v0']:.6g}; R0F = {pa['r0']:.6g}; BFIT = {pa['b']:.6g}\n")
        for k in FITS:
            f.write(f"EOPT_{k} = {by[k]['value']:.6g}\n")
    # readable table
    L = []
    L.append("Fits (refined optimum = lowest mean over 10 replicates of 10^4 simulated animals)")
    L.append("Search box (all fits): " + ", ".join(
        f"{PARAM_TEX[k]} [{v[0]:g}, {v[1]:g}]" for k, v in table[0]["free"].items()))
    L.append("")
    head = f"{'':22s}" + "".join(f"{t['fit']:>12s}" for t in table)
    L.append(head)
    L.append(f"{'objective':22s}" + "".join(f"{t['objective']:>12s}" for t in table))
    L.append(f"{'objective mean':22s}" + "".join(f"{t['value']:12.3f}" for t in table))
    L.append(f"{'  sd (10 reps)':22s}" + "".join(f"{t['sd']:12.3f}" for t in table))
    L.append(f"{'  sigma_E pooled':22s}" + "".join(f"{t['sigma_pooled']:12.4f}" for t in table))
    L.append(f"{'  best single eval':22s}" + "".join(f"{t['best_single']:12.3f}" for t in table))
    L.append(f"{'evaluations (ok)':22s}" + "".join(f"{t['n_evals']:>7d}({t['n_ok']})".rjust(12) for t in table))
    L.append("errors at the optimum (mean of 10 reps):")
    for e in ERRS:
        L.append(f"{'  ' + e:22s}" + "".join(f"{t['errors'][e]:12.3f}" for t in table))
    L.append("parameters:")
    for k in FREE_ORDER:
        cells = []
        for t in table:
            v = t["params"][k]
            b = t["bounds"].get(k)
            mark = "" if b is None else ("*" if b else "")
            if b is None:
                mark = " (fix)"
            cells.append(f"{v:.4g}{mark}".rjust(12))
        L.append(f"{'  ' + PARAM_TEX[k]:22s}" + "".join(cells))
    L.append("* = within 2 % of the search range from a bound; (fix) = not free in this fit.")
    L.append("Bounds hit: " + "; ".join(
        f"{t['fit']}: " + (", ".join(f"{PARAM_TEX[k]} {b}" for k, b in t["bounds"].items() if b) or "none")
        for t in table))
    (TABLES / "fits.txt").write_text("\n".join(L) + "\n")
    # Table 1 of the paper (row order and 3 decimals)
    T1 = [("n0", "N_0", ""), ("v0", "v_0", "day^-1"), ("r0", "R_max0", "mm"), ("nata", "rho", "day^-1 mm^-1"),
          ("b", "beta", "mm^-1"), ("tc", "t_imm", "day"), ("fvimm", "f_imm", "")]
    T = ["Table 1 (new values). Columns: range min, range max, fit A, fit B | no immunity, t_imm = 21 d (suppl.)"]
    for k, name, unit in T1:
        lo, hi = table[0]["free"][k]
        cells = [f"{t['params'][k]:.3f}" if k in t["free"] else
                 ("—" if (t["fit"] == "noimm" and k == "tc") else f"{t['params'][k]:g} (fixed)") for t in table]
        T.append(f"{name:7s} {unit:13s} {lo:g} & {hi:g} & {cells[0]} & {cells[1]}   | {cells[2]} & {cells[3]}")
    T.append("Objective at the optimum (mean ± SD of 10 replicates): " + ", ".join(
        f"{t['fit']} {t['objective']} {t['value']:.3f} ± {t['sd']:.3f}" for t in table))
    (TABLES / "table1.txt").write_text("\n".join(T) + "\n")
    print("\n".join(T))
    print("\n".join(L))
    return table


# ------------------------------------------------------------------------------------------------
def best_so_far(vals):
    out, b = [], math.inf
    for v in vals:
        b = min(b, v)
        out.append(b)
    return out


def minerr():
    """Best-so-far objective of the global search and of each local search, per fit."""
    d = DATA / "minerr"
    d.mkdir(parents=True, exist_ok=True)
    summary = []
    for fit in FITS:
        w = WORK / fit
        best = json.loads((w / "best.json").read_text())
        rows = read_csv(w / "global" / "evals.csv")
        vals = [float(r["objective"]) if r["status"] == "ok" else math.inf for r in rows]
        bsf = best_so_far(vals)
        n_g = len(vals)
        with open(d / f"{fit}_global.dat", "w") as f:
            f.write("# eval  objective  best_so_far\n")
            for i, (v, b) in enumerate(zip(vals, bsf), 1):
                f.write(f"{i} {fmt(v)} {fmt(b)}\n")
        marks = {n: bsf[n - 1] for n in (500, 1000, 2000, 3000, 4000, 5000) if n <= n_g}
        loc = []
        for i in range(1, 6):
            rows = read_csv(w / f"local{i}" / "evals.csv")
            vals = [float(r["objective"]) if r["status"] == "ok" else math.inf for r in rows]
            bsf_l = best_so_far(vals)
            with open(d / f"{fit}_local{i}.dat", "w") as f:
                f.write(f"# local search {i}; x = global evals ({n_g}) + local eval index\n")
                f.write("# eval  objective  best_so_far\n")
                for j, (v, b) in enumerate(zip(vals, bsf_l), 1):
                    f.write(f"{n_g + j} {fmt(v)} {fmt(b)}\n")
            loc.append(bsf_l[-1])
        (d / f"{fit}_refined.dat").write_text(f"# refined optimum mean sd\n{best['value']} {best['value_sd']}\n")
        summary.append(f"{fit}: global {n_g} evals, best-so-far " +
                       ", ".join(f"@{n} {v:.3f}" for n, v in marks.items()) +
                       f"; global best {bsf[-1]:.3f}; local bests " + ", ".join(f"{v:.3f}" for v in loc) +
                       f"; refined {best['value']:.3f} ± {best['value_sd']:.3f}")
    (TABLES / "minerr.txt").write_text("\n".join(summary) + "\n")
    print("\n".join(summary))


# ------------------------------------------------------------------------------------------------
def config_case_names(run):
    """Case names the current runs/configs/<run>.toml generates (cases dropped from a config are
    ignored even if an older work dir still holds them)."""
    sys.path.insert(0, str(ROOT / "runs"))
    import bm
    return {c["name"] for c in bm.generate_cases(bm.load_config(run), WORK)}


def grid_position(run):
    """(t_imm, f_imm) of each grid case by name, as on the grid axes: the f_imm = 1 cases run with
    t_imm past t_max (tc_no_immunity), but keep their place on the t_imm axis."""
    sys.path.insert(0, str(ROOT / "runs"))
    import bm
    cfg = bm.load_config(run)
    cfg.pop("tc_no_immunity", None)
    return {c["name"]: (float(c["params"]["tc"]), round(float(c["params"]["fvimm"]), 6))
            for c in bm.generate_cases(cfg, WORK)}


def cases(run):
    """Cases of a `cases` run that are in its current config -> data/curves/<run>/<name>.dat;
    returns the cases_out rows."""
    w = WORK / run
    names = config_case_names(run)
    rows = [r for r in read_csv(w / "cases_out.csv") if r["name"] in names]
    for r in rows:
        p = w / "results" / "cases" / f"{int(r['id']):04d}_{r['name']}.txt"
        arr = read_results(p) if (r["status"] == "ok" and p.exists()) else nan_curve()
        write_curve(arr, DATA / "curves" / run / f"{r['name']}.dat", f"{run}/{r['name']} status {r['status']}")
    return rows


def write_family(tag, run, rows, par):
    """gnuplot include data/fam_<tag>.gp: arrays FF (curve files) and FV (parameter values)."""
    n = len(rows)
    with open(DATA / f"fam_{tag}.gp", "w") as f:
        f.write(f"# family {tag}: {run}, parameter {par}\nNFAM = {n}\n")
        f.write(f"array FF[{n}] = [" + ", ".join(f'"curves/{run}/{r["name"]}.dat"' for r in rows) + "]\n")
        f.write(f"array FV[{n}] = [" + ", ".join(f"{float(r[par]):.6g}" for r in rows) + "]\n")


def value_at(run, name, metric, day):
    a = np.loadtxt(DATA / "curves" / run / f"{name}.dat")
    i = METRICS.index(metric)
    return float(a[day - 1, 1 + 2 * i]), float(a[day - 1, 2 + 2 * i])


def experiments():
    rows = cases("experiments")
    for prefix, par in (("fimm", "fvimm"),):
        sel = [r for r in rows if r["name"].startswith(prefix + "_")]
        write_family(prefix, "experiments", sel, par)
    L = ["Experiments at the fit-A optimum (10^4 simulated animals each); values at day 77 (mean ± SD)"]
    for name in ("control", "immunosuppressed", "dseg"):
        r = next(r for r in rows if r["name"] == name)
        L.append(f"{name:18s} E_A {float(r['E_A']):8.3f}  E_B {float(r['E_B']):8.3f}  " + "  ".join(
            f"{m.upper()} {v:.4g}±{s:.3g}" for m, (v, s) in ((m, value_at('experiments', name, m, 77)) for m in METRICS)))
    for day in (21, 53, 77):
        c = value_at("experiments", "control", "vfc", day)[0]
        L.append(f"VFC day {day}: control {c:.4f}, immunosuppressed {value_at('experiments', 'immunosuppressed', 'vfc', day)[0]:.4f}, "
                 f"dseg {value_at('experiments', 'dseg', 'vfc', day)[0]:.4f}")
    c, im, ds = (value_at("experiments", n, "vfc", 77)[0] for n in ("control", "immunosuppressed", "dseg"))
    L.append(f"Immunosuppressed / control VFC at day 77: {im / c:.2f}x; control is {100 * (1 - c / im):.1f} % lower")
    L.append(f"Segment-limited / control VFC at day 77: {ds / c:.3f} ({100 * (ds / c - 1):+.1f} %)")
    sel = [r for r in rows if r["name"].startswith("fimm_")]
    L.append("f_imm sweep (t_imm fixed): " + ", ".join(
        f"{float(r['fvimm']):.3f}:{value_at('experiments', r['name'], 'vfc', 77)[0]:.4f}" for r in sel))
    (TABLES / "experiments.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def immgrid():
    """t_imm x f_imm grid: day-77 values of each metric, gnuplot grid format (blank line per tc)."""
    rows = cases("immgrid")
    pos = grid_position("immgrid")
    tcs = sorted({pos[r["name"]][0] for r in rows})
    fs = sorted({pos[r["name"]][1] for r in rows})
    by = {pos[r["name"]]: r for r in rows}
    with open(DATA / "immgrid.dat", "w") as f:
        f.write("# t_imm f_imm E_A E_B " + " ".join(f"{m}77 sd" for m in METRICS) +
                " status (NaN = case aborted: one simulated animal reached a model safety limit, MAX_* in main.h)\n")
        for tc in tcs:
            for fv in fs:
                r = by[(tc, fv)]
                if r["status"] == "ok":
                    vals = " ".join(f"{fmt(v)} {fmt(s)}" for v, s in (value_at("immgrid", r["name"], m, 77) for m in METRICS))
                    ea, eb = r["E_A"], r["E_B"]
                else:
                    vals = " ".join("NaN NaN" for _ in METRICS)
                    ea = eb = "NaN"
                f.write(f"{tc:.6g} {fv:.6g} {ea} {eb} {vals} {r['status']}\n")
            f.write("\n")
    bad = [(r["name"], float(r["tc"]), float(r["fvimm"])) for r in rows if r["status"] != "ok"]
    # summary: VFC at day 77 along lines of the grid
    fit = json.loads((WORK / "fitA" / "best.json").read_text())["params"]
    L = [f"immgrid: {len(tcs)} t_imm x {len(fs)} f_imm, t_imm {tcs[0]:g}-{tcs[-1]:g}, f_imm {fs[0]:g}-{fs[-1]:g}; "
         f"aborted: {bad}", f"fit-A point: t_imm {fit['tc']:.2f}, f_imm {fit['fvimm']:.3f}"]
    g = np.loadtxt(DATA / "immgrid.dat", usecols=range(16))
    vfc = g[:, 4 + 2 * METRICS.index("vfc")]
    L.append("VFC77 table (rows t_imm, cols f_imm):")
    L.append("t_imm\\f " + " ".join(f"{x:7.2f}" for x in fs))
    for i, tc in enumerate(tcs):
        L.append(f"{tc:7.1f} " + " ".join(f"{v:7.4f}" for v in vfc[i * len(fs):(i + 1) * len(fs)]))
    (TABLES / "immgrid.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L[:2]))


def error_A(curve, exp):
    """E_A of model/src/exp_err.c for one curve (day, mean/sd x 6) against targets `exp` (layout of
    read_targets): squared z-scores of NC, NM and VFC at the three scan days, weights 1, 1, 5."""
    w = np.array([1.0, 1.0, 5.0])
    rows = exp[:, 0].astype(int) - 1
    e = 0.0
    for m in ("nc", "nm", "vfc"):
        i = 1 + 2 * METRICS.index(m)
        e += float((w * ((curve[rows, i] - exp[:, i]) / exp[:, i + 1]) ** 2).sum())
    return e


def vaccinated():
    """Intradermal-BCG-vaccinated macaques (Suppl. Section S5): E_A of every cell of the immgrid run
    (curves written by immgrid()) against the unvaccinated and the vaccinated targets; best cell when
    only f_imm, only t_imm, or both may change (the other parameters at fit A)."""
    exp_ridv = read_targets(TARGETS_RIDV, DAYS_RIDV)
    write_curve(exp_ridv, DATA / "exp_ridv.dat",
                "intradermal-BCG-vaccinated macaques (n = 5), mean and sample SD; scans at days 21, 49, 77")
    exp_ruv = read_targets(TARGETS, DAYS)
    names = config_case_names("immgrid")
    rows = [r for r in read_csv(WORK / "immgrid" / "cases_out.csv") if r["name"] in names]
    pos = grid_position("immgrid")
    rows.sort(key=lambda r: pos[r["name"]])
    tc = np.array([pos[r["name"]][0] for r in rows])
    fv = np.array([pos[r["name"]][1] for r in rows])
    E = {"ruv": [], "ridv": []}
    for r in rows:
        c = np.loadtxt(DATA / "curves" / "immgrid" / f"{r['name']}.dat")
        ok = r["status"] == "ok"
        E["ruv"].append(error_A(c, exp_ruv) if ok else np.nan)
        E["ridv"].append(error_A(c, exp_ridv) if ok else np.nan)
    E = {k: np.array(v) for k, v in E.items()}
    with open(DATA / "vaccinated_err.dat", "w") as f:
        f.write("# immgrid cells: t_imm f_imm E_A(unvaccinated targets) E_A(intradermal-BCG targets); NaN = aborted\n")
        for i in range(len(rows)):
            if i and tc[i] != tc[i - 1]:
                f.write("\n")
            f.write(f"{tc[i]:g} {fv[i]:g} {fmt(E['ruv'][i])} {fmt(E['ridv'][i])}\n")
    fit = json.loads((WORK / "fitA" / "best.json").read_text())["params"]
    ut, uf = np.unique(tc), np.unique(fv)
    t0, f0 = ut[np.argmin(abs(ut - fit["tc"]))], uf[np.argmin(abs(uf - fit["fvimm"]))]
    hyp = {"F": tc == t0, "T": fv == f0, "B": np.ones(len(rows), bool)}   # f_imm only, t_imm only, both
    L = [f"Vaccinated comparison on the immgrid run ({len(ut)} t_imm x {len(uf)} f_imm), other parameters at fit A",
         f"Cell nearest fit A: t_imm {t0:g}, f_imm {f0:.2f}", ""]
    with open(DATA / "vaccinated.gp", "w") as f:
        f.write("# best cells (lowest E_A) per hypothesis against the intradermal-BCG targets: curve file,\n"
                "# t_imm, f_imm, E_A; and the lowest E_A of each panel of the error grid\n")
        for k in ("ruv", "ridv"):
            j = int(np.nanargmin(E[k]))
            f.write(f"EMIN_{k.upper()} = {E[k][j]:.6g}\nEMAX_{k.upper()} = {np.nanmax(E[k]):.6g}\n")
            L.append(f"{k.upper()}: lowest E_A {E[k][j]:.3f} at t_imm {tc[j]:g}, f_imm {fv[j]:.2f}; "
                     f"E_A at the fit-A cell {E[k][(tc == t0) & (fv == f0)][0]:.3f}")
        for h, sel in hyp.items():
            idx = np.where(sel)[0]
            j = idx[np.nanargmin(E["ridv"][idx])]
            f.write(f'RIDV_{h}_FILE = "curves/immgrid/{rows[j]["name"]}.dat"; RIDV_{h}_T = {tc[j]:g}; '
                    f'RIDV_{h}_F = {fv[j]:g}; RIDV_{h}_E = {E["ridv"][j]:.6g}\n')
            L.append(f"RIDV, {dict(F='only f_imm free', T='only t_imm free', B='both free')[h]}: "
                     f"t_imm {tc[j]:g}, f_imm {fv[j]:.2f}, E_A {E['ridv'][j]:.3f}")
    (TABLES / "vaccinated.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def small_runs():
    """immsplit, thresholds: curves + summary tables."""
    L = []
    for run in ("immsplit",):
        rows = cases(run)
        L.append(f"{run}:")
        for r in rows:
            vals = "  ".join(f"{m.upper()}77 {value_at(run, r['name'], m, 77)[0]:.4g}" for m in METRICS) \
                if r["status"] == "ok" else r["status"]
            L.append(f"  {r['name']:12s} E_A {float(r['E_A']):8.3f} E_B {float(r['E_B']):8.3f}  {vals}")
    rows = cases("thresholds")
    L.append("thresholds (micronodule min MA x consolidation min MA, mm), errors at the fit-A optimum:")
    L.append(f"  {'micro':>6s} {'consol':>6s} " + " ".join(f"{e:>8s}" for e in ERRS))
    for r in rows:
        defs = dict(x.lstrip("-D").split("=") for x in r["defs"].split())
        mi, co = float(defs["MICRONOD_MIN_MAINAXIS"]), float(defs["CONS_MIN_MAINAXIS"])
        L.append(f"  {mi:6.1f} {co:6.1f} " + " ".join(f"{float(r[e]):8.3f}" for e in ERRS))
    (TABLES / "small_runs.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def sweeps():
    """Single-parameter sweeps (+-25 %, 20 values): curves + objective vs parameter."""
    rows = cases("sweeps")
    fitA = json.loads((WORK / "fitA" / "best.json").read_text())
    p0, sig = fitA["params"], fitA["refine"]["sigma_pooled"]
    E0 = fitA["value"]
    L = [f"Sweeps around fit A (E_A* = {E0:.4f}, sigma_E = {sig:.4f}); range of each parameter with "
         f"E_A <= E_A* + 2 sigma_E = {E0 + 2 * sig:.4f} (single evaluations, 10^4 animals each):"]
    with open(DATA / "sweeps_err.dat", "w") as f:
        f.write(f"# E_A* {E0} sigma_E {sig}\n# one block per parameter (index): rel value " + " ".join(ERRS) + "\n")
        for par in ["v0", "fvimm", "r0", "b", "nata", "rmin", "n0", "tc"]:
            sel = [r for r in rows if r["name"].startswith(f"sweep_{par}_")]
            write_family(f"sweep_{par}", "sweeps", sel, par)
            f.write(f"# {par}\n")
            ok_rel = []
            for r in sel:
                v = float(r[par])
                rel = v / p0[par]
                f.write(f"{rel:.6g} {v:.6g} " + " ".join(r[e] for e in ERRS) + "\n")
                if float(r["E_A"]) <= E0 + 2 * sig:
                    ok_rel.append(rel)
            f.write("\n\n")
            ea = [float(r["E_A"]) for r in sel]
            L.append(f"  {PARAM_TEX[par]:7s} (opt {p0[par]:.4g}): E_A from {min(ea):.3f} to {max(ea):.3f} over ±25 %; "
                     + (f"within 2σ: rel {min(ok_rel):.3f}–{max(ok_rel):.3f}" if ok_rel else "no sweep point within 2σ"))
    (TABLES / "sweeps.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


# ------------------------------------------------------------------------------------------------
def prcc_matrix(X, Y):
    """PRCC of each column of Y on each column of X (rank-transform, regress out other ranks)."""
    n, P = X.shape
    RX = np.apply_along_axis(stats.rankdata, 0, X)
    RY = np.apply_along_axis(stats.rankdata, 0, Y)
    out = np.zeros((P, Y.shape[1]))
    for i in range(P):
        A = np.column_stack([np.ones(n), np.delete(RX, i, axis=1)])
        rx = RX[:, i] - A @ np.linalg.lstsq(A, RX[:, i], rcond=None)[0]
        ry = RY - A @ np.linalg.lstsq(A, RY, rcond=None)[0]
        den = np.sqrt((rx ** 2).sum() * (ry ** 2).sum(axis=0))
        with np.errstate(invalid="ignore", divide="ignore"):
            out[i] = np.where(den > 0, (rx[:, None] * ry).sum(axis=0) / den, 0.0)
    return out


def prcc():
    rows = cases("prcc")
    pars = ["v0", "fvimm", "r0", "b", "nata", "rmin", "n0", "tc"]
    ok = [r for r in rows if r["status"] == "ok"]
    X = np.array([[float(r[p]) for p in pars] for r in ok])
    curves = [np.loadtxt(DATA / "curves" / "prcc" / f"{r['name']}.dat") for r in ok]
    n, P = X.shape
    dof = n - 2 - (P - 1)
    tcrit = stats.t.ppf(0.975, dof)
    rcrit = tcrit / math.sqrt(tcrit ** 2 + dof)
    # Bonferroni over 8 parameters x 6 outputs (final-day figure)
    tb = stats.t.ppf(1 - 0.05 / (2 * P * len(METRICS)), dof)
    rcrit_b = tb / math.sqrt(tb ** 2 + dof)
    days = curves[0][:, 0].astype(int)
    with open(DATA / "prcc_daily.dat", "w") as f:
        f.write(f"# PRCC per day; n = {n} LHS points (uniform ±25 % around fit A), {P} parameters; "
                f"|PRCC| < {rcrit:.4f} not significant (alpha 0.05, two-sided)\n")
        f.write("# one block per output (index 0..5 = " + " ".join(METRICS) + "); columns: day " + " ".join(pars) + "\n")
        final = {}
        for k, m in enumerate(METRICS):
            Y = np.array([c[:, 1 + 2 * k] for c in curves]).T  # (days, n)
            R = prcc_matrix(X, Y.T)  # (P, days)
            final[m] = R[:, list(days).index(77)]
            f.write(f"# {m}\n")
            for j, d in enumerate(days):
                f.write(f"{d} " + " ".join(f"{R[i, j]:.5f}" for i in range(P)) + "\n")
            f.write("\n\n")
    se = lambda r: math.sqrt((1 - r * r) / dof)
    with open(DATA / "prcc_final.dat", "w") as f:
        f.write(f"# PRCC at day 77; rcrit {rcrit:.5f} (alpha 0.05), Bonferroni(48) {rcrit_b:.5f}\n")
        f.write("# one block per output (index 0..5 = " + " ".join(METRICS) + "); columns: i param prcc se\n")
        for m in METRICS:
            for i, p in enumerate(pars):
                f.write(f"{i} {p} {final[m][i]:.5f} {se(final[m][i]):.5f}\n")
            f.write("\n\n")
    (DATA / "prcc_rcrit.dat").write_text(f"{rcrit:.6f} {rcrit_b:.6f} {n}\n")
    L = [f"PRCC: n = {n} of {len(rows)} LHS points ok; rcrit(0.05) = {rcrit:.4f}; Bonferroni(48) = {rcrit_b:.4f}",
         "PRCC at day 77 (rows outputs, cols " + " ".join(pars) + "):"]
    for m in METRICS:
        L.append(f"  {m:4s} " + " ".join(f"{x:+.3f}" for x in final[m]))
    (TABLES / "prcc.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


# ------------------------------------------------------------------------------------------------
def all_evals(fit):
    rows = read_csv(WORK / fit / "global" / "evals.csv")
    for i in range(1, 6):
        rows += read_csv(WORK / fit / f"local{i}" / "evals.csv")
    return [r for r in rows if r["status"] == "ok"]


def tradeoff():
    """Near-optimal parameter sets of fit A: all evaluations (Suppl. trade-off figure) and the refined
    candidates within 2 sigma_E of the optimum (figure + Suppl. table of equivalent sets)."""
    best = json.loads((WORK / "fitA" / "best.json").read_text())
    E0, sig = best["value"], best["refine"]["sigma_pooled"]
    with open(DATA / "fitA_evals.dat", "w") as f:
        f.write("# all ok evaluations of fit A (global + locals): " + " ".join(FREE_ORDER) + " E_A\n")
        for r in sorted(all_evals("fitA"), key=lambda r: -float(r["E_A"])):  # best last (drawn on top)
            f.write(" ".join(r[k] for k in FREE_ORDER) + f" {r['E_A']}\n")
    c = [x for x in best["refine"]["candidates"] if x["n"] == x["n_runs"]]
    near = sorted((x for x in c if x["mean"] <= E0 + 2 * sig), key=lambda x: x["params"]["v0"])
    with open(DATA / "fitA_refined_near.dat", "w") as f:
        f.write("# fit A refined candidates with mean E_A (10 reps) <= E_A* + 2 sigma_E: "
                + " ".join(FREE_ORDER) + " mean_E_A\n")
        for x in near:
            f.write(" ".join(repr(x["params"][k]) for k in FREE_ORDER) + f" {x['mean']}\n")
    L = [f"Fit A: refined optimum E_A* = {E0:.4f}, sigma_E = {sig:.4f}",
         f"{len(near)} of {len(c)} refined candidates with mean E_A <= E_A* + 2 sigma_E = {E0 + 2 * sig:.4f}"
         " (sorted by v_0):",
         "  " + " ".join(f"{PARAM_TEX[k]:>8s}" for k in FREE_ORDER) + f" {'mean E_A':>9s}"]
    L += ["  " + " ".join(f"{x['params'][k]:8.4g}" for k in FREE_ORDER) + f" {x['mean']:9.4f}" for x in near]
    (TABLES / "tradeoff.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))


PARTS = {"exp": experimental, "animals": per_animal, "fits": fits, "minerr": minerr,
         "experiments": experiments, "immgrid": immgrid, "vaccinated": vaccinated, "small": small_runs, "sweeps": sweeps,
         "prcc": prcc, "tradeoff": tradeoff}

if __name__ == "__main__":
    DATA.mkdir(exist_ok=True)
    TABLES.mkdir(exist_ok=True)
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("parts", nargs="*", metavar="part", help=f"default: all, in this order: {' '.join(PARTS)}")
    todo = ap.parse_args().parts or list(PARTS)
    if set(todo) - set(PARTS):
        ap.error(f"unknown part(s) {sorted(set(todo) - set(PARTS))}; parts: {' '.join(PARTS)}")
    for k in todo:
        print(f"--- {k}")
        PARTS[k]()
