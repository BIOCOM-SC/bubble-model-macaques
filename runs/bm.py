"""
Shared helpers for the Bubble Model run pipeline (configs, parameters, model runs, manifests).

Standard library only (runs on the cluster with a plain python3). Reading the TOML configs needs
Python >= 3.11 (tomllib) or the `tomli` package; the scripts that run inside a work directory
(driver.py) only read the resolved JSON config written by launch.py.
"""
import csv
import fcntl
import hashlib
import json
import os
import random
import shlex
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

RUNS = Path(__file__).resolve().parent          # bubble-model-macaques/runs
ROOT = RUNS.parent                              # bubble-model-macaques
MODEL = ROOT / "model"
LUNG = ROOT / "lung" / "macaque"
CONFIGS = RUNS / "configs"
WORK = RUNS / "work"                            # git-ignored: all run outputs

# Order of the 19 parameters in config/params_bm.txt (read by src/params.c with fixed fscanf formats)
PARAM_NAMES = ["v0", "fvimm", "r0", "fr0dev", "frimm", "frimmdev", "b", "nata", "f", "tmax",
               "dt", "n", "a", "rmin", "n0", "fndev", "tc", "ftcdev", "tmin"]

# Column order of results/error.txt (src/exp_err.c, calculate_and_write_errors).
# NB: the stdout line printed by the model labels nmi/nmd the wrong way round; error.txt is right.
ERROR_NAMES = ["nc", "nm", "nmi", "nmd", "vfc", "E_A", "E_B"]

DOMAIN_FLAGS = {"lobe": "-dlobe", "segment": "-dseg"}

# Status of a model run
OK, LIMIT, CRASH = "ok", "limit", "crash"

# Persistent model process (`bubble-model -server <dir>`, started by launch.py): directory of the
# running server, inherited by Dakota and driver.py. Unset: every evaluation starts the model.
SERVER_ENV = "BM_SERVER"
SERVER_TIMEOUT = 1800       # s without an answer -> the evaluation counts as a crash (MAX_SIMU_TIME = 600 s)


# ------------------------------------------------------------------------------------------------
# configuration
# ------------------------------------------------------------------------------------------------
def load_toml(path):
    try:
        import tomllib
    except ModuleNotFoundError:
        try:
            import tomli as tomllib
        except ModuleNotFoundError:
            sys.exit("Reading TOML needs Python >= 3.11 or the 'tomli' package "
                     "(on the cluster: load a recent python module).")
    with open(path, "rb") as f:
        return tomllib.load(f)


def config_names():
    return sorted(p.stem for p in CONFIGS.glob("*.toml") if not p.stem.startswith("_"))


def load_config(name):
    """Config `name` merged over configs/_base.toml (top-level keys; tables are merged)."""
    path = CONFIGS / f"{name}.toml"
    if not path.exists():
        sys.exit(f"No config {path}. Available: {', '.join(config_names())}")
    base = load_toml(CONFIGS / "_base.toml")
    cfg = load_toml(path)
    out = dict(base)
    for k, v in cfg.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = {**out[k], **v}
        else:
            out[k] = v
    out["name"] = name
    if out.get("type") not in ("fit", "cases", "group"):
        sys.exit(f"{name}: type must be 'fit', 'cases' or 'group'")
    return out


def load_site(required=True):
    """HPC site settings (account, queues, modules, Dakota path). Never part of the public folder."""
    cands = [os.environ.get("BM_HPC_SITE"), ROOT / "hpc_site.toml", ROOT.parent / "hpc_site.toml"]
    for c in cands:
        if c and Path(c).exists():
            site = load_toml(c)
            site["_path"] = str(Path(c).resolve())
            return site
    if required:
        sys.exit("hpc_site.toml not found (copy hpc_site.example.toml, or set BM_HPC_SITE)")
    return None


# ------------------------------------------------------------------------------------------------
# parameters
# ------------------------------------------------------------------------------------------------
def apply_tied(p, tied, explicit=()):
    """Tied parameters (e.g. frimm = fvimm) follow their source unless set explicitly."""
    for dst, src in tied.items():
        if dst not in explicit:
            p[dst] = p[src]
    return p


def write_params(path, p):
    missing = [k for k in PARAM_NAMES if k not in p]
    if missing:
        raise ValueError(f"missing parameters: {missing}")
    with open(path, "w") as f:
        f.write("\n".join(f"{k} = {float(p[k]):.12e}" for k in PARAM_NAMES))
        f.write("\n")


def read_dakota_params(path):
    """Dakota (standard format) parameters file -> ({descriptor: value}, eval_id)."""
    vals, eval_id = {}, None
    for line in Path(path).read_text().splitlines():
        tok = line.split()
        if len(tok) < 2:
            continue
        try:
            x = float(tok[0])
        except ValueError:
            continue
        if tok[1] == "eval_id":
            eval_id = int(x)
        elif tok[1] in PARAM_NAMES:
            vals[tok[1]] = x
    return vals, eval_id


def base_params(cfg, root, base=None):
    """Fixed parameters of a run: configs/_base.toml, or the optimum of a fit ("fitA", ...)."""
    base = base or cfg.get("base", "default")
    if base == "default":
        p = dict(load_toml(CONFIGS / "_base.toml")["defaults"])
    else:
        best = Path(root) / base / "best.json"
        if not best.exists():
            sys.exit(f"{cfg['name']}: base '{base}' needs {best} (run the fit first)")
        b = json.loads(best.read_text())
        if load_config(base).get("fit", {}).get("local") and "local" not in b:
            sys.exit(f"{cfg['name']}: {best} is from the global stage only (the fit has a local "
                     f"stage: run '{base}' again to finish it)")
        p = dict(b["params"])
    return p


# ------------------------------------------------------------------------------------------------
# model build and run
# ------------------------------------------------------------------------------------------------
def defs_string(defs):
    return " ".join(f"-D{k}={v}" for k, v in sorted(defs.items()))


def binary_name(defs):
    if not defs:
        return "bubble-model"
    tag = hashlib.sha1(defs_string(defs).encode()).hexdigest()[:8]
    return f"bubble-model_{tag}"


def build(workdir, defs, mode):
    """Compile the model into <workdir>/bin/ (make local | make hpc) with compile-time defs."""
    exe = Path(workdir) / "bin" / binary_name(defs)
    exe.parent.mkdir(exist_ok=True)
    target = "hpc" if mode == "hpc" else "local"
    cmd = ["make", "-s", "-C", str(MODEL), target, f"EXEC={exe}", f"DEFS={defs_string(defs)}"]
    print("build:", shlex.join(cmd), flush=True)
    subprocess.run(cmd, check=True)
    return {"exe": str(exe), "defs": defs, "command": shlex.join(cmd), "sha256": sha256(exe)}


def model_command(exe, launcher, n_sims, domain, seed):
    return shlex.split(launcher) + [str(exe), "-noplot", DOMAIN_FLAGS[domain], "-N", str(int(n_sims)),
                                    "-seed", str(int(seed))]


def seed_for(*parts):
    """Seed of one evaluation: 64-bit hash of its identity, e.g. (run name, Dakota eval_id),
    (run, "refine", point, replicate), (run, "case", case id). Distinct evaluations get independent
    seeds; re-evaluating the same one (restart) gets the same seed, hence the same result."""
    return int.from_bytes(hashlib.sha256("|".join(map(str, parts)).encode()).digest()[:8], "little")


def read_errors(workdir):
    """results/error.txt -> (status, {name: value}). The model writes '100, 100, ...' (integers)
    when a run hits MAX_LESIONS / MAX_* / MAX_SIMU_TIME, and nothing if it crashed."""
    f = Path(workdir) / "results" / "error.txt"
    if not f.exists():
        return CRASH, None
    lines = [l for l in f.read_text().splitlines() if l.strip()]
    if not lines:
        return CRASH, None
    tok = [t.strip() for t in lines[-1].split(",")]
    if len(tok) != len(ERROR_NAMES):
        return CRASH, None
    if all(t == "100" for t in tok):
        return LIMIT, None
    return OK, dict(zip(ERROR_NAMES, map(float, tok)))


def run_model(workdir, exe, params, n_sims, domain, launcher, log, seed):
    """Write params, run one model evaluation (n_sims simulations, given seed) in workdir, return
    outcome. Uses the running model server for `exe` if there is one (launch.py), otherwise starts
    the model (through the launcher)."""
    workdir = Path(workdir)
    write_params(workdir / "config" / "params_bm.txt", params)
    err = workdir / "results" / "error.txt"
    if err.exists():
        err.unlink()
    server = server_for(exe)
    t0 = time.time()
    with open(log, "a") as lf:
        if server:
            lf.write(f"\n### {time.strftime('%Y-%m-%d %H:%M:%S')}  server {server}: -N {int(n_sims)} "
                     f"-seed {int(seed)} {DOMAIN_FLAGS[domain]} (model output in the server's log)\n")
            lf.flush()
            answer = server_request(server, f"workdir {workdir.resolve()}\nn {int(n_sims)}\n"
                                            f"seed {int(seed)}\ndomain {domain}\n")
            rc = 0 if answer == "ok" else 1
            if rc:
                lf.write(f"### server: {answer}\n")
        else:
            cmd = model_command(exe, launcher, n_sims, domain, seed)
            lf.write(f"\n### {time.strftime('%Y-%m-%d %H:%M:%S')}  {shlex.join(cmd)}\n")
            lf.flush()
            rc = subprocess.run(cmd, cwd=workdir, stdout=lf, stderr=subprocess.STDOUT).returncode
    status, errors = read_errors(workdir)
    if rc != 0 and status == OK:
        status = CRASH
    return {"status": status, "errors": errors, "seconds": round(time.time() - t0, 2), "returncode": rc,
            "seed": int(seed)}


# ------------------------------------------------------------------------------------------------
# persistent model process (client side; the server is started and stopped by launch.py)
# Files in the server directory: exe (binary it runs), pid (launcher process), exited (written by
# launch.py when the process ends), lock, req.<id> / done.<id> (protocol: see src/main.c).
# ------------------------------------------------------------------------------------------------
def server_for(exe):
    """Directory of the running server for this binary, or None."""
    d = os.environ.get(SERVER_ENV)
    if not d:
        return None
    d = Path(d)
    try:
        return d if Path((d / "exe").read_text().strip()) == Path(exe).resolve() else None
    except OSError:
        return None


def server_alive(d):
    if (d / "exited").exists():
        return False
    try:
        os.kill(int((d / "pid").read_text()), 0)
    except (OSError, ValueError):
        return False
    return True


def server_request(d, text, timeout=SERVER_TIMEOUT):
    """Send one request (one evaluation at a time: lock file), wait for the answer: 'ok', 'error',
    'stopped', or 'server died' / 'no answer after ... s' if it never comes."""
    with open(d / "lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rid = f"{os.getpid()}_{time.time_ns()}"
        tmp = d / f"req.{rid}.tmp"
        tmp.write_text(text)
        tmp.rename(d / f"req.{rid}")
        done = d / f"done.{rid}"
        t0 = last = time.time()
        while not done.exists():
            time.sleep(0.002)
            now = time.time()
            if now - last > 0.5:
                last = now
                if not server_alive(d):
                    return "server died"
                if now - t0 > timeout:
                    return f"no answer after {timeout} s"
        answer = done.read_text().strip()
        done.unlink()
        return answer


def append_row(path, header, row):
    """Append one row to a CSV record, writing the header first. Refuses a file with other
    columns (e.g. written by an older version of the pipeline)."""
    path = Path(path)
    new = not path.exists()
    if not new:
        with open(path, newline="") as f:
            first = next(csv.reader(f), None)
        if first != header:
            sys.exit(f"{path}: columns differ from this version of the pipeline (no 'seed' column? "
                     f"written before the random-number change); move the run directory away")
    with open(path, "a", newline="") as f:
        w = csv.writer(f)
        if new:
            w.writerow(header)
        w.writerow(row)


# ------------------------------------------------------------------------------------------------
# work directory and manifest
# ------------------------------------------------------------------------------------------------
def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def tree_hash(files):
    h = hashlib.sha256()
    for f in sorted(files):
        h.update(str(f.relative_to(ROOT)).encode())
        h.update(sha256(f).encode())
    return h.hexdigest()


def git_info():
    try:
        rev = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True,
                             text=True, check=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain", "--", "model", "runs",
                                "data", "lung"], capture_output=True, text=True).stdout.strip()
        return {"commit": rev, "dirty": bool(dirty)}
    except (OSError, subprocess.CalledProcessError):
        return {"commit": None, "dirty": None}


def make_workdir(workdir, cfg, targets):
    """Self-contained work dir: config/{params_bm.txt, experimental_results.txt, macaque -> lung},
    results/, bin/. The model is run with workdir as cwd (it uses relative paths)."""
    workdir = Path(workdir)
    (workdir / "config").mkdir(parents=True, exist_ok=True)
    (workdir / "results").mkdir(exist_ok=True)
    link = workdir / "config" / "macaque"
    if link.is_symlink() and link.resolve() != LUNG.resolve():
        link.unlink()                       # stale link (work dir moved or copied)
    if not link.exists():
        link.symlink_to(os.path.relpath(LUNG, link.parent))
    shutil.copyfile(targets, workdir / "config" / "experimental_results.txt")
    (workdir / "config.json").write_text(json.dumps(cfg, indent=2) + "\n")


def write_manifest(workdir, cfg, targets, builds, mode, extra=None):
    src = [*MODEL.glob("src/*.c"), *MODEL.glob("include/*.h"), MODEL / "Makefile"]
    man = {
        "run": cfg["name"],
        "created": time.strftime("%Y-%m-%d %H:%M:%S %z"),
        "host": socket.gethostname(),
        "mode": mode,
        "git": git_info(),
        "model_source_sha256": tree_hash(src),
        "targets": {"file": str(Path(targets).resolve().relative_to(ROOT)), "sha256": sha256(targets)},
        "builds": builds,
        "python": sys.version.split()[0],
        "env": {k: os.environ[k] for k in ("OMP_NUM_THREADS", "SLURM_JOB_ID", "SLURM_NNODES",
                                           "SLURM_NTASKS", "SLURM_CPUS_PER_TASK") if k in os.environ},
        "config": cfg,
    }
    man.update(extra or {})
    path = Path(workdir) / "manifest.json"
    # keep the manifests of earlier attempts (restarts)
    if path.exists():
        n = len(list(Path(workdir).glob("manifest.*.json"))) + 1
        path.rename(Path(workdir) / f"manifest.{n}.json")
    path.write_text(json.dumps(man, indent=2) + "\n")


# ------------------------------------------------------------------------------------------------
# case generators (type = "cases")
# ------------------------------------------------------------------------------------------------
def linspace(a, b, n):
    n = int(n)
    return [a] if n == 1 else [a + (b - a) * i / (n - 1) for i in range(n)]


def _value(v, name, base):
    return base[name] if v == "base" else float(v)


def _spec_values(spec, name, base):
    """Values of one swept parameter: [a, b, n] linspace (a/b may be "base"), {rel = [lo, hi, n]}
    relative to the base value, or {values = [...]}."""
    if isinstance(spec, dict) and "linspace" in spec:
        spec = spec["linspace"]
    if isinstance(spec, list):
        a, b, n = spec
        return linspace(_value(a, name, base), _value(b, name, base), n)
    if "rel" in spec:
        lo, hi, n = spec["rel"]
        return [base[name] * f for f in linspace(lo, hi, n)]
    return [float(x) for x in spec["values"]]


def _lhs(n, k, rng):
    """Latin hypercube in [0,1]^k: one point per stratum and dimension, random within stratum."""
    cols = []
    for _ in range(k):
        perm = list(range(n))
        rng.shuffle(perm)
        cols.append([(perm[i] + rng.random()) / n for i in range(n)])
    return [[cols[j][i] for j in range(k)] for i in range(n)]


def generate_cases(cfg, root):
    """Expand the generators of a cases config into a list of cases
    {id, name, params (all 19), domain, defs}."""
    tied = cfg.get("tied", {})
    domain0 = cfg.get("domain", "lobe")
    fixed = cfg.get("params", {})
    cases = []

    def add(name, over, base_name=None, domain=None, defs=None):
        base = base_params(cfg, root, base_name)
        p = {**base, **fixed}
        explicit = set(fixed)
        for k, v in over.items():
            p[k] = _value(v, k, base) if isinstance(v, str) else float(v)
            explicit.add(k)
        apply_tied(p, tied, explicit)
        # No immune response (f_imm = 1): optionally move t_imm past t_max, so that the maximum
        # radii of existing lesions are not redrawn at t_imm (activate_immune, src/lesions.c).
        if "tc_no_immunity" in cfg and p["fvimm"] == 1 and p["frimm"] == 1:
            p["tc"] = float(cfg["tc_no_immunity"])
        cases.append({"id": len(cases) + 1, "name": name, "params": p,
                      "domain": domain or domain0, "defs": {**cfg.get("defs", {}), **(defs or {})}})

    for c in cfg.get("case", []):
        add(c["name"], c.get("params", {}), c.get("base"), c.get("domain"), c.get("defs"))

    for s in cfg.get("sweep", []):
        base = base_params(cfg, root, s.get("base"))
        for i, x in enumerate(_spec_values(s, s["param"], base)):
            add(f"{s['name']}_{i + 1:02d}", {s["param"]: x}, s.get("base"))

    for g in cfg.get("grid", []):
        base = base_params(cfg, root, g.get("base"))
        axes = [(k, _spec_values(v, k, base)) for k, v in g.get("params", {}).items()]
        daxes = [(k, [float(x) for x in v]) for k, v in g.get("defs", {}).items()]
        combos = [({}, {})]
        for k, vals in axes:
            combos = [({**pp, k: x}, dd) for pp, dd in combos for x in vals]
        for k, vals in daxes:
            combos = [(pp, {**dd, k: x}) for pp, dd in combos for x in vals]
        for i, (pp, dd) in enumerate(combos):
            add(f"{g['name']}_{i + 1:03d}", pp, g.get("base"), defs=dd)

    if "lhs" in cfg:
        # rel = [lo, hi]: each parameter uniform in [lo, hi] x its base value
        L = cfg["lhs"]
        base = base_params(cfg, root, L.get("base"))
        rng = random.Random(int(L["seed"]))
        names = L["params"]
        lo, hi = L["rel"]
        box = {k: [base[k] * lo, base[k] * hi] for k in names}
        for i, u in enumerate(_lhs(int(L["n"]), len(names), rng)):
            over = {k: box[k][0] + (box[k][1] - box[k][0]) * x for k, x in zip(names, u)}
            add(f"lhs_{i + 1:04d}", over, L.get("base"))

    return cases
