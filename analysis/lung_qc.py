"""
Quality checks of the virtual macaque lung: where the bronchial-tree terminals lie relative
to the lobe meshes the model uses, and the distributions of the terminal quantities the model uses.

    python analysis/lung_qc.py      -> analysis/data/lung_terminals.dat, analysis/tables/lung_qc.txt

Inside test: ray casting, as the model (geometry.c is_inside_lobe: one ray along +x, odd number of
crossings = inside). Also with two more ray directions (majority vote of 3) to check that the
fixed ray does not misclassify points.
Model constants (model/include/main.h): lung volume 319 000 mm^3, trachea diameter 8.1 mm;
radius of the volume irrigated by a terminal R_i = (3 V / 4 pi)^(1/3) d_i / d_trachea (aux.c).
Terminals outside their own lobe: distance to the own lobe mesh (point-triangle), and how the model's
placement (lesions.c: terminal drawn by p_initial_infection, position uniform in the ball of radius R_i,
both redrawn until inside the terminal's lobe) treats them: fraction of the ball inside the lobe
(acceptance, Monte Carlo with the model's +x ray test) and the resulting share of initial lesions.
dlobe in terminals.txt (distance to the nearest vertex of the original lobe surfaces, unsigned) is a
model input (cap on the lesion radius) and is reported as such.
"""
import argparse
import math
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
LUNG = HERE.parent / "lung" / "macaque"
DATA = HERE / "data"
TABLES = HERE / "tables"
V_LUNG, D_TRACHEA = 319000.0, 8.1
LOBES = ["RUL", "RML", "RLL", "LIC", "LUL", "LML", "LLL"]  # files L1..L7 (p_init shares = week-3 lobe frequencies)
K_MC, SEED_MC = 200, 7    # Monte Carlo points per terminal for the acceptance


def read_mesh(path):
    """Lobe/segment file: N, then per face 3 vertices + the first vertex repeated (closing)."""
    nums = path.read_text().split()
    n = int(nums[0])
    v = np.array(nums[1:], dtype=float).reshape(n, 4, 3)
    return v[:, :3, :]


def crossings(points, tri, ray):
    """Number of crossings of the ray from each point with the triangles (Moller-Trumbore, t > eps)."""
    eps = 1e-9
    v0, v1, v2 = tri[:, 0], tri[:, 1], tri[:, 2]
    e1, e2 = v1 - v0, v2 - v0
    h = np.cross(ray, e2)                      # (F,3)
    a = np.einsum("ij,ij->i", e1, h)           # (F,)
    ok = np.abs(a) > eps
    e1, e2, v0, h, a = e1[ok], e2[ok], v0[ok], h[ok], a[ok]
    f = 1.0 / a
    out = np.zeros(len(points), dtype=int)
    for s in range(0, len(points), 512):
        P = points[s:s + 512]
        S = P[:, None, :] - v0[None, :, :]                      # (P,F,3)
        u = f * np.einsum("pfj,fj->pf", S, h)
        q = np.cross(S, e1[None, :, :])
        w = f * (q @ ray)
        t = f * np.einsum("pfj,fj->pf", q, e2)
        hit = (u >= 0) & (u <= 1) & (w >= 0) & (u + w <= 1) & (t > eps)
        out[s:s + 512] = hit.sum(axis=1)
    return out


def distance_to_mesh(points, tri):
    """Minimum distance from each point to a triangle mesh (closest point on each triangle)."""
    a, b, c = tri[:, 0], tri[:, 1], tri[:, 2]
    ab, ac = b - a, c - a
    n = np.cross(ab, ac)
    nn = np.einsum("ij,ij->i", n, n)
    d00, d01, d11 = (np.einsum("ij,ij->i", ab, ab), np.einsum("ij,ij->i", ab, ac), np.einsum("ij,ij->i", ac, ac))
    den = d00 * d11 - d01 * d01
    out = np.empty(len(points))
    for i, p in enumerate(points):
        ap = p - a
        q = p - (np.einsum("ij,ij->i", ap, n) / nn)[:, None] * n      # projection on the plane
        aq = q - a
        d20, d21 = np.einsum("ij,ij->i", aq, ab), np.einsum("ij,ij->i", aq, ac)
        v, w = (d11 * d20 - d01 * d21) / den, (d00 * d21 - d01 * d20) / den
        d = np.where((v >= 0) & (w >= 0) & (v + w <= 1), np.linalg.norm(p - q, axis=1), np.inf)
        for s, e in ((a, b), (b, c), (c, a)):                          # closest point on the edges
            se = e - s
            t = np.clip(np.einsum("ij,ij->i", p - s, se) / np.einsum("ij,ij->i", se, se), 0, 1)
            d = np.minimum(d, np.linalg.norm(p - (s + t[:, None] * se), axis=1))
        out[i] = d.min()
    return out


def acceptance(pos, Ri, lobe, meshes):
    """Fraction of the ball of radius R_i around each terminal inside its own lobe (model test, +x ray)."""
    rng = np.random.default_rng(SEED_MC)
    u = rng.normal(size=(len(pos), K_MC, 3))
    u /= np.linalg.norm(u, axis=2, keepdims=True)
    pts = pos[:, None, :] + u * (Ri[:, None] * np.cbrt(rng.random((len(pos), K_MC))))[:, :, None]
    acc = np.zeros(len(pos))
    for k, m in enumerate(meshes):
        sel = lobe == k + 1
        hit = crossings(pts[sel].reshape(-1, 3), m, np.array([1.0, 0, 0])) % 2 == 1
        acc[sel] = hit.reshape(-1, K_MC).mean(axis=1)
    return acc


def main():
    argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter).parse_args()
    TABLES.mkdir(exist_ok=True)
    rows = LUNG.joinpath("terminals.txt").read_text().split()
    n = int(rows[0])
    T = np.array([r.split(",") for r in rows[1:n + 1]], dtype=float)
    pos, seg, dseg, lobe, dlobe, diam, pinit = T[:, :3], T[:, 3].astype(int), T[:, 4], T[:, 5].astype(int), \
        T[:, 6], T[:, 7], T[:, 8]
    meshes = [read_mesh(LUNG / "lobes" / f"L{k}.txt") for k in range(1, 8)]
    rays = [np.array([1.0, 0, 0]), np.array([0.267, 0.534, 0.802]), np.array([-0.577, 0.577, -0.577])]
    inside = np.zeros((n, 7, 3), dtype=bool)
    for k, m in enumerate(meshes):
        for r, ray in enumerate(rays):
            inside[:, k, r] = crossings(pos, m, ray / np.linalg.norm(ray)) % 2 == 1
    model_own = inside[np.arange(n), lobe - 1, 0]             # model's test (+x ray), own lobe
    vote = inside.sum(axis=2) >= 2                            # (n,7) majority of 3 rays
    vote_own = vote[np.arange(n), lobe - 1]
    in_any = vote.any(axis=1)
    disagree = (inside[:, :, 0] != vote).any(axis=1).sum()
    w = pinit / pinit.sum()
    Ri = (3 * V_LUNG / (4 * math.pi)) ** (1 / 3) * diam / D_TRACHEA

    out_own = ~model_own
    dmesh = np.zeros(n)                                       # distance to own lobe mesh (outside terminals)
    for k, m in enumerate(meshes):
        sel = out_own & (lobe == k + 1)
        dmesh[sel] = distance_to_mesh(pos[sel], m)
    acc = acceptance(pos, Ri, lobe, meshes)
    seeded = w * acc / (w * acc).sum()                       # share of initial lesions seeded by each terminal
    do = dmesh[out_own]
    L = ["Virtual macaque lung: terminals of the bronchial tree vs lobe meshes",
         f"terminals: {n}; lobe meshes: " + ", ".join(f"L{k + 1} {len(m)} faces" for k, m in enumerate(meshes)),
         f"outside own lobe (model test, +x ray): {out_own.sum()} ({100 * out_own.mean():.1f} %)",
         f"outside own lobe (majority of 3 rays): {(~vote_own).sum()} ({100 * (~vote_own).mean():.1f} %); "
         f"terminals where the +x ray disagrees with the majority for some lobe: {disagree}",
         f"outside all lobes (majority): {(~in_any).sum()} ({100 * (~in_any).mean():.1f} %)",
         f"outside own lobe but inside another lobe: {((~vote_own) & in_any).sum()}",
         f"share of initial-infection probability on terminals outside their own lobe: {100 * w[out_own].sum():.1f} %",
         f"dlobe of terminals outside their own lobe: median "
         f"{np.median(dlobe[out_own]):.2f} mm, 90th percentile {np.percentile(dlobe[out_own], 90):.2f} mm, max {dlobe[out_own].max():.2f} mm",
         f"distance to own lobe mesh of terminals outside their own lobe: median {np.median(do):.2f} mm, "
         f"within 1 mm {100 * (do < 1).mean():.0f} %, within 2 mm {100 * (do < 2).mean():.0f} %, "
         f"90th percentile {np.percentile(do, 90):.2f} mm, max {do.max():.2f} mm; closer than their R_i: "
         f"{100 * (do < Ri[out_own]).mean():.0f} %",
         f"acceptance of a lesion position drawn around the terminal (fraction of its R_i ball inside its lobe, "
         f"{K_MC} points per terminal): terminals inside their lobe mean {acc[~out_own].mean():.2f}, "
         f"outside mean {acc[out_own].mean():.2f} ({(acc[out_own] == 0).sum()} of {out_own.sum()} never accepted)",
         f"share of initial lesions seeded by terminals outside their own lobe: {100 * seeded[out_own].sum():.1f} % "
         f"(p_initial_infection {100 * w[out_own].sum():.1f} %)",
         "dlobe (terminals.txt; distance to the nearest vertex of the original lobe surfaces, unsigned) is the "
         "model's cap on the lesion radius (-dlobe), not the distance to the lobe meshes above.",
         "Handling in the model: initial lesions are placed at a terminal drawn by p_initial_infection, displaced"
         " uniformly within R_i, and re-drawn until the position is inside the terminal's lobe (lesions.c);"
         " reinfections are re-drawn up to 10 times, then dropped (reinfection.c; not counted).",
         "",
         f"terminal diameter d_i (mm): median {np.median(diam):.3f}, IQR {np.percentile(diam, 25):.3f}–{np.percentile(diam, 75):.3f}, "
         f"range {diam.min():.3f}–{diam.max():.3f}",
         f"irrigated radius R_i (mm): median {np.median(Ri):.2f}, IQR {np.percentile(Ri, 25):.2f}–{np.percentile(Ri, 75):.2f}, "
         f"range {Ri.min():.2f}–{Ri.max():.2f}",
         f"distance to segment wall (dseg, mm): median {np.median(dseg):.2f}, IQR {np.percentile(dseg, 25):.2f}–{np.percentile(dseg, 75):.2f}",
         f"distance to lobe wall (dlobe, mm): median {np.median(dlobe):.2f}, IQR {np.percentile(dlobe, 25):.2f}–{np.percentile(dlobe, 75):.2f}",
         "per lobe (file L1..L7): terminals, outside own (model test), share of p_init, share of initial lesions seeded:"]
    for k in range(7):
        sel = lobe == k + 1
        L.append(f"  L{k + 1} {LOBES[k]}: {sel.sum():5d} terminals, {out_own[sel].sum():4d} outside ({100 * out_own[sel].mean():.1f} %), "
                 f"p_init {100 * w[sel].sum():.1f} %, seeded {100 * seeded[sel].sum():.1f} %")
    with open(DATA / "lung_lobes.dat", "w") as f:
        f.write("# i lobe terminals outside_own(model) frac_outside p_init_share frac_outside_all frac_in_other_lobe\n")
        for k in range(7):
            sel = lobe == k + 1
            f.write(f"{k} {LOBES[k]} {sel.sum()} {out_own[sel].sum()} {out_own[sel].mean():.5f} {w[sel].sum():.5f} "
                    f"{(out_own & ~in_any)[sel].mean():.5f} {(out_own & in_any)[sel].mean():.5f}\n")
    (TABLES / "lung_qc.txt").write_text("\n".join(L) + "\n")
    print("\n".join(L))
    with open(DATA / "lung_terminals.dat", "w") as f:
        f.write("# x y z seg dseg lobe dlobe diameter p_init R_i inside_own(model) inside_own(vote) inside_any(vote) "
                "dist_own_lobe_mesh(outside only, else 0) acceptance\n")
        for i in range(n):
            f.write(" ".join(f"{v:.6g}" for v in T[i]) + f" {Ri[i]:.5f} {int(model_own[i])} {int(vote_own[i])} {int(in_any[i])}"
                    f" {dmesh[i]:.4f} {acc[i]:.3f}\n")


if __name__ == "__main__":
    main()
