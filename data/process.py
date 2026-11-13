#!/usr/bin/env python3
"""
Process the public CT raw data (Nogueira et al. 2023, Pathogens 12:236,
Supplementary File S2) into the experimental targets used by the Bubble Model.
VFC uses the lung volume of each animal at CT1 (raw/lung_volumes.csv).

Cohorts (--cohort):
  RUV   unvaccinated rhesus macaques (RUV01-RUV11), excluding RUV05 (miliary TB): the fit targets.
        Time points: CT1, CT2, CT3 -> days 21, 53, 77 (weeks 3, 7 or 8, 11).  -> data/processed/
  RIDV  rhesus macaques vaccinated with intradermal BCG (RIDV01-RIDV05), same corrections; used to
        compare earlier vs stronger immunity (Suppl. Section S5). CT1-CT3 at weeks 3, 7, 11 -> days
        21, 49, 77 (CT4, week 16, is beyond the simulated period).                -> data/processed/ridv/

Outputs (in data/processed/ or data/processed/ridv/):
  experimental_results.txt   BM input: 12 rows (mean, SD) for NC, NM, NMD, NMI, VFC, MAC
  per_animal.csv             one row per animal and time point
  consolidations.csv         one row per tracked lesion and time point
  lobe_frequency_week3.csv   initial-infection lobe weights
  summary.txt                human-readable table + derived constants

Usage:
  python process.py                 # corrected processing (default), unvaccinated cohort
  python process.py --cohort RIDV   # intradermal-BCG-vaccinated cohort
"""
import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RAW = HERE / "raw" / "NogueiraSuppl_File_2_RawData.xlsx"
LUNGVOL = HERE / "raw" / "lung_volumes.csv"
OUT = HERE / "processed"

ANIMALS = ["RUV01", "RUV02", "RUV03", "RUV04", "RUV06",
           "RUV07", "RUV08", "RUV09", "RUV10", "RUV11"]   # RUV05 excluded (miliary TB)
EXCLUDED = {"RUV05": "miliary TB, outside the scope of the BM (as in Nogueira et al. 2023)"}
N_CT = 3
DAYS = [21, 53, 77]
CONSOLIDATION_MA = 4.5          # mm, main axis threshold (Nogueira et al. 2023)
LOBES = ["RUL", "RML", "RLL", "LIC", "LUL", "LML", "LLL"]   # order used by the BM
LOBE_ALIASES = {"ICL": "LIC", "LII": "LLL", "LSI": "LUL"}   # LII/LSI: Spanish labels (RIDV03 sheet;
                                                               # segment numbers confirm the lobe)
COHORTS = {
    "RUV": dict(animals=ANIMALS, days=DAYS, out=OUT, excluded=EXCLUDED, label="unvaccinated rhesus"),
    "RIDV": dict(animals=[f"RIDV0{i}" for i in range(1, 6)], days=[21, 49, 77], out=OUT / "ridv",
                 excluded={}, label="intradermal-BCG-vaccinated rhesus"),
}
LABEL = COHORTS["RUV"]["label"]

# Raw-sheet column indices (0-based), consolidation sheets
C_NUMBER, C_LOBE, C_SEGMENT = 1, 2, 3
C_DAUGHTER, C_VOLUME, C_MA, C_SA = 8, 12, 16, 20      # each followed by CT1..CT4
# Isolated-micronodule sheets: counts per segment in columns 3..6 (CT1..CT4)
C_ISOL = 3


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def parse_count(v, plus_value):
    """Daughter counts may be written as '10+' in the raw file."""
    if isinstance(v, str):
        digits = re.sub(r"[^0-9.]", "", v)
        if not digits:
            return np.nan
        val = float(digits)
        return plus_value if v.strip().endswith("+") and plus_value is not None else val
    return v


def read_raw(plus_value):
    book = pd.read_excel(RAW, header=None, sheet_name=None)
    lesions, isolated = [], {}
    for a in ANIMALS:
        s = book[a].iloc[2:]
        s = s[pd.to_numeric(s[C_NUMBER], errors="coerce").notna()]
        for _, r in s.iterrows():
            lobe = LOBE_ALIASES.get(str(r[C_LOBE]).strip(), str(r[C_LOBE]).strip())
            for t in range(N_CT):
                lesions.append(dict(
                    animal=a, lesion=int(r[C_NUMBER]), lobe=lobe, segment=r[C_SEGMENT], ct=t + 1,
                    day=DAYS[t],
                    main_axis=pd.to_numeric(r[C_MA + t], errors="coerce"),
                    short_axis=pd.to_numeric(r[C_SA + t], errors="coerce"),
                    volume=pd.to_numeric(r[C_VOLUME + t], errors="coerce"),
                    daughters_raw=r[C_DAUGHTER + t],
                    daughters=pd.to_numeric(parse_count(r[C_DAUGHTER + t], plus_value), errors="coerce"),
                ))
        im = book[a + " Isol Mic"].iloc[2:, C_ISOL:C_ISOL + N_CT].apply(pd.to_numeric, errors="coerce")
        isolated[a] = im.fillna(0).sum().to_numpy()
    les = pd.DataFrame(lesions)
    for c in ["lobe"]:
        bad = set(les[c]) - set(LOBES)
        if bad:
            sys.exit(f"Unknown lobe labels in raw data: {bad}")
    return les, isolated


def read_lung_volumes():
    v = pd.read_csv(LUNGVOL).set_index("animal").lung_volume_cm3
    return {a: float(v[a]) for a in ANIMALS}


# --------------------------------------------------------------------------- #
# Corrected processing
# --------------------------------------------------------------------------- #
def process(les, isolated, lungvol, ddof=1):
    les = les.copy()
    les["present"] = les.main_axis.notna()
    les["is_consolidation"] = les.main_axis >= CONSOLIDATION_MA
    # tracked lesions that are present but still below the consolidation threshold
    les["is_small_tracked"] = les.present & ~les.is_consolidation
    rows = []
    for a in ANIMALS:
        for t in range(N_CT):
            L = les[(les.animal == a) & (les.ct == t + 1)]
            cons = L[L.is_consolidation]
            nmd = L.daughters.fillna(0).sum()
            nmi = isolated[a][t] + L.is_small_tracked.sum()
            rows.append(dict(
                animal=a, ct=t + 1, day=DAYS[t],
                NC=int(len(cons)), NMD=nmd, NMI=nmi, NM=nmd + nmi,
                VFC=cons.volume.sum() / lungvol[a],
                MAC=cons.main_axis.mean() if len(cons) else np.nan,   # undefined without consolidations
                consolidation_volume_cm3=cons.volume.sum(),
                lung_volume_cm3=lungvol[a],
                n_isolated_sheet=isolated[a][t], n_small_tracked=int(L.is_small_tracked.sum()),
            ))
    pa = pd.DataFrame(rows)
    stats = {}
    for k in ["NC", "NM", "NMD", "NMI", "VFC"]:
        g = pa.groupby("ct")[k]
        stats[k] = (g.mean().to_numpy(), g.std(ddof=ddof).to_numpy())
    # MAC: mean main axis of each animal's consolidations, then mean and SD over animals with >= 1
    # consolidation (n = 9, 10, 10), as the model computes it (one value per simulation)
    g = pa.groupby("ct")["MAC"]
    stats["MAC"] = (g.mean().to_numpy(), g.std(ddof=ddof).to_numpy())

    # lobe distribution of consolidations at week 3, averaged over animals with >=1 consolidation
    w3 = les[(les.ct == 1) & les.is_consolidation]
    per_animal = []
    for a in ANIMALS:
        x = w3[w3.animal == a]
        if len(x):
            per_animal.append([np.mean(x.lobe == lb) for lb in LOBES])
    lobe_freq = np.mean(per_animal, axis=0)

    # SA/MA ratio over all tracked lesions with both axes at the same time point
    both = les.dropna(subset=["main_axis", "short_axis"])
    sa_ma = float((both.short_axis / both.main_axis).mean())
    return pa, les, stats, lobe_freq, sa_ma


# --------------------------------------------------------------------------- #
# Output
# --------------------------------------------------------------------------- #
def write_outputs(pa, les, stats, lobe_freq, sa_ma, lungvol, args):
    OUT.mkdir(parents=True, exist_ok=True)
    with open(OUT / "experimental_results.txt", "w") as f:       # format read by exp_err.c
        for k in ["NC", "NM", "NMD", "NMI", "VFC", "MAC"]:
            m, s = stats[k]
            f.write(",".join(f"{x:.10g}" for x in m) + "\n")
            f.write(",".join(f"{x:.10g}" for x in s) + "\n")
    pa.to_csv(OUT / "per_animal.csv", index=False)
    les.drop(columns=["daughters_raw"]).to_csv(OUT / "consolidations.csv", index=False)
    pd.DataFrame({"lobe": LOBES, "frequency": lobe_freq}).to_csv(OUT / "lobe_frequency_week3.csv", index=False)
    with open(OUT / "summary.txt", "w") as f:
        f.write(summary_text(stats, lobe_freq, sa_ma, lungvol, args))


def summary_text(stats, lobe_freq, sa_ma, lungvol, args):
    s = ["Experimental targets (mean +- SD), %s, n = %d animals per time point" % (LABEL, len(ANIMALS)),
         f"Settings: '10+' daughters read as {args.plus_value}; SD ddof = {args.ddof}",
         f"Excluded: {EXCLUDED}", ""]
    s.append(f"{'':5s}" + "".join(f"{'day ' + str(d):>22s}" for d in DAYS))
    for k, (m, sd) in stats.items():
        s.append(f"{k:5s}" + "".join(f"{a:12.4g} +- {b:<7.4g}" for a, b in zip(m, sd)))
    s += ["", "Lobe frequency of consolidations at week 3: " +
          ", ".join(f"{l} {v:.4f}" for l, v in zip(LOBES, lobe_freq)),
          f"Mean SA/MA ratio: {sa_ma:.4f}", "", "Lung volumes (cm3, CT1):"]
    s += [f"  {a}: {v:.1f}" for a, v in lungvol.items()]
    return "\n".join(s) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plus-value", type=float, default=10, help="value assigned to '10+' daughter counts")
    ap.add_argument("--ddof", type=int, default=1, help="SD degrees of freedom for per-animal metrics")
    ap.add_argument("--cohort", choices=list(COHORTS), default="RUV", help="animals to process (default RUV)")
    args = ap.parse_args()
    global ANIMALS, DAYS, OUT, EXCLUDED, LABEL       # the functions below read these module globals
    c = COHORTS[args.cohort]
    ANIMALS, DAYS, OUT, EXCLUDED, LABEL = c["animals"], c["days"], c["out"], c["excluded"], c["label"]

    les, isolated = read_raw(args.plus_value)
    lungvol = read_lung_volumes()

    pa, les2, stats, lobe_freq, sa_ma = process(les, isolated, lungvol, ddof=args.ddof)
    write_outputs(pa, les2, stats, lobe_freq, sa_ma, lungvol, args)
    print()
    print(summary_text(stats, lobe_freq, sa_ma, lungvol, args))
    return 0


if __name__ == "__main__":
    sys.exit(main())
