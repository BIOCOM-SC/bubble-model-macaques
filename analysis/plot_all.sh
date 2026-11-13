#!/bin/sh
# Rebuild all analysis tables and figures.
#   ./plot_all.sh            collect.py + lung_qc.py (need runs/work/), then all figures
#   ./plot_all.sh --plots    figures only, from the committed analysis/data/
#   ./plot_all.sh --plots fits experiments   only these scripts
# gnuplot writes EPS (postscript terminal); PDF via epstopdf and PNG
# previews via ghostscript.
set -e
cd "$(dirname "$0")"
PY=${PYTHON:-python3}   # e.g. PYTHON=/abs/path/.venv/bin/python (the script cd's into analysis/)
if [ "$1" != "--plots" ]; then
    $PY collect.py
    $PY lung_qc.py
else
    shift
fi
SCRIPTS=${*:-"fits experiments immgrid ridv_vs_ruv ridv_heatmap ridv_fits families alternatives immsplit thresholds sweeps_err prcc minerr tradeoff per_animal lung_qc"}
mkdir -p out
cd plots
for f in $SCRIPTS; do gnuplot $f.plt; done
cd ../out
for e in *.eps; do
    b=${e%.eps}
    if [ ! -f "$b.png" ] || [ "$e" -nt "$b.png" ]; then
        epstopdf "$e" --outfile="$b.pdf"
        gs -q -dSAFER -dBATCH -dNOPAUSE -dEPSCrop -sDEVICE=png16m -r300 -dTextAlphaBits=4 -dGraphicsAlphaBits=4 -sOutputFile="$b.png" "$e"
    fi
done
echo "figures in analysis/out/ (eps, pdf, png)"
