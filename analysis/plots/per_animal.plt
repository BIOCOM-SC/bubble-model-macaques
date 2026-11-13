# Supplementary: paired data, one line per animal (CT at days 21, 53, 77; n = 10
# animals, RUV05 excluded: miliary TB), with the mean ± SD used as fit targets (black).
# MAC: mean main axis of each animal's consolidations (missing if the animal had none).
# Legend in a strip at the top, so no curve overlaps it.
load "common.gp"
eval term(6, 4.4)
set output out("per_animal")
NAMES = system("cat ".D."per_animal_names.txt")
array AC[10] = ["#0072bd", "#d95319", "#edb120", "#7e2f8e", "#77ac30", "#4dbeee", "#a2142f", \
                "#999999", "#000080", "#007F7F"]
array YM[6] = [30, 160, 0.15, 100, 100, 30]
set style line 1 lw 3 lc 7 pt 0

set multiplot
# legend strip
set origin 0, 0.9
set size 1, 0.1
set xrange [0:1]
set yrange [0:1]
unset xtics; unset ytics
set border 0
set key at screen 0.5, 0.95 center horizontal maxrows 2 box font ",11" samplen 2 spacing 1.1
plot for [i=1:10] NaN w lp pt 7 ps 0.6 lw 2 lc rgb AC[i] title word(NAMES, i), \
     NaN w p pt 5 lc -1 title "Mean {/Symbol \261} SD"
unset key

set border
set xtics; set ytics
set title font TFONT
set grid lc rgb "gray" lw 1
array ORX[6] = [0, 0.333, 0.666, 0, 0.333, 0.666]
array ORY[6] = [0.45, 0.45, 0.45, 0, 0, 0]
do for [k=1:6] {
    c = PCOL[k]
    set origin ORX[k], ORY[k]
    set size 0.333, 0.45
    set title PTITLE[k]
    set xrange [0:80]
    set xtics (21, 53, 77)
    set xlabel XLAB
    set ylabel PYLAB[k] offset 1,0
    set yrange [0:YM[k]]
    eval panel_label(PLAB[k])
    plot for [i=1:10] D."per_animal.dat" index (i-1) u 1:c w lp pt 7 ps 0.6 lw 2 lc rgb AC[i] notitle, \
         D."exp.dat" u ($1+1.5):c w p pt 5 lc -1 notitle, \
         D."exp.dat" u ($1+1.5):c:(column(c+1)) w errorbars ls 1 notitle
}
unset multiplot
