# Supplementary: progress of the searches.
# Best-so-far objective of the global DIRECT search (~5000 evaluations) followed by the 5 local
# DIRECT searches (~500 evaluations each, drawn after the global search); dashed: refined optimum
# (lowest mean over 10 replicates). Grey dots: single evaluations (10^4 simulated animals each).
load "common.gp"
eval term(6, 5.4)
set output out("evolution_min_error")
array FIT[4] = ["fitA", "fitB", "noimm", "timm21"]
array FT[4]  = ["Fit A (E_A)", "Fit B (E_B)", "No immunity (E_A)", "t_{imm} = 21 d (E_A)"]
array YR0[4] = [0.2, 1.0, 5.2, 3.3]
array YR1[4] = [0.7, 2.5, 7.0, 5.0]
PLX = -0.12
array LC[5]  = ["#0072B2", "#D55E00", "#E69F00", "#009E73", "#CC79A7"]
set title font TFONT
set grid lw 1
set multiplot
# legend strip at the top, so no curve overlaps it
set origin 0, 0.92
set size 1, 0.08
set xrange [0:1]
set yrange [0:1]
unset xtics; unset ytics
set border 0
set key at screen 0.5, 0.96 center horizontal maxrows 1 box font ",11" samplen 2
plot NaN w l lw 4 lc rgb "black" title "Global", \
     for [j=1:5] NaN w l lw 4 lc rgb LC[j] title "Local ".j, \
     NaN w l dt 2 lw 2 lc rgb "#444444" title "Refined"
unset key
set border
set xtics; set ytics
array ORX[4] = [0, 0.5, 0, 0.5]
array ORY[4] = [0.46, 0.46, 0, 0]
set xrange [0:5700]
set xtics 1000
set xlabel "Function evaluations"
set ylabel "Minimum error (unitless)"
do for [i=1:4] {
    f = D."minerr/".FIT[i]
    R = real(system("awk 'NR==2{print $1}' ".f."_refined.dat"))
    set origin ORX[i], ORY[i]
    set size 0.5, 0.46
    set title FT[i]
    set yrange [YR0[i]:YR1[i]]
    eval panel_label(PLAB[i])
    plot f."_global.dat" u 1:2 w p pt 7 ps 0.15 lc rgb "#BBBBBB" notitle, \
         for [j=1:5] f."_local".j.".dat" u 1:2 w p pt 7 ps 0.15 lc rgb "#BBBBBB" notitle, \
         f."_global.dat" u 1:3 w l lw 4 lc rgb "black" notitle, \
         for [j=1:5] f."_local".j.".dat" u 1:3 w l lw 4 lc rgb LC[j] notitle, \
         R w l dt 2 lw 2 lc rgb "#444444" notitle
}
unset multiplot
