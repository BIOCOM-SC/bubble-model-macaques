# Supplementary: virtual lung QC. A: distance of each terminal to the nearest segment wall
# and lobe wall. B: terminal diameters. C: per lobe, share of terminals that
# lie outside their own lobe mesh: outside all lobe meshes, or inside another lobe. D: distance of those
# terminals to their own lobe mesh.
load "common.gp"
eval term(6, 4.8)
set output out("lung_qc")
F = D."lung_terminals.dat"   # x y z seg dseg lobe dlobe diam p_init R_i in_own in_own_vote in_any dist_mesh acceptance
set title font TFONT
set multiplot layout 2,2 margins 0.09,0.98,0.09,0.92 spacing 0.12,0.16
set grid lw 1
set style fill solid border -1
PLX = -0.17

binwidth = 1.0
bin1(x) = binwidth * floor(x / binwidth) + binwidth / 5
bin2(x) = binwidth * floor(x / binwidth) + binwidth * 3.0 / 5
set title "Distance to closest wall"
eval panel_label("A")
set xlabel "Distance to closest wall (mm)"
set ylabel "Count"
set xrange [0:15]
set xtics 2
set boxwidth binwidth * 2.0/5
set key top right font ",10"
plot F using (bin1(column(5))):(1) smooth frequency with boxes title 'Segment walls' lc rgb "#0072bd", \
     F using (bin2(column(7))):(1) smooth frequency with boxes title 'Lobe walls' lc rgb "#d95319"

bw = 0.025
bin(x) = bw*floor(x/bw) + bw/2
set title "Terminal size"
eval panel_label("B")
set xlabel "Terminal diameter d_i (mm)"
set xrange [0:0.6]
set xtics 0.1 nomirror
set boxwidth bw*0.9
unset key
plot F u (bin($8)):(1) smooth freq w boxes lc rgb "#77ac30" notitle
set xtics mirror

set title "Terminals outside own lobe"
eval panel_label("C")
set xlabel ""
set ylabel "Terminals of the lobe (%)" offset 1,0
set xrange [-0.6:6.6]
set yrange [0:25]
set xtics ("RUL" 0, "RML" 1, "RLL" 2, "LIC" 3, "LUL" 4, "LML" 5, "LLL" 6)
set boxwidth 0.6
set key top left font ",10"
# stacked: outside all lobes (bottom) + inside another lobe (top) = outside own lobe
plot D."lung_lobes.dat" u 1:(100*($7+$8)) w boxes lc rgb "#77ac30" title "Inside another lobe", \
     D."lung_lobes.dat" u 1:(100*$7) w boxes lc rgb "#0072bd" title "Outside all lobes"

bwd = 0.5
bind(x) = bwd*floor(x/bwd) + bwd/2
set title "Outside terminals: distance to own lobe"
eval panel_label("D")
set xlabel "Distance to own lobe surface (mm)"
set ylabel "Count" offset 0,0
set xrange [0:14]
set yrange [0:*]
set xtics 2
set boxwidth bwd*0.9
unset key
plot "< awk '!/^#/ && $11 == 0' ".F u (bind($14)):(1) smooth freq w boxes lc rgb "#d95319" notitle
unset multiplot
