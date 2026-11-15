# Supplementary: parameter trade-offs of fit A.
# Single evaluations with E_A < THR x E_A* (global + local searches), projected on N0-v0 (A),
# N0-R_max0 (B) and v0-R_max0 (C), coloured by E_A (best on top). Open circles: refined candidates
# whose mean E_A over 10 replicates is within 2 sigma_E of the optimum (statistically equivalent).
# Star: refined fit-A optimum.
load "common.gp"
eval term(6, 2.6)
set output out("tradeoff")
F = D."fitA_evals.dat"           # v0 fvimm r0 b nata tc n0 E_A
G = D."fitA_refined_near.dat"    # v0 fvimm r0 b nata tc n0 mean_E_A
load D."fit_values.gp"
THR = 1.2
set palette defined (0 "#0033AA", 1 "#3377FF", 2 "#99AAFF", 3 "#C4C4C4", 4 "#F09A5C", 5 "#FF4415", 6 "#BB0007")
set title font TFONT
set multiplot layout 1,3 margins 0.08,0.88,0.18,0.86 spacing 0.09,0
set grid lc rgb "gray" lw 1
PLX = -0.2
set cbrange [E0*0.95:THR*E0]
set cblabel "E_A" offset 0,0
set cbtics 0.02
unset colorbox
STAR = "w p pt 3 ps 2.6 lw 3 lc rgb 'black' notitle"
STARW = "w p pt 3 ps 2.6 lw 8 lc rgb 'white' notitle"
NEAR = "w p pt 6 ps 1.2 lw 2 lc rgb 'black' notitle"

set xlabel "N_0"
set xrange [9.2:10.6]
set xtics 9.2,0.4,10.6
set title "N_0 vs v_0"
eval panel_label("A")
set yrange [1.2:1.4]
set ytics 0.05
set ylabel "v_0 (day^{-1})" offset 1,0
plot F u ($8 < THR*E0 ? $7 : NaN):1:8 w p pt 7 ps 0.6 lc palette notitle, \
     G u 7:1 @NEAR, \
     "+" u (N0F):(V0F) every ::0::0 @STARW, \
     "+" u (N0F):(V0F) every ::0::0 @STAR

set title "N_0 vs R_{max 0}"
eval panel_label("B")
set yrange [5.2:6.0]
set ytics 0.2
set ylabel "R_{max 0} (mm)" offset 1,0
plot F u ($8 < THR*E0 ? $7 : NaN):3:8 w p pt 7 ps 0.6 lc palette notitle, \
     G u 7:3 @NEAR, \
     "+" u (N0F):(R0F) every ::0::0 @STARW, \
     "+" u (N0F):(R0F) every ::0::0 @STAR

set title "v_0 vs R_{max 0}"
eval panel_label("C")
set colorbox user origin 0.90, 0.18 size 0.015, 0.68
set xrange [1.2:1.4]
set xtics 0.05
set xlabel "v_0 (day^{-1})"
set yrange [5.2:6.0]
set ytics 0.2
set ylabel "R_{max 0} (mm)" offset 1,0
plot F u ($8 < THR*E0 ? $1 : NaN):3:8 w p pt 7 ps 0.6 lc palette notitle, \
     G u 1:3 @NEAR, \
     "+" u (V0F):(R0F) every ::0::0 @STARW, \
     "+" u (V0F):(R0F) every ::0::0 @STAR
unset multiplot
