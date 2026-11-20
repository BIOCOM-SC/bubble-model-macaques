# Supplementary: objective E_A along each single-parameter sweep (+-25 % around fit A, 20 values,
# one evaluation of 10^4 simulated animals each) compared with the noise of the objective:
# grey band = E_A* +- 2 sigma_E (refined fit-A optimum and the SD of its 10 replicates).
load "common.gp"
eval term(6, 3.4)
set output out("sweeps_err")
F = D."sweeps_err.dat"
E0 = real(system("awk 'NR==1{print $3}' ".F))
SIG = real(system("awk 'NR==1{print $5}' ".F))
array PN[8] = ["v_0", "f_{imm}", "R_{max 0}", "{/Symbol b}", "{/Symbol r}", "r_{min}", "N_0", "t_{imm}"]
array LB[8] = ["A", "B", "C", "D", "E", "F", "G", "H"]
set title font TFONT
set multiplot layout 2,4
set grid lc rgb "gray" lw 1
set logscale y
set yrange [0.2:10]
set ytics (0.2, 0.5, 1, 2, 5, 10)
set xrange [0.73:1.27]
set xtics 0.75,0.25,1.25
set style fill solid 0.5 noborder
set xlabel "value / fit A"
PLX = -0.3
do for [i=1:8] {
    set title PN[i]
    eval panel_label(LB[i])
    if (i == 1 || i == 5) { set ylabel "E_A" offset 1,0 } else { unset ylabel }
    plot "+" u 1:(E0-2*SIG):(E0+2*SIG) w filledcurves lc rgb "#999999" notitle, \
         F index (i-1) u 1:8 w lp pt 7 ps 0.6 lw 3 lc rgb C_BLUE notitle
}
unset multiplot
