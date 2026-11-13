# Supplementary: PRCC (1000 LHS points, uniform +-25 % around fit A, 10^4 simulated animals per
# point). prcc_final: day 77, bars +- 1 SE; prcc: time course. Grey lines: |PRCC| below which the coefficient is not significant (alpha = 0.05).
# Time course: days < 10 shaded (outputs zero at many LHS points before detection, PRCC undefined).
load "common.gp"
RC = real(system("awk '{print $1}' ".D."prcc_rcrit.dat"))
array POUT[6] = [0, 1, 4, 2, 3, 5]    # block index of NC, NM, VFC, NMD, NMI, MAC (panel order)
array PTP[6]  = ["Consolidations (NC)", "Micronodules (NM)", "Volume occupied by consol. (VFC)", \
                 "Daughter micronodules (NMD)", "Isolated micronodules (NMI)", "Main axis of consolidations (MAC)"]

# ---------------------------------------------------------------- day 77
eval term(6, 4)
set output out("prcc_final")
prcc = D."prcc_final.dat"
set style line 1 lw 3 lc rgb "#999999"  # Grey
set style line 2 lw 2 lc rgb "#0072bd"  # Blue
set xtics ("v_0" 0, "f_{imm}" 1, "R_{max 0}" 2, "{/Symbol b}" 3, \
           "{/Symbol r}" 4, "r_{min}" 5, "N_0" 6, "t_{imm}" 7)
set xtics rotate by 45 right offset 0.5,0
set title font TFONT
set grid lw 1
set style histogram errorbars gap 2 lw 1
set style data histograms
set style fill solid 1.0 border -1
set boxwidth 2.2
set yrange [-1:1]
set ylabel "PRCC" offset 1,0
set multiplot layout 2, 3
do for [k=1:6] {
    set title PTP[k]
    eval panel_label(PLAB[k])
    set arrow 1 from graph 0, first RC to graph 1, first RC nohead ls 1 front
    set arrow 2 from graph 0, first -RC to graph 1, first -RC nohead ls 1 front
    plot prcc i POUT[k] using 3:4 notitle ls 2
}
unset multiplot
unset output

# ---------------------------------------------------------------- time course
reset
eval term(6, 4)
set output out("prcc")
prcc = D."prcc_daily.dat"     # one block per output: day, then the 8 parameters
set style line 1 lw 3 lc rgb "#999999"  # Grey
set style line 2 lw 2 lc rgb "#8B4513"  # Saddle Brown
set style line 3 lw 2 lc rgb "#DC143C"  # Red
set style line 4 lw 2 lc rgb "#E69F00"  # Orange
set style line 5 lw 2 lc rgb "#66C266"  # Green
set style line 6 lw 2 lc rgb "#CC79A7"  # Purple
set style line 7 lw 2 lc rgb "#56B4E9"  # Light Blue
set style line 8 lw 2 lc rgb "#000080"  # Navy
set style line 9 lw 2 lc rgb "#007F7F"  # Teal

set multiplot
# title and legend strip
set origin 0, 0.9
set size 1, 0.1
set xrange [0:1]
set yrange [0:1]
unset xtics
unset ytics
set border 0
set key at screen 0.5, 0.95 center horizontal
set key box spacing 1.75 font ",13"
plot NaN ls 2 lw 4 title "v_0", NaN ls 3 lw 4 title "f_{imm}", NaN ls 4 lw 4 title "R_{max 0}", \
     NaN ls 5 lw 4 title "{/Symbol b}", NaN ls 6 lw 4 title "{/Symbol r}", NaN ls 7 lw 4 title "r_{min}", \
     NaN ls 8 lw 4 title "N_0", NaN ls 9 lw 4 title "t_{imm}"
unset title
set title font TFONT
unset key
set grid lw 1
set border
set xtics
set ytics
array ORX[6] = [0, 0.33, 0.66, 0, 0.33, 0.66]
array ORY[6] = [0.45, 0.45, 0.45, 0, 0, 0]
do for [k=1:6] {
    set origin ORX[k], ORY[k]
    set size 0.33, 0.45
    set title PTP[k]
    set xlabel XLAB
    set ylabel "PRCC" offset 1,0
    set xrange [1:80]
    set yrange [-1:1]
    eval panel_label(PLAB[k])
    set arrow 1 from 1, RC to 80, RC nohead ls 1 front
    set arrow 2 from 1, -RC to 80, -RC nohead ls 1 front
    # days < 10: outputs still zero (or nearly) at many LHS points (lesions seeded at t_min = 7,
    # micronodules detectable from day 8, consolidations from day 9-10) -> PRCC undefined/unreliable
    set object 1 rect from 1, -1 to 10, 1 fc rgb "#E0E0E0" fs solid 1.0 noborder behind
    plot for [j=1:8] prcc i POUT[k] using 1:(column(j+1)) w l ls (j+1) notitle
}
unset multiplot
