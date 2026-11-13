# Fig. 4 (main, replaces the t_imm line plot): "stronger vs earlier" immunity.
# NC, NM and VFC at day 77 over the t_imm x f_imm grid (t_imm 0-45 d, step 1; f_imm 0.10-1.00, step
# 0.05; 10^4 simulated animals per cell) at the fit-A values of the other parameters (run `immgrid`,
# also used for the vaccinated macaques). Star = fit-A optimum; x = cell aborted, if any (one
# simulated animal reached a model safety limit, MAX_* in main.h). Level curves: gnuplot contours
# of the cell means, linear interpolation between neighbouring grid points (no smoothing).
load "common.gp"
set terminal unknown
F = D."immgrid.dat"
# columns: 1 t_imm, 2 f_imm, 3 E_A, 4 E_B, 5 NC, 7 NM, 9 NMD, 11 NMI, 13 VFC, 15 MAC, 17 status
load D."fit_values.gp"

set contour base
unset surface
set view map
set cntrparam linear
set table D."immgrid_cont_nc.dat"
set cntrparam levels discrete 8, 12, 16
splot F u 1:2:5 w l
unset table
set table D."immgrid_cont_nm.dat"
set cntrparam levels discrete 30, 40, 50, 60
splot F u 1:2:7 w l
unset table
set table D."immgrid_cont_vfc.dat"
set cntrparam levels discrete 0.01, 0.02, 0.05, 0.1, 0.2
splot F u 1:2:13 w l
unset table
unset contour
reset
eval term(6, 2.4)
set output out("immgrid")

set palette defined (0 "#0033AA", 1 "#0044CC", 2 "#0055DD", 3 "#1166EE", 4 "#3377FF", 5 "#5588FF", \
    6 "#7799FF", 7 "#99AAFF", 8 "#B3BBFF", 9 "#C4C4C4", 10 "#D6BBA0", 11 "#E8AA7F", 12 "#F09A5C", \
    13 "#F87B3C", 14 "#FF5C1F", 15 "#FF4415", 16 "#FF2A10", 17 "#EE120C", 18 "#DD0008", 19 "#BB0007")
array GC[3] = [5, 7, 13]
array CMIN[3]; array CMAX[3]     # data ranges, before any axis range is set
do for [k=1:2] { stats F u GC[k] nooutput; CMIN[k] = STATS_min; CMAX[k] = STATS_max }
set title font TFONT
set xrange [-1:46]
set yrange [0.07:1.03]
set xtics 0,10,40 out nomirror
set ytics 0.2,0.2,1 out nomirror
set cbtics out
set xlabel "t_{imm} (days)"
set ylabel "f_{imm}" offset 1,0
set key off
unset grid
PLX = -0.3
array PT[3] = ["Consolidations (NC)", "Micronodules (NM)", "Vol. occupied by consol. (VFC)"]
array CF[3] = ["immgrid_cont_nc.dat", "immgrid_cont_nm.dat", "immgrid_cont_vfc.dat"]
array OX[3] = [0, 0.335, 0.67]
array CBT[2] = [2, 10]
NCB = 20      # colour bar drawn as NCB discrete steps (the cells keep the continuous palette)
set multiplot
do for [k=1:3] {
    set origin OX[k], 0
    set size 0.33, 1
    set lmargin at screen OX[k] + 0.065
    set rmargin at screen OX[k] + 0.255
    set tmargin at screen 0.86
    set bmargin at screen 0.2
    unset colorbox
    set border 31
    set xrange [-1:46]; set yrange [0.07:1.03]; unset logscale y
    set xtics 0,10,40 out nomirror; set ytics 0.2,0.2,1 out nomirror; unset y2tics
    set xlabel "t_{imm} (days)"; set ylabel "f_{imm}" offset 1,0
    set title PT[k] font TFONT
    eval panel_label(PLAB[k])
    if (k == 3) { set logscale cb; set cbrange [0.0003:0.3] } \
    else { unset logscale cb; set cbrange [*:*] }
    if (k < 3) { set cbrange [floor(CMIN[k] / CBT[k]) * CBT[k]:ceil(CMAX[k] / CBT[k]) * CBT[k]] }
    plot F u 1:2:(0.5):(0.025):(column(GC[k])) w boxxyerror fs solid 1.0 noborder lc palette, \
         D.CF[k] u 1:2 w l lw 2 lc rgb "black", \
         D.CF[k] u 1:2:(sprintf("%g",$3)) every 30::15 w labels font ",11" tc rgb "black", \
         F u 1:(strcol(17) eq "ok" ? 1/0 : $2) w p pt 2 ps 1.3 lc rgb "black", \
         "+" u (TFIT):(FFIT) every ::0::0 w p pt 3 ps 2.6 lw 3 lc rgb "black"
    # colour bar: NCB steps between the colour-axis limits of the panel (log steps for VFC); linear
    # limits = multiples of the tick step CBT around the data range
    CB0 = GPVAL_CB_MIN; CB1 = GPVAL_CB_MAX
    set cbrange [CB0:CB1]
    set lmargin at screen OX[k] + 0.265
    set rmargin at screen OX[k] + 0.277
    unset title; unset label 1; unset xlabel; unset ylabel; unset xtics; unset ytics
    set xrange [0:1]
    set y2range [CB0:CB1]
    set yrange [CB0:CB1]
    if (k == 3) { set logscale y; set logscale y2; set y2tics (0.001, 0.003, 0.01, 0.03, 0.1, 0.3) out nomirror } \
    else { unset logscale y; unset logscale y2; set y2tics CBT[k] out nomirror }
    lo(i) = (k == 3) ? CB0 * (CB1 / CB0) ** (i / NCB) : CB0 + (CB1 - CB0) * i / NCB
    plot sample [i=0:NCB-1:1] "+" u (0.5):(lo(i)):(0):(1):(lo(i)):(lo(i + 1)):(sqrt(lo(i) * lo(i + 1))) \
         w boxxyerror fs solid 1.0 noborder lc palette axes x1y1
    unset logscale y; unset logscale y2
}
unset multiplot
