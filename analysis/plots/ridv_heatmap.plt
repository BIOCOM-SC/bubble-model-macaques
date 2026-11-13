# Suppl. Fig. 12: E_A over the immgrid grid (t_imm 0-45 d x f_imm 0.10-1.00, other parameters at fit
# A; the grid of Fig. 4) against the unvaccinated (A) and the intradermal-BCG-vaccinated (B) targets
# (data/vaccinated_err.dat, collect.py vaccinated). Star = fit A; open circle = cell with the lowest
# error. Log colour scale, per panel; colour bar in NCB steps (the cells keep the continuous palette).
load "common.gp"
load D."fit_values.gp"
load D."vaccinated.gp"
DX = 0.5; DY = 0.025
eval term(4.4, 2.4)
set output out("ridv_heatmap")
array EC[2] = [3, 4]     # columns of vaccinated_err.dat
array PT[2] = ["Unvaccinated (n=10)", "Intradermal BCG (n=5)"]
array OX[2] = [0, 0.5]
array EM[2] = [EMIN_RUV, EMIN_RIDV]; array EX[2] = [EMAX_RUV, EMAX_RIDV]
NCB = 20      # colour bar drawn as NCB discrete (log) steps; the cells keep the continuous palette
set palette defined (0 "#0033AA", 1 "#3377FF", 2 "#99AAFF", 3 "#C4C4C4", 4 "#F09A5C", 5 "#FF4415", 6 "#BB0007")
set title font TFONT
set xrange [-1:46]
set yrange [0.07:1.03]
set xtics 0,10,40 out nomirror
set ytics 0.2,0.2,1 out nomirror
set logscale cb
set xlabel "t_{imm} (days)"
set ylabel "f_{imm}" offset 1,0
set key off
unset grid
PLX = -0.3
set multiplot
do for [k=1:2] {
    F = D."vaccinated_err.dat"
    EMIN = EM[k]
    set origin OX[k], 0
    set size 0.5, 1
    set lmargin at screen OX[k] + 0.09
    set rmargin at screen OX[k] + 0.38
    set tmargin at screen 0.86
    set bmargin at screen 0.2
    unset colorbox
    set xrange [-1:46]; set yrange [0.07:1.03]; unset logscale y
    set xtics 0,10,40 out nomirror; set ytics 0.2,0.2,1 out nomirror; unset y2tics
    set xlabel "t_{imm} (days)"; set ylabel "f_{imm}" offset 1,0
    set title PT[k] font TFONT
    eval panel_label(PLAB[k])
    set cbrange [EMIN:EX[k]]
    plot F u 1:2:(DX):(DY):(column(EC[k])) w boxxyerror fs solid 1.0 noborder lc palette, \
         F u 1:(column(EC[k]) == EMIN ? $2 : 1/0) w p pt 6 ps 2 lw 3 lc rgb "black", \
         "+" u (TFIT):(FFIT) every ::0::0 w p pt 3 ps 2.6 lw 3 lc rgb "black"
    # colour bar: NCB log steps between the colour-axis limits of the panel
    CB0 = EMIN; CB1 = EX[k]
    set lmargin at screen OX[k] + 0.395
    set rmargin at screen OX[k] + 0.413
    unset title; unset label 1; unset xlabel; unset ylabel; unset xtics; unset ytics
    set xrange [0:1]
    set logscale y; set logscale y2
    set yrange [CB0:CB1]; set y2range [CB0:CB1]
    set y2tics out nomirror (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
    lo(i) = CB0 * (CB1 / CB0) ** (i / NCB)
    plot sample [i=0:NCB-1:1] "+" u (0.5):(lo(i)):(0):(1):(lo(i)):(lo(i + 1)):(sqrt(lo(i) * lo(i + 1))) \
         w boxxyerror fs solid 1.0 noborder lc palette
    unset logscale y; unset logscale y2
}
unset multiplot
