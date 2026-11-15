# Suppl. Fig. 11: experimental mean +- SD, unvaccinated (n = 10, data/exp.dat) vs intradermal-BCG-
# vaccinated (n = 5, data/exp_ridv.dat) macaques. Points at the scan days (21, 53 or 49, 77); at days
# 21 and 77 the two groups are drawn one day apart so that their error bars do not overlap.
load "common.gp"
eval term(6, 4)
set output out("ridv_vs_ruv")
array YMAX[6] = [30, 160, 0.15, 80, 60, 30]
set title font TFONT
set multiplot layout 2,3 columns rowsfirst
set grid lc rgb "gray" lw 1
set xrange [0:80]
set xtics 21, 14, 77     # scans: days 21, 53 (unvaccinated) or 49 (intradermal), 77
set xlabel XLAB
do for [k=1:6] {
    c = PCOL[k]
    set title PTITLE[k]
    set ylabel PYLAB[k] offset 1,0
    set yrange [0:YMAX[k]]
    eval panel_label(PLAB[k])
    if (k == 1) { set key top left Left reverse } else { unset key }
    plot D."exp.dat" u ($1 == 53 ? $1 : $1-1):c:(column(c+1)) w yerrorlines lw 2 pt 5 lc rgb C_RED title "Unvaccinated (n=10)", \
         D."exp_ridv.dat" u ($1 == 49 ? $1 : $1+1):c:(column(c+1)) w yerrorlines lw 2 pt 7 lc rgb C_PURPLE title "Intradermal BCG (n=5)"
}
unset multiplot
