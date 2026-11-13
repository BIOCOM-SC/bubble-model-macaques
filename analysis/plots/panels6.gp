# 2 x 3 panels (NC, NM, VFC / NMD, NMI, MAC): model curves (mean, ± SD error bars every SKIPEB
# days) and data.
# Set before `load "panels6.gp"` (after common.gp and the terminal/output):
#   NMOD                      number of model curves
#   array MF[NMOD]            curve files (relative to D)
#   array MT[NMOD]            legend titles
#   array MC[NMOD]            colours
#   array MDT[NMOD]           dashtypes (1 = solid)
#   array MB[NMOD]            1 = draw ± SD error bars, 0 = mean only
#   optional: SHOWDATA (0/1), array YMAX[6], EXPF (data file in D, default exp.dat)

if (!exists("SHOWDATA")) { SHOWDATA = 1 }
if (!exists("EXPF")) { EXPF = "exp.dat" }
if (!exists("YMAX")) { array YMAX[6]; do for [k=1:6] { YMAX[k] = PYMAX[k] } }
if (!exists("KEYFIRST")) { KEYFIRST = 0 }   # 1 = legend only in the first panel
if (!exists("KEYTITLE")) { KEYTITLE = "" }  # optional legend title

set title font TFONT
set key top left Left reverse
set style line 1 lw 3 lc 7 pt 0     # data error bars (black)
set multiplot layout 2,3 columns rowsfirst
set grid lc rgb "gray" lw 1
set xrange [0:80]
set xlabel XLAB
do for [k=1:6] {
    c = PCOL[k]
    set title PTITLE[k]
    set ylabel PYLAB[k] offset 1,0
    set yrange [0:YMAX[k]]
    eval panel_label(PLAB[k])
    if (KEYFIRST && k > 1) { unset key } else { set key top left Left reverse title KEYTITLE left }
    plot for [j=1:NMOD] D.MF[j] u 1:c w l lw 4 lc rgb MC[j] dt MDT[j] title MT[j], \
         for [j=1:NMOD] D.MF[j] every SKIPEB u 1:(MB[j] ? column(c) : NaN):(column(c+1)) \
             w errorbars lw 2 lt -1 pt 0 lc rgb MC[j] notitle, \
         D.EXPF u 1:(SHOWDATA ? column(c) : NaN):(column(c+1)) w p pt 5 lc -1 title DATAT, \
         D.EXPF u 1:(SHOWDATA ? column(c) : NaN):(column(c+1)) w errorbars ls 1 notitle
}
unset multiplot
