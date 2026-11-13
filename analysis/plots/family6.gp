# 2 x 3 panels (NC, NM, VFC / NMD, NMI, MAC) for a one-parameter family of 20 model runs, in
# 20 fixed colours (deep blue -> grey -> deep red, by run index), colour bar on the right. Added: fit A (black, dashed) and the data.
# Set before `load "family6.gp"`: load D."fam_<tag>.gp" (NFAM, FF[], FV[]), CBLAB, optional
# CBTICS (step), CBFMT, YMAX[6], SHOWDATA, REF (curve file, "" = none), REFT.

if (!exists("SHOWDATA")) { SHOWDATA = 1 }
if (!exists("YMAX")) { array YMAX[6]; do for [k=1:6] { YMAX[k] = PYMAX[k] } }
if (!exists("REF")) { REF = "curves/fitA.dat" }
if (!exists("REFT")) { REFT = "fit A" }
if (!exists("CBFMT")) { CBFMT = "%g" }

array FC[20] = ["#0033AA", "#0044CC", "#0055DD", "#1166EE", "#3377FF", "#5588FF", "#7799FF", \
                "#99AAFF", "#B3BBFF", "#C4C4C4", "#D6BBA0", "#E8AA7F", "#F09A5C", "#F87B3C", \
                "#FF5C1F", "#FF4415", "#FF2A10", "#EE120C", "#DD0008", "#BB0007"]
set palette defined (0 "#0033AA", 1 "#0044CC", 2 "#0055DD", 3 "#1166EE", 4 "#3377FF", 5 "#5588FF", \
    6 "#7799FF", 7 "#99AAFF", 8 "#B3BBFF", 9 "#C4C4C4", 10 "#D6BBA0", 11 "#E8AA7F", 12 "#F09A5C", \
    13 "#F87B3C", 14 "#FF5C1F", 15 "#FF4415", 16 "#FF2A10", 17 "#EE120C", 18 "#DD0008", 19 "#BB0007")
set style line 1 lw 3 lc 7 pt 0     # data error bars (black)

set multiplot
set grid lc rgb "gray" lw 1

# colour bar (own panel on the right)
set cbrange [FV[1]:FV[NFAM]]
if (exists("CBTICS")) { set cbtics CBTICS }
set format cb CBFMT
set colorbox vertical user origin 0.905, 0.05 size 0.02, 0.9
set cblabel CBLAB font TFONT offset 0,0
set border 0
set origin 0.9, 0
set size 0.1, 1
set xrange [0:1]
set yrange [0:1]
unset xtics; unset ytics
unset key
plot for [j=1:20] NaN w l lc palette frac (j-1)/20. notitle
unset colorbox

set title font TFONT offset 0,0
set border
set xtics; set ytics
set grid lw 1
set key top left Left reverse
array ORX[6] = [0, 0.30, 0.60, 0, 0.30, 0.60]
array ORY[6] = [0.5, 0.5, 0.5, 0, 0, 0]
do for [k=1:6] {
    c = PCOL[k]
    set origin ORX[k], ORY[k]
    set size 0.30, 0.5
    set title PTITLE[k]
    set xlabel XLAB
    set ylabel PYLAB[k] offset 1,0
    set xrange [0:80]
    set yrange [0:YMAX[k]]
    eval panel_label(PLAB[k])
    if (k == 1) { set key } else { unset key }
    plot for [j=1:NFAM] D.FF[j] u 1:c w l lw 3 lc rgb FC[j] notitle, \
         D.(REF ne "" ? REF : "exp.dat") u 1:(REF ne "" ? column(c) : NaN) w l lw 3 dt 2 lc rgb "black" title REFT, \
         D."exp.dat" u 1:(SHOWDATA ? column(c) : NaN) w p pt 5 lc -1 title DATAT, \
         D."exp.dat" u 1:(SHOWDATA ? column(c) : NaN):(column(c+1)) w errorbars ls 1 notitle
}
unset multiplot
