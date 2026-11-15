# Shared settings for all figures (postscript eps, Helvetica 12, bold titles 13).
# Run gnuplot from analysis/plots/ (see ../plot_all.sh, which also converts EPS -> PDF/PNG):
#   gnuplot fits.plt
# Data: ../data/   Output: ../out/<name>.eps

D = "../data/"
O = "../out/"

# terminal: call as  eval term(width_in, height_in)  then  set output out("name")
term(w, h) = sprintf("set terminal postscript eps color font 'Helvetica,12' enhanced size %g,%g", w, h)
out(name) = O.name.".eps"
TFONT = "Helvetica Bold,13"

# colours
C_BLUE = "#0072bd"; C_ORANGE = "#d95319"; C_GREEN = "#77ac30"; C_RED = "#a2142f"
C_YELLOW = "#edb120"; C_PURPLE = "#7e2f8e"; C_CYAN = "#4dbeee"; C_GREY = "#999999"

# The six outputs. Columns in every curve file (and data/exp.dat): day, then mean/SD pairs
# NC 2/3, NM 4/5, NMD 6/7, NMI 8/9, VFC 10/11, MAC 12/13.
# Panel order (2 x 3): NC, NM, VFC / NMD, NMI, MAC.
array PCOL[6]   = [2, 4, 10, 6, 8, 12]
array PTITLE[6] = ["Consolidations (NC)", "Micronodules (NM)", "Vol. occupied by consol. (VFC)", \
                   "Daughter micronodules (NMD)", "Isolated micronodules (NMI)", "Main axis of consolidations (MAC)"]
array PYLAB[6]  = ["Number", "Number", "Lung volume fraction", "Number", "Number", "Main axis (mm)"]
array PYMAX[6]  = [30, 160, 0.15, 80, 60, 30]     # common scale of all six-panel figures
array PLAB[6]   = ["A", "B", "C", "D", "E", "F"]

XLAB = "time (days)"
DATAT = "Unvaccinated (n=10)"
SKIPEB = 12        # model error bars (± SD) every 12 days

# panel label (reviewer request): call  eval panel_label("A")
PLX = -0.2
panel_label(s) = sprintf("set label 1 '%s' at graph %g, graph 1.09 font 'Helvetica Bold,14'", s, PLX)
