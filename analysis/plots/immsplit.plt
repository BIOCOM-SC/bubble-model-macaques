# Supplementary: which effect of immunity matters? At the fit-A optimum, immunity acting on
# both lesion growth speed and maximum radius (control), on growth speed only (r_max = R_max0),
# or on maximum radius only (v = v_0).
load "common.gp"
eval term(6, 4)
set output out("immsplit")
NMOD = 3
array MF[3]  = ["curves/immsplit/both.dat", "curves/immsplit/speed_only.dat", "curves/immsplit/radius_only.dat"]
array MT[3]  = ["Speed and radius (fit A)", "Growth speed only", "Maximum radius only"]
array MC[3]  = [C_BLUE, C_YELLOW, C_PURPLE]
array MDT[3] = [1, 1, 1]
array MB[3]  = [1, 1, 1]
array YMAX[6] = [30, 160, 0.4, 140, 100, 30]   # common scale, VFC up to 0.4, NMD up to 140
load "panels6.gp"
