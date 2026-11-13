# Supplementary: detection-threshold sensitivity at the fit-A optimum (no refit), means of
# 10^4 simulated animals. Two figures, one threshold varied at a time (the other at its default):
#   thresholds_consol : consolidation = main axis >= 3.5 / 4.5 / 5.5 mm (micronodule > 1 mm)
#   thresholds_micro  : micronodule = main axis > 0.5 / 1.0 / 1.5 mm (consolidation >= 4.5 mm)
load "common.gp"
NMOD = 3
array MDT[3] = [1, 1, 1]
array MB[3]  = [0, 0, 0]
KEYFIRST = 1

eval term(6, 4)
set output out("thresholds_consol")
array MF[3] = ["curves/thresholds/thr_004.dat", "curves/thresholds/thr_005.dat", "curves/thresholds/thr_006.dat"]
array MT[3] = ["3.5 mm", "4.5 mm (default)", "5.5 mm"]
KEYTITLE = "Consolidation: main axis {/Symbol \263}"
array MC[3] = ["#99AAFF", C_BLUE, "#0033AA"]
load "panels6.gp"
unset output

eval term(6, 4)
set output out("thresholds_micro")
array MF[3] = ["curves/thresholds/thr_002.dat", "curves/thresholds/thr_005.dat", "curves/thresholds/thr_008.dat"]
array MT[3] = ["0.5 mm", "1 mm (default)", "1.5 mm"]
KEYTITLE = "Micronodule: main axis >"
array MC[3] = ["#F09A5C", C_ORANGE, "#BB0007"]
load "panels6.gp"
