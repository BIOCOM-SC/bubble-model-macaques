# Fig. 3 (main): in-silico experiments at the fit-A optimum: immunosuppressed (f_imm = 1) and
# segment-limited growth (-dseg), with the control (fit A) and the data.
load "common.gp"
eval term(6, 4)
set output out("experiments")
NMOD = 3
array MF[3]  = ["curves/experiments/control.dat", "curves/experiments/dseg.dat", \
                "curves/experiments/immunosuppressed.dat"]
array MT[3]  = ["Control scenario", "Segment limited", "Immunosuppressed"]
array MC[3]  = [C_BLUE, C_GREEN, C_RED]
array MDT[3] = [2, 1, 1]
array MB[3]  = [0, 1, 1]
array YMAX[6] = [30, 160, 0.5, 80, 60, 35]   # common scale; VFC fits the immunosuppressed mean (its SD bars are cut), MAC its SD bars
load "panels6.gp"
