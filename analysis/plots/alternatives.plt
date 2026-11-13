# Supplementary: best fits of the alternative models vs fit A and the data.
#   no immunity (f_imm = 1, 5 free parameters) and early immunity (t_imm = 21 d fixed, 6 free).
load "common.gp"
load D."fit_values.gp"
eval term(6, 4)
set output out("alternatives")
NMOD = 3
array MF[3]  = ["curves/fitA.dat", "curves/noimm.dat", "curves/timm21.dat"]
array MT[3]  = [sprintf("fit A (E_A = %.3f)", EOPT_fitA), sprintf("No immunity (E_A = %.2f)", EOPT_noimm), \
                sprintf("t_{imm} = 21 d (E_A = %.2f)", EOPT_timm21)]
array MC[3]  = [C_BLUE, C_RED, C_GREEN]
array MDT[3] = [1, 1, 1]
array MB[3]  = [1, 1, 1]
load "panels6.gp"
