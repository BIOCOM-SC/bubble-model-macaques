# Suppl. Fig. 13: best fits to the intradermal-BCG-vaccinated macaques (data/exp_ridv.dat) when only
# f_imm, only t_imm, or both may change, the other parameters at fit A: lowest-E_A cells of the immgrid
# run (data/vaccinated.gp, collect.py vaccinated). Same layout and style as Fig. 2 (panels6.gp).
load "common.gp"
load D."vaccinated.gp"
eval term(6, 4)
set output out("ridv_fits")
NMOD = 3
array MF[3]  = [RIDV_F_FILE, RIDV_T_FILE, RIDV_B_FILE]
array MT[3]  = [sprintf("only f_{imm} free (%g, %.2f)", RIDV_F_T, RIDV_F_F), \
                sprintf("only t_{imm} free (%g, %.2f)", RIDV_T_T, RIDV_T_F), \
                sprintf("both free (%g, %.2f)", RIDV_B_T, RIDV_B_F)]
array MC[3]  = [C_GREEN, C_ORANGE, C_BLUE]
array MDT[3] = [1, 1, 1]
array MB[3]  = [0, 0, 1]
KEYFIRST = 1
KEYTITLE = "Best fit with (t_{imm}, f_{imm}):"
DATAT = "Intradermal BCG"
EXPF = "exp_ridv.dat"
load "panels6.gp"
