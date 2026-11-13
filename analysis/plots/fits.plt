# Fig. 2 (main): model at the fit-A and fit-B optima vs data.
load "common.gp"
eval term(6, 4)
set output out("fits")
NMOD = 2
array MF[2]  = ["curves/fitA.dat", "curves/fitB.dat"]
array MT[2]  = ["fit A", "fit B"]
array MC[2]  = [C_BLUE, C_ORANGE]
array MDT[2] = [1, 1]
array MB[2]  = [1, 1]
load "panels6.gp"
