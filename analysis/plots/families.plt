# Supplementary: one-parameter families at the fit-A optimum (10^4 simulated animals per curve).
#   exp_fimm_all : f_imm sweep fit-A value .. 1 (no immunity)
#   sweep_<p>    : single-parameter sweeps, +-25 % around fit A, 20 values
load "common.gp"

do for [tag in "fimm sweep_v0 sweep_fvimm sweep_r0 sweep_b sweep_nata sweep_rmin sweep_n0 sweep_tc"] {
    reset
    load D."fam_".tag.".gp"
    undefine YMAX CBTICS CBFMT
    CBFMT = "%.2f"
    if (tag eq "fimm")       { name = "exp_fimm_all"; CBLAB = "f_{imm}"; \
        array YMAX[6] = [30, 160, 0.4, 100, 100, 30] }
    if (tag eq "sweep_v0")    { name = tag; CBLAB = "v_0 (day^{-1})" }
    if (tag eq "sweep_fvimm") { name = "sweep_fimm"; CBLAB = "f_{imm}"; CBFMT = "%.3f" }
    if (tag eq "sweep_r0")    { name = "sweep_rmax"; CBLAB = "R_{max 0} (mm)" }
    if (tag eq "sweep_b")     { name = tag; CBLAB = "{/Symbol b} (mm^{-1})"; CBFMT = "%.3f" }
    if (tag eq "sweep_nata")  { name = "sweep_rho"; CBLAB = "{/Symbol r} (day^{-1} mm^{-1})"; CBFMT = "%.3f" }
    if (tag eq "sweep_rmin")  { name = tag; CBLAB = "r_{min} (mm)"; CBFMT = "%.3f" }
    if (tag eq "sweep_n0")    { name = tag; CBLAB = "N_0"; CBFMT = "%.1f" }
    if (tag eq "sweep_tc")    { name = "sweep_timm"; CBLAB = "t_{imm} (days)"; CBFMT = "%.1f" }
    if (strstrt(tag, "sweep_") == 1) { \
        CBTICS = (FV[NFAM] - FV[1]) / 11 }
    eval term(6, 4)
    set output out(name)
    load "family6.gp"
    unset output
}
