(* Free impurity with Tz=0; the +/-1 channels are spectators. *)
def3ch[1];
MyAssert[SYMTYPE === "QSTZ" && Ninit == 0];
qstzImpurityCharges = {d[] -> 0};
snegrealconstants[eps, V];
H0 = H0 /. f[ca_, ch_, sp_] :> f[ca, 1 - ch, sp];
Himp = H1 = eps number[d[]];
Hhyb = Hc = V hop[f[0], d[]];
H = Expand[H0 + Himp + Hc];
