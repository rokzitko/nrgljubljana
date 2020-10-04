def1ch[1];

(* Channel 1 and channel 2 have different hybridizations (set by parameters Gamma1 and Gamma2).
   The difference of phase of phi is taken into account by performing a rotation in the Nambu
   space by phi/2. This factor then enters the hopping matrix element between channel 1 and
   the impurity orbital. *)
   
Hc = gammaPolch1 hopphi[f[0], d[], phi/2];

H = H0 + H1 + Hc;

(* Recall: the gaps are set by bcsgap1 and bcsgap2. *)     
