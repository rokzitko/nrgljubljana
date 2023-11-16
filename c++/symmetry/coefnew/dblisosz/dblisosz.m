(* 
SYMTYPE=DBLISOSZ  SU(2)_{iso} x SU(2)_{iso} x U(1)_{spin}
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2009-2023
*)

PREFIX = "dblisosz-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
basisops1 = {First[basisops]};
basisops2 = {Last[basisops]};
basis = quickDBLSZ[basisops1, basisops2, quickISOSZbasis];

snegrealconstants[n, NN];

(* IMPORTANT: the following choice of nnop corresponds
to a particular partition of the lattice sites into two
sets! *)

nnop[a[]] = n;
nnop[b[]] = n;

Tminus1 = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops1];
Tminus2 = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops2];

mult[{i1_, i2_, _}] := (2i1+1)(2i2+1);

basisprep[];

(* Simplify n: simply putting n \[Element] Integers interferes with 
snegrealconstants[n] *)

simpln[expr_] := 
  Simplify[(expr /. n->nnn), nnn \[Element] Integers] /. nnn->n;

simplNN[expr_] := 
  Simplify[(expr /. NN->nnn), nnn \[Element] Integers] /. nnn->NN;

(* Reduce isospin by one unit *)

isospin1down[op_] := Module[{tmp},
 tmp = nc[Tminus1, op];
 tmp = zeroonvac[tmp];
 If[tmp =!= 0,
  tmp = tmp / Sqrt[simpln[scalarproductop[tmp,tmp]]],
  0
 ]
];

isospin2down[op_] := Module[{tmp},
 tmp = nc[Tminus2, op];
 tmp = zeroonvac[tmp];
 If[tmp =!= 0,
  tmp = tmp / Sqrt[simpln[scalarproductop[tmp,tmp]]],
  0
 ]
];

(* Problem dependent utility functions *)

(* This will apply to f on site N+1 (that's being added), 
therefore n -> (NN+1). *)

fixisospin1[i_, iz_, op_] := Nest[isospin1down, op, i-iz];
fixisospin2[i_, iz_, op_] := Nest[isospin2down, op, i-iz];

fixstate[i1_, iz1_, i2_, iz2_, state_] :=
  (simpln @ fixisospin1[i1, iz1, fixisospin2[i2, iz2, state]]) /. n -> NN+1;

Invar[{diffi1_, diffi2_, diffsz_}] := 
"Invar(" <> ToString[2 diffi1] <> ", " <> ToString[2 diffi2] <> 
", " <> ToString[2 diffsz] <> ");";

InvarQN[{i1_, i2_, sz_}] := 
"Invar(" <> ToString[2i1+1] <> ", " <> ToString[2i2+1] <> 
", " <> ToString[2 sz] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr //. rule, Iq1+Iz1 \[Element] Integers && 
                                Iq2+Iz2 \[Element] Integers &&
                                2Sz \[Element] Integers &&
                                NN \[Element] Integers] //. rule;

rules4 = 2Iq1 \[Element] Integers && Iq1 > 4 && Iq1+Iz1 \[Element] Integers &&
         2Iq2 \[Element] Integers && Iq2 > 4 && Iq2+Iz2 \[Element] Integers &&
         2Sz \[Element] Integers; 

simpl4[expr_] := ratsimpl[ simplNN @ FullSimplify[expr, rules4], rules4];

simpl5[expr_] := simpl4[expr /. {Iz1->Iq1, Iz2->Iq2}]; 
(* Take maximal Iz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"Iq1"->"ISO(ii1)"] ;
  str = StringReplace[str,"Iq2"->"ISO(ii2)"] ;
  str
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"Iq1"->"ISO(ii1p)"] ;
  str = StringReplace[str,"Iq2"->"ISO(ii2p)"] ;
  str
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"Iq1"->"ISO(ii11)"];
  str = StringReplace[str,"Iq2"->"ISO(ii21)"];
  str 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

(* Recall: newstates2[qn,d,state] where qn are the preserved
quantum numbers, d the additional quantum numbers. see ../nrgcoef.m *)

newstates1[qn:{i1_, i2_, sz_}, state_] := 
  Do[newstates2[qn, {diffi1, diffi2, -sz}, state], {diffi1, -i1, i1}, {diffi2, -i2, i2}];

(* Produce the expression for the new state *)

getexpr1[qn:{i1_, i2_, sz_}, d:{diffi1_, diffi2_, diffsz_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{Iq1+diffi1, Iz1-iz1}, {i1, iz1}, {Iq1, Iz1}] *
    CG[{Iq2+diffi2, Iz2-iz2}, {i2, iz2}, {Iq2, Iz2}] ],
    dirpdt[{{Iq1+diffi1, Iq2+diffi2, Sz+diffsz}, {Iz1-iz1, Iz2-iz2}, r}, fixstate[i1, iz1, i2, iz2, state]] ],
  {iz1, -i1, i1}, {iz2, -i2, i2}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
(* Called from matrixel[]. *)
(* Result must be non-zero, since the term is found to be present !! *)
CGmatrixel[qn1:{i11_, i21_, sz1_}, diff1:{iz11_, iz21_}, 
           qn2:{i12_, i22_, sz2_}, diff2:{iz12_, iz22_},
           f[izn_, ch_, szn_]] := Module[{opiz, cg1, cg2, opsz, faktor},
   opiz = isoof[izn];
   cg1 = 1;
   cg2 = 1;
   If[ch == 0, (* cf. ch2op numbering! *)
     cg1 = CG[{i12, iz12}, {1/2, opiz}, {i11, iz11}];
   ];
   If[ch == 1,
     cg2 = CG[{i22, iz22}, {1/2, opiz}, {i21, iz21}];
   ];
   opsz = spinofop[izn, szn];
   faktor = If[opiz == 1/2, 1, (-1)^NN (2opsz)]; (* !!! Restored !!! *)
   faktor cg1 cg2 ];

(* irrule[] and ircg[] are called from recalcf[]. The results are then used
as parameters to ireducel[] in nrgcoef.m. irrule is used to prepare the 
bra when computing brakets. ircg is used as the denominator in the resulting
expression. *)

irrule[opiz1_, opiz2_, opsz_] = {
  Sz -> Sz + opsz,
  Iq1 -> Iq1 + opiz1, Iz1 -> Iz1 + opiz1,
  Iq2 -> Iq2 + opiz2, Iz2 -> Iz2 + opiz2
};

ircg[opiz1_, 0, opsz_] := 
  CG[{Iq1, Iz1}, {1/2, opiz1}, {Iq1+opiz1, Iz1+opiz1}] *
  If[opiz1 == 1/2, 1, (-1)^(NN+1) (2opsz)];

ircg[0, opiz2_, opsz_] :=
  CG[{Iq2, Iz2}, {1/2, opiz2}, {Iq2+opiz2, Iz2+opiz2}] *
  If[opiz2 == 1/2, 1, (-1)^(NN+1) (2opsz)];

ircg[a_, b_, c_] := (Print["ircg error:", a, " ",b, " ",c]; Exit[1]);

outmake[{a_, b_, c_}] := 
"Invar(" <> toc[2(a-Iq1)+ii11] <> ", " <> toc[2(b-Iq2)+ii21] <> 
", " <> toc[2(c-Sz)+ssz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i11_, i21_, sz1_}, diff1:{iz11_, iz21_}, 
            qn2:{i12_, i22_, sz2_}, diff2:{iz12_, iz22_}, 
            {mi1_, mi2_, ms_}] := Module[{opiz1, opiz2},
  opiz1 = iz11-iz12;
  opiz2 = iz21-iz22;
  opsz = sz1-sz2;
  CG[{i12, iz12}, {mi1, opiz1}, {i11, iz11}] *
  CG[{i22, iz22}, {mi2, opiz2}, {i21, iz21}] *
  If[mi == 1/2 && opiz1 == -1/2, 1 (2opsz), 1] *
  If[mi == 1/2 && opiz2 == -1/2, 1 (2opsz), 1]  
]; (* XXX: Are the two previous lines ok? *)

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai1_, deltai2_, deltasz_}] :=
  (Ip1 = Iq1+deltai1 ; Ip2 = Iq2+deltai2 ; Szp = Sz+deltasz);

RULE1recalcop[{mi1_, mi2_, ms_}] := {
 Iz1 -> Iq1, 
 Iz2 -> Iq2
};
RULE2recalcop[{mi1_, mi2_, ms_}] := {
 Iq1 -> Ip1, Iz1 -> Iq1-mi1,
 Iq2 -> Ip2, Iz2 -> Iq2-mi2,
 Sz -> Szp 
}; (* XXX: Is this correct??? *)

CG1recalcop[{mi1_, mi2_, ms_}] := 
  CG[{Ip1, Iq1-mi1}, {mi1, mi1}, {Iq1, Iq1}] *
  CG[{Ip2, Iq2-mi2}, {mi2, mi2}, {Iq2, Iq2}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely
symmetry reasons. *)

diffsD1 = {
  {{1/2, 0,  1/2},   FNDOUBLET <> "p0p.dat"},
  {{1/2, 0, -1/2},   FNDOUBLET <> "p0m.dat"},
  {{-1/2, 0,  1/2},  FNDOUBLET <> "m0p.dat"},
  {{-1/2, 0, -1/2},  FNDOUBLET <> "m0m.dat"}
};

diffsD2 = {
  {{0, 1/2,  1/2},   FNDOUBLET <> "0pp.dat"},
  {{0, 1/2, -1/2},   FNDOUBLET <> "0pm.dat"},
  {{0, -1/2,  1/2},  FNDOUBLET <> "0mp.dat"},
  {{0, -1/2, -1/2},  FNDOUBLET <> "0mm.dat"}
};

diffsS = {
  {{0, 0, 0}, FNSINGLET <> ".dat"} 
};

diffsT = {
  {{0, 0,  0}, FNTRIPLET <> "s.dat"},
  {{0, 0,  1}, FNTRIPLET <> "p.dat"},
  {{0, 0, -1}, FNTRIPLET <> "m.dat"}
};

checksignfn[{mi1_, mi2_, ms_}] := ! (mi1 \[Element] Integers && mi2 \[Element] Integers);

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

(* NOTE: matrixel[] calls CGmatirxel[]. *)

dorecalcfLOOP["DBLISOSZ"] := Module[{ch, op, nn},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];

    If[ch == 0,
      qn1 = {1/2, 0}; (* CR *)
      qn2 = {-1/2, 0}; (* AN *)
    ];
    If[ch == 1,
      qn1 = {0, 1/2};
      qn2 = {0, -1/2};
    ];

    (* First doublet: [f^dag_UP, f_DO] szop=1/2 pair *)
    makerecalcf[op[CR, UP], Append[qn1, 1/2], PREFIX <> "-type1-isoup-" <> tos[op]]; (* type1 = spinup *)
    makerecalcf[op[AN, DO], Append[qn2, 1/2], PREFIX <> "-type1-isodown-" <> tos[op]]; (* type1 = spinup *)

    (* Second doublet: [f^dag_DO, -f_UP] szop=-1/2 pair *)
    makerecalcf[op[CR, DO],  Append[qn1, -1/2], PREFIX <> "-type2-isoup-" <> tos[op]]; (* type2 = spindo *)
    makerecalcf[op[AN, UP], Append[qn2, -1/2], PREFIX <> "-type2-isodown-" <> tos[op]]; (* type2 = spindo *)
    (* XXX: NOTE: in DBLSU2, we have -op[AN,UP] !! *)
  ];
];

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},  
 op = ch2op[ch];
 opspinz = spinz[ op[] ];
 matrixop1[i1, i2, opspinz]
];

(* ADDED 16.12.2009 *)
orthogonalitytest[] := Module[{res},
  Print["orthogonalitytest[]"];
  For[in1 = 1, in1 <= nrstates, in1++,
    For[in2 = 1, in2 <= nrstates, in2++,
      res = matrixop1[in1, in2, 1];
      If[res =!= 0,
        Print[in1, " ", in2, " ", res];
      ];
    ];
  ];
  Print["DONE"];
];

(* NOTE: we also enforce agreement in Iz of the other element (i.e. the one
not related to the operator f being added. Therefore we need to compare
newd[], which contains the quantum number differences. *)

MATRIXNICHECKSYM[i1_, i2_, _] := ! (
  ((Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] == 1/2) && (Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] == 0) &&
  (Abs[ newd [i1] [[1]] - newd [i2] [[1]] ] == 1/2) && (Abs[ newd [i1] [[2]] - newd [i2] [[2]] ] == 0))
  ~ Xor ~
  ((Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] == 1/2) && (Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] == 0) &&
  (Abs[ newd [i1] [[2]] - newd [i2] [[2]] ] == 1/2) && (Abs[ newd [i1] [[1]] - newd [i2] [[1]] ] == 0))
);

doall[] := Module[{},
  donew[];
  orthogonalitytest[];
  dooffdiagNEW[];
  dorecalcf["DBLISOSZ"];

  (* ATTENTION: upper=False here! These are not used to construct the Hamiltonian... *)
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-H-spinz.dat", True];

  recalcop[{1/2, 0, 1/2}, diffsD1];
  recalcop[{0, 1/2, 1/2}, diffsD2];
  recalcop[{0, 0, 0}, diffsS];
  recalcop[{0, 0, 1}, diffsT];
];
