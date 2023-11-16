(* 
SYMTYPE=DBLSU2
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2009-2023
*)

PREFIX = "dblsu2-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
Print["basisops=", basisops];

basisops1 = {First[basisops]};
basisops2 = {Last[basisops]};
basis = quickDBL[basisops1, basisops2, quickSU2basis];
Print["basis=", basis];

snegrealconstants[n, NN];

(* IMPORTANT: the following choice of nnop corresponds
to a particular partition of the lattice sites into two
sets! *)

nnop[a[]] = n;
nnop[b[]] = n;

Tminus1 = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops1];
Print["Tminus1=", Tminus1];
Tminus2 = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops2];
Print["Tminus2=", Tminus2];

mult[{i1_, i2_}] := (2i1+1)(2i2+1);

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

Invar[{diffi1_, diffi2_}] := 
  "Invar(" <> ToString[2 diffi1] <> ", " <> ToString[2 diffi2] <> ");";

InvarQN[{i1_, i2_}] := 
  "Invar(" <> ToString[2i1+1] <> ", " <> ToString[2i2+1] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule, Iq1+Iz1 \[Element] Integers && 
                                Iq2+Iz2 \[Element] Integers &&
                                NN \[Element] Integers];

rules4 = 2Iq1 \[Element] Integers && Iq1 > 4 && Iq1+Iz1 \[Element] Integers &&
         2Iq2 \[Element] Integers && Iq2 > 4 && Iq2+Iz2 \[Element] Integers;

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

newstates1[qn:{i1_, i2_}, state_] := 
  Do[newstates2[qn, {diffi1, diffi2}, state], {diffi1, -i1, i1}, {diffi2, -i2, i2}];

(* Produce the expression for the new state *)

getexpr1[qn:{i1_, i2_}, d:{diffi1_, diffi2_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{Iq1+diffi1, Iz1-iz1}, {i1, iz1}, {Iq1, Iz1}] *
    CG[{Iq2+diffi2, Iz2-iz2}, {i2, iz2}, {Iq2, Iz2}] ],
    dirpdt[{{Iq1+diffi1, Iq2+diffi2}, {Iz1-iz1, Iz2-iz2}, r}, fixstate[i1, iz1, i2, iz2, state]] ],
  {iz1, -i1, i1}, {iz2, -i2, i2}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
(* Called from matrixel[]. *)
(* Result must be non-zero, since the term is found to be present !! *)
CGmatrixel[qn1:{i11_, i21_}, diff1:{iz11_, iz21_}, 
           qn2:{i12_, i22_}, diff2:{iz12_, iz22_},
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

irrule[opiz1_, opiz2_] = {
  Iq1 -> Iq1 + opiz1, Iz1 -> Iz1 + opiz1,
  Iq2 -> Iq2 + opiz2, Iz2 -> Iz2 + opiz2
};

ircg[opiz1_, 0] := 
  CG[{Iq1, Iz1}, {1/2, opiz1}, {Iq1+opiz1, Iz1+opiz1}] *
  If[opiz1 == 1/2, 1, (-1)^(NN+1)];

ircg[0, opiz2_] := 
  CG[{Iq2, Iz2}, {1/2, opiz2}, {Iq2+opiz2, Iz2+opiz2}] *
  If[opiz2 == 1/2, 1, (-1)^(NN+1)];

ircg[a_, b_] := (Print["ircg error:", a, " ",b]; Exit[1]);

outmake[{a_, b_}] := 
  "Invar(" <> toc[2(a-Iq1)+ii11] <> ", " <> toc[2(b-Iq2)+ii21] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i11_, i21_}, diff1:{iz11_, iz21_}, 
            qn2:{i12_, i22_}, diff2:{iz12_, iz22_}, 
            {mi1_, mi2_}] := Module[{opiz1, opiz2},
  opiz1 = iz11-iz12;
  opiz2 = iz21-iz22;
  CG[{i12, iz12}, {mi1, opiz1}, {i11, iz11}] *
  CG[{i22, iz22}, {mi2, opiz2}, {i21, iz21}]
];

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai1_, deltai2_}] :=
  (Ip1 = Iq1+deltai1 ; Ip2 = Iq2+deltai2);

RULE1recalcop[{mi1_, mi2_}] := {
 Iz1 -> Iq1, 
 Iz2 -> Iq2
};
RULE2recalcop[{mi1_, mi2_}] := {
 Iq1 -> Ip1, Iz1 -> Iq1-mi1,
 Iq2 -> Ip2, Iz2 -> Iq2-mi2
};

CG1recalcop[{mi1_, mi2_}] := 
  CG[{Ip1, Iq1-mi1}, {mi1, mi1}, {Iq1, Iq1}] *
  CG[{Ip2, Iq2-mi2}, {mi2, mi2}, {Iq2, Iq2}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely
symmetry reasons. *)

(* This HAS TO BE USED! In this way we avoid generating terms which appear
   because we are sloppy and do not explicitly attach a "type" label
   to the f operators. *)

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

diffsD1 = {
  {{1/2, 0},   FNDOUBLET <> "p0.dat"},
  {{-1/2, 0},  FNDOUBLET <> "m0.dat"}
};

diffsD2 = {
  {{0, 1/2},   FNDOUBLET <> "0p.dat"},
  {{0, -1/2},  FNDOUBLET <> "0m.dat"}
};

diffsS = {
  {{0, 0}, FNSINGLET <> ".dat"} 
};

checksignfn[{mi1_, mi2_}] := ! (mi1 \[Element] Integers && mi2 \[Element] Integers);

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

(* NOTE: matrixel[] calls CGmatirxel[]. *)

(* Combination 1: related to type 1 doublet *)
matrixnihop1[i1_, i2_, ch_] := Module[{op, expr = 0},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];
  Print[newqn[i1]," ",newqn[i2]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR, ch, UP], op[AN, UP]]];
  expr += matrixel[i1, i2, nc[op[CR, DO], f[AN, ch, DO]]];
  simpl4 @ expr
];

(* Combination 2: related to type 2 doublet *)
matrixnihop2[i1_, i2_, ch_] := Module[{op, expr = 0},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];
  Print[newqn[i1]," ",newqn[i2]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR, ch, DO], op[AN, DO]]];
  expr += matrixel[i1, i2, nc[op[CR, UP], f[AN, ch, UP]]];
  simpl4 @ expr
];

dorecalcfLOOP["DBLSU2"] := Module[{ch, op, nn},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];

    (* makerecalcf[] -- First argument: operator. Second argument: izop. *)

    If[ch == 0,
      qn1 = {1/2, 0};
      qn2 = {-1/2, 0};
    ];
    If[ch == 1, 
      qn1 = {0, 1/2};
      qn2 = {0, -1/2};
    ];

    (* First doublet: [f^dag_UP, f_DO] szop=1/2 pair *)
    makerecalcf[op[CR, UP], qn1, PREFIX <> "-type1-isoup-" <> tos[op]];
    makerecalcf[op[AN, DO], qn2, PREFIX <> "-type1-isodown-" <> tos[op]];

    (* Second doublet: [f^dag_DO, -f_UP] szop=-1/2 pair *)
    makerecalcf[op[CR, DO],  qn1, PREFIX <> "-type2-isoup-" <> tos[op]];
    makerecalcf[-op[AN, UP], qn2, PREFIX <> "-type2-isodown-" <> tos[op]];
  ];
];

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},  
 op = ch2op[ch];
 opspinz = spinz[ op[] ];
 matrixop1[i1, i2, opspinz]
];

matrixopspiny[i1_, i2_, ch_] := Module[{op, opspiny},
 op = ch2op[ch];
 opspiny = spiny[ op[] ];
 matrixop1[i1, i2, opspiny]
];

matrixopspinx[i1_, i2_, ch_] := Module[{op, opspinx},
 op = ch2op[ch];
 opspinx = spinx[ op[] ];
 matrixop1[i1, i2, opspinx]
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

doall[] := Module[{},
  donew[];
  orthogonalitytest[];
  dogeneralloop[matrixnihop1, "OFFDIAG_1", PREFIX <> "-offdiag-1.dat", True];
  dogeneralloop[matrixnihop2, "OFFDIAG_2", PREFIX <> "-offdiag-2.dat", True];
  dorecalcf["DBLSU2"];

  (* ATTENTION: upper=False here! These are not used to construct
  the Hamiltonian... *)
  dogeneralloop[matrixopspinx, "SPINX", PREFIX <> "-spinx.dat", False];
  dogeneralloop[matrixopspiny, "SPINY", PREFIX <> "-spiny.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];

  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-H-spinz.dat", True];

  recalcop[{1/2, 0}, diffsD1];
  recalcop[{0, 1/2}, diffsD2];
  recalcop[{0, 0}, diffsS];
];
