(*
SYMTYPE=DBLQSZ  U(1)_{charge} x U(1)_{charge} x U(1)_{spin}
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2022-2023
*)

debug=10;

PREFIX = "dblqsz-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
Print["basisops=", basisops];

basisops1 = {First[basisops]};
basisops2 = {Last[basisops]};
basis = quickDBLSZ[basisops1, basisops2, qszbasis];
Print["basis=", basis];

mult[{_, _, _}] := 1;

basisprep[];

(* Problem dependent utility functions *)

Invar[{diffq1_, diffq2_, diffsz_}] := 
"Invar(" <> ToString[diffq1] <> ", " <> ToString[diffq2] <>
", " <> ToString[2 diffsz] <> ");";

InvarQN[{q1_, q2_, sz_}] := 
"Invar(" <> ToString[q1] <> ", " <> ToString[q2] <>
", " <> ToString[2 sz] <> ");";

simpl2[expr_] := Simplify[expr, 2Sz \[Element] Integers];

rules4 = 2Sz \[Element] Integers;

simpl4[expr_] := FullSimplify[expr, rules4];

simpl5[expr_] := simpl4[expr];

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

(* Recall: newstates2[qn,d,state] where qn are the preserved
quantum numbers, d the additional quantum numbers. see ../nrgcoef.m *)

newstates1[qn:{q1_, q2_, sz_}, state_] := newstates2[qn, {-q1, -q2, -sz}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{q1_, q2_, sz_}, d:{diffq1_, diffq2_, diffsz_}, state_] := 
  nc[coef[1], dirpdt[{{Q1+diffq1, Q2+diffq2, Sz+diffsz}, {}, r}, state]];

CGmatrixel[qn1:{q11_, q21_, sz1_}, diff1:{}, qn2:{q12_, q22_, sz2_}, diff2:{},
  f[izn_, _, szn_]] := 1;

(* irrule[] and ircg[] are called from recalcf[]. The results are then used
as parameters to ireducel[] in nrgcoef.m. irrule is used to prepare the 
bra when computing brakets. ircg is used as the denominator in the resulting
expression. *)

irrule[opq1_, opq2_, opsz_] = { Q1 -> Q1+opq1, Q2 -> Q2+opq2, Sz -> Sz+opsz };
ircg[opq1_, opq2_, opsz_] = 1;

outmake[{a_, b_, c_}] :=
"Invar(" <> toc[a /. Q1->q11] <> ", " <> toc[b /. Q2->q21] <>
", " <> toc[2(c-Sz)+ssz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[___] := 1;

(* Left is (Q1, Q2, Sz), right is (Q1', Q2', Sz')=(Q1p, Q2p, Szp). *)
QNrecalcop[delta:{deltaq1_, deltaq2_, deltasz_}] :=
  (Q1p = Q1+deltaq1 ; Q2p = Q2+deltaq2 ; Szp = Sz+deltasz);

RULE1recalcop[{__}] := {};
(* Just renamed according to QNrecalcop[] definitions! *)
RULE2recalcop[{mq1_, mq2_, ms_}] := { Q1 -> Q1p, Q2 -> Q2p,  Sz -> Szp };

CG1recalcop[{mq1_, mq2_, ms_}] := 1;

diffsD1 = {
  {{-1, 0,  1/2},  FNDOUBLET <> "m0p.dat"},
  {{-1, 0, -1/2},  FNDOUBLET <> "m0m.dat"}
};

diffsD2 = {
  {{0, -1,  1/2},  FNDOUBLET <> "0mp.dat"},
  {{0, -1, -1/2},  FNDOUBLET <> "0mm.dat"}
};

diffsS = {
  {{0, 0, 0}, FNSINGLET <> ".dat"}
};

diffsT = {
  {{0, 0,  0}, FNTRIPLET <> "s.dat"},
  {{0, 0,  1}, FNTRIPLET <> "p.dat"},
  {{0, 0, -1}, FNTRIPLET <> "m.dat"}
};

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

(* NOTE: matrixel[] calls CGmatirxel[]. *)

dorecalcfLOOP["DBLQSZ"] := Module[{ch, op, qn},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];

    If[ch == 0,
      qn = {1, 0}; (* CR *)
    ];
    If[ch == 1,
      qn = {0, 1};
    ];

    makerecalcf[op[CR, UP], Append[qn, 1/2],  PREFIX <> "-spinup-" <> tos[op]];
    makerecalcf[op[CR, DO], Append[qn, -1/2], PREFIX <> "-spindown-" <> tos[op]];
  ];
];

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
 op = ch2op[ch];
 opspinz = spinz[ op[] ];
 matrixop1[i1, i2, opspinz]
];

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

checksignfn[{_, _, ms_}] := ! (ms \[Element] Integers);

doall[] := Module[{},
  donew[];
  orthogonalitytest[];
  dooffdiagNEW[];
  dodiag[];
  dorecalcf["DBLQSZ"];

  (* ATTENTION: upper=False here! These are *not* used to construct
  the Hamiltonian... *)
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-H-spinz.dat", True];

  recalcop[{-1, 0, 1/2}, diffsD1];
  recalcop[{0, -1, 1/2}, diffsD2];
  recalcop[{0, 0, 0}, diffsS];
  recalcop[{0, 0, 1}, diffsT];
];
