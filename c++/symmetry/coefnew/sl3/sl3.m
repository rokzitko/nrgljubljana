(*
SYMTYPE=SL3
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2010-2023
*)

If[channels != 3,
  MyError["3 channels required!"];
];

PREFIX = "sl3-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
op1 = basisops[[1]];
op2 = basisops[[2]];
op3 = basisops[[3]];

bz1 = qbasis[{op1}];
bz2 = qbasis[{op2}];
bz3 = qbasis[{op3}];

bz12 = basistensorproduct[bz1, bz2, {#1[[1]], #2[[1]]} &];
basis = basistensorproduct[bz12, bz3, {#1[[1]], #1[[2]], #2[[1]]} &];

basis = spinlessbasis[basis];
mult[{q1_,q2_,q3_}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffq1_, diffq2_, diffq3_}] :=
"Invar(" <> 
ToString[diffq1] <> ", " <> 
ToString[diffq2] <> ", " <> 
ToString[diffq3] <> ");";

InvarQN[{q1_, q2_, q3_}] := 
"Invar(" <> 
ToString[q1] <> ", " <> 
ToString[q2] <> ", " <> 
ToString[q3] <> ");";

(* For expressions with S and Sz. *)
simpl2[expr_] := Simplify[expr];

(* For expressions with S only. *)
simpl4[expr_] := FullSimplify[expr];
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

newstates1[qn:{q1_, q2_, q3_}, state_] := newstates2[qn, {-q1, -q2, -q3}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{q1_, q2_, q3_}, d:{diffq1_, diffq2_, diffq3_}, state_] := 
  nc[coef[1], dirpdt[{{Q1+diffq1, Q2+diffq2, Q3+diffq3}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
(* First index: channel; Second index: which of the two *)
CGmatrixel[qn1:{q11_, q21_, q31_}, diff1:{}, 
           qn2:{q12_, q22_, q32_}, diff2:{},
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)

dorecalcfLOOP["SL3"] := Module[{ch, op, opq1, opq2, opq3},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    opq1 = opq2 = opq3 = 0;
    If[ch == 0, opq1 = 1];
    If[ch == 1, opq2 = 1];
    If[ch == 2, opq3 = 1];
    makerecalcf[op[CR, UP], {opq1, opq2, opq3}, FNRECALCF <> "-" <> tos[op] ];
  ];
];

irrule[opq1_, opq2_, opq3_] = {Q1 -> Q1+opq1, Q2 -> Q2+opq2, Q3 -> Q3+opq3};
ircg[opq1_, opq2_, opq3_] = 1;

outmake[{a_, b_, c_}] :=
"Invar(" <> toc[a /. Q1 -> q11] <> ", " <>
            toc[b /. Q2 -> q21] <> ", " <>
            toc[c /. Q3 -> q31] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q11_, q21_, q31_}, diff1:{},
            qn2:{q12_, q22_, q32_}, diff2:{}, 
            {}] := 1;

QNrecalcop[delta:{deltaq1_, deltaq2_, deltaq3_}] :=
  (Q1p = Q1+deltaq1; Q2p = Q2+deltaq2; Q3p = Q3+deltaq3);

RULE1recalcop[{}] := {};
RULE2recalcop[{}] := {Q1 -> Q1p, Q2 -> Q2p, Q3 -> Q3p};

CG1recalcop[{}] := 1;

diffsD = {
  {{-1, 0, 0}, FNDOUBLET<>".dat"}
};

diffsS = { {{0, 0, 0}, FNSINGLET<>".dat"} };

matrixopqtot[i1_, i2_, ch_] := Module[{op, opqtot},
  op = ch2op[ch];
  opqtot = number[ op[] ] - 3;
  matrixop1[i1, i2, opqtot]
];

matrixopN1[i1_, i2_, ch_] := Module[{op, opqtot},
  If[ch != 0, Return[0]];
  op = ch2op[0];
  opqtot = number[ op[] ] - 1; 
  matrixop1[i1, i2, opqtot]
];

matrixopN2[i1_, i2_, ch_] := Module[{op, opqtot},
  If[ch != 1, Return[0]];
  op = ch2op[1];
  opqtot = number[ op[] ] - 1; 
  matrixop1[i1, i2, opqtot]
];

matrixopN3[i1_, i2_, ch_] := Module[{op, opqtot},
  If[ch != 2, Return[0]];
  op = ch2op[2];
  opqtot = number[ op[] ] - 1;
  matrixop1[i1, i2, opqtot]
];

doall[] := Module[{},
  donew[];
  dooffdiagNEW[]; (* No change required: spin-down terms will just drop out! *)
  dodiag[];
  dorecalcf["SL3"];
  dogeneralloop[matrixopqtot,  "QTOT",  PREFIX <> "-qtot.dat",  True];
  dogeneralloop[matrixopN1, "N1", PREFIX <> "-N1.dat", True];
  dogeneralloop[matrixopN2, "N2", PREFIX <> "-N2.dat", True];
  dogeneralloop[matrixopN3, "N3", PREFIX <> "-N3.dat", True];
  checksignfn[{}] := True; (* Always proceeds *)
  recalcop[{}, diffsD];
  checksignfn[{}] := False;
  recalcop[{}, diffsS];
];
