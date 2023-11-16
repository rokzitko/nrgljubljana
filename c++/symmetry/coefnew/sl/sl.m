(*
SYMTYPE=SL
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "sl-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = qbasis[ Take[allops, channels] ];
basis = spinlessbasis[basis];
mult[{q_}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffq_}] :=
  "Invar(" <> ToString[diffq] <> ");";

InvarQN[{q_}] := 
  "Invar(" <> ToString[q] <> ");";

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

newstates1[qn:{q_}, state_] := newstates2[qn, {-q}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{q_}, d:{diffq_}, state_] := 
  nc[coef[1], dirpdt[{{Q+diffq}, {}, r}, state]];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely 
   symmetry reasons. *)
MATRIXNICHECKSYM[i1_, i2_, _] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1;

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{q1_}, diff1:{}, qn2:{q2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)

dorecalcfLOOP["SL"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {}, FNRECALCF <> "-" <> tos[op] ];
  ];
];

irrule[] = {Q -> Q+1};
ircg[] = 1;

outmake[{a_}] := 
  "Invar(" <> toc[a /. Q -> q1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q1_}, diff1:{},
            qn2:{q2_}, diff2:{}, 
            {}] := 1;
 
QNrecalcop[delta:{deltaq_}] :=
  (Qp = Q+deltaq);

RULE1recalcop[{}] := {};
RULE2recalcop[{}] := {Q -> Qp};

CG1recalcop[{}] := 1;

diffsS = { {{0}, FNSINGLET<>".dat"} };

diffsD = {
  {{-1}, FNDOUBLET<>".dat"}
};

matrixopqdiff[i1_, i2_, ch_] := Module[{op, opqdiff},
  op = ch2op[ch];
  opqdiff = number[ op[] ]  If[ch == 1, +1, -1, MyError[]];
  matrixop1[i1, i2, opqdiff]
];

matrixopqtot[i1_, i2_, ch_] := Module[{op, opqtot},
  op = ch2op[ch];
  opqtot = number[ op[] ] - 1; (* Subtract 1 to conform to the definition of Q in the NRG code! *)
  matrixop1[i1, i2, opqtot]
];

matrixopN1[i1_, i2_, ch_] := Module[{op, opqtot},
  op = ch2op[0];
  opqtot = number[ op[] ] - 1; (* Subtract 1 to conform to the definition of Q in the NRG code! *)
  matrixop1[i1, i2, opqtot]
];

matrixopN2[i1_, i2_, ch_] := Module[{op, opqtot},
  op = ch2op[1];
  opqtot = number[ op[] ] - 1; (* Subtract 1 to conform to the definition of Q in the NRG code! *)
  matrixop1[i1, i2, opqtot]
];

matrixopN3[i1_, i2_, ch_] := Module[{op, opqtot},
  If[ch <= -1 || ch >= 3, Exit[]];
  op = ch2op[2];
  opqtot = number[ op[] ] - 1; (* Subtract 1 to conform to the definition of Q in the NRG code! *)
  matrixop1[i1, i2, opqtot]
];

doall[] := Module[{},
  donew[];
  dooffdiagNEW[]; (* No change required: spin-down terms will just drop out! *)
  dodiag[];
  dorecalcf["SL"];
  (* Global spin-difference operator for two-channel problems *)
  If[channels == 2,
      dogeneralloop[matrixopqdiff, "QDIFF", PREFIX <> "-qdiff.dat", True];
      dogeneralloop[matrixopqtot,  "QTOT",  PREFIX <> "-qtot.dat",  True];
  ];
  If[channels == 3,
      dogeneralloop[matrixopN1, "N1", PREFIX <> "-N1.dat", True];
      dogeneralloop[matrixopN2, "N2", PREFIX <> "-N2.dat", True];
      dogeneralloop[matrixopN3, "N3", PREFIX <> "-N3.dat", True];
  ];
  checksignfn[{}] := True;
  recalcop[{}, diffsD];
  checksignfn[{}] := False;
  recalcop[{}, diffsS];
];
