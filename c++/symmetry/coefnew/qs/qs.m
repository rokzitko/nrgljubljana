(*
SYMTYPE=QS
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

(* CHANGE LOG *)
(* 15.12.2010 - Q1, Q2 *)

symtype = "QS";
PREFIX = "qs-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = qsbasis[ Take[allops, channels] ];
mult[{q_, s_}] := 2s+1;
basisprep[];

(* Problem dependent utility functions *)
fixstate[s_, sz_, state_] := fixspin[s, sz, state];

Invar[{diffq_, diffs_}] := 
  "Invar(" <> ToString[diffq] <> ", " <> ToString[2diffs] <> ");";

InvarQN[{q_, s_}] := 
  "Invar(" <> ToString[q] <> ", " <> ToString[2s+1] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

(* For expressions with S and Sz. *)
simpl2[expr_] := Simplify[expr /. rule, S+Sz \[Element] Integers && S - Sz \[Element] Integers &&
                                2S  \[Element] Integers && S >= 4 &&    (* WAS S >= 0 *)
                                2Sz \[Element] Integers && S >= Sz && S+Sz > 2 && S > Sz+2];

(* For expressions with S only. *)
simpl4[expr_] := FullSimplify[expr, 2S \[Element] Integers && S >= 4];  (* WAS S >=- 0 *)
simpl5[expr_] := simpl4[expr /. Sz->S]; (* Take maximal Sz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss)"] 
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"S"->"S(ssp)"] 
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss1)"] 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

newstates1[qn:{q_, s_}, state_] := 
  Do[newstates2[qn, {-q, diffs}, state], {diffs, -s, s}];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, s_}, d:{diffq_, diffs_}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}],
    dirpdt[{{Q+diffq, S+diffs}, {Sz-sz}, r}, fixstate[s, sz, state]]], 
    {sz, -s, s}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely 
   symmetry reasons. *)
MATRIXNICHECKSYM[i1_, i2_, _] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1 ||
                                 Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] != 1/2;   

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{q1_, s1_}, diff1:{sz1_}, qn2:{q2_, s2_}, diff2:{sz2_}, 
  f[izn_, _, szn_]] := CG[{s2, sz2}, {1/2, spinofop[izn,szn]}, {s1, sz1}];

(* recalcf[] stuff *)
irrule[opsz_] = {S -> S + opsz, Sz -> Sz + opsz, Q -> Q+1};
ircg[opsz_] = CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

outmake[{a_, b_}] := 
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-S)+ss1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q1_, s1_}, diff1:{sz1_}, 
            qn2:{q2_, s2_}, diff2:{sz2_}, 
            {m_}] :=
  CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];

QNrecalcop[delta:{deltaq_, deltas_}] := 
  (Qp = Q+deltaq; Sp = S+deltas);

RULE1recalcop[{m_}] := {Sz -> S};
RULE2recalcop[{m_}] := {Q -> Qp, S -> Sp, Sz -> S-m};

(* Prime: bra. *)
(* QS se nanasajo na bra. *)
(* Ket in op se sklopita v bra. *)
CG1recalcop[{m_}] := CG[{Sp, S-m}, {m, m}, {S, S}];

diffsD = {
  {{-1, 1/2}, FNDOUBLET<>"p.dat"},
  {{-1, -1/2}, FNDOUBLET<>"m.dat"}
};

diffsS = { {{0, 0}, FNSINGLET<>".dat"} };

diffsT = { 
  {{0, 0}, FNTRIPLET<>"s.dat"}, 
  {{0, 1}, FNTRIPLET<>"p.dat"},
  {{0, -1}, FNTRIPLET<>"m.dat"}
};
  
matrixopqdiff[i1_, i2_, ch_] := Module[{op, opqdiff},
  op = ch2op[ch];
  opqdiff = number[ op[] ]  If[ch == 1, +1, -1, MyError[]];
  matrixop1[i1, i2, opqdiff]
];

matrixopq1[i1_, i2_, ch_] := Module[{op, opq1},
  op = ch2op[ch];
  opq1 = number[ op[] ];
  If[ch == 0, matrixop1[i1, i2, opq1], 0]
];

matrixopq2[i1_, i2_, ch_] := Module[{op, opq2},
  op = ch2op[ch];
  opq2 = number[ op[] ];
  If[ch == 1, matrixop1[i1, i2, opq2], 0]
];

matrixopqtot[i1_, i2_, ch_] := Module[{op, opqtot},
  op = ch2op[ch];
  opqtot = number[ op[] ] - 1; (* Subtract 1 to conform to the definition of Q in the NRG code! *)
  matrixop1[i1, i2, opqtot]
];

do2x2loop[function_, prefix_, filename_, upper_:False, INCLUDEZEROONDIAG_:False] := Module[
  {code, in1, startin2, in2, mx, koef, string},
  MyPrint["do2x2loop[", function, " ", prefix, " ", filename, "]"];
  code = {};
  For[in1 = 1, in1 <= nrstates, in1++,
    startin2 = If[upper, in1, 1];
    For[in2 = startin2, in2 <= nrstates, in2++,
      MyPrint[2, in1, " ", in2];
      mx = function[in1, in2];
      If[mx =!= 0 || (INCLUDEZEROONDIAG && in1 == in2),
        koef = koefstr1[mx];
        string = prefix <> "(" <>
          tos[in1] <> ", " <> tos[in2] <> ", " <>
          koef <> ");";
        Print[string];
        AppendTo[code,string];
      ];
    ];
  ];
  store[code, filename];
];

matrixrunghop[i1_, i2_] := Module[{op1, op2, op},
  op1 = ch2op[0];
  op2 = ch2op[1];
  op = hop[ op1[], op2[]];
  matrixop1[i1, i2, op]
];

matrixnimix[i1_, i2_, ch_] := Module[{expr=0, op1},
  op1 = ch2op[1-ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op1[AN,UP]] ];
  expr += matrixel[i1, i2, nc[f[CR,ch,DO], op1[AN,DO]] ];
  res = simpl4 @ expr;
  Print[res];
  res
];

doall[] := Module[{},
  donew[];
  dogeneralloop[matrixnimix,      "OFFDIAG_MIX",       PREFIX <> "-offdiag-mix.dat"];
  do2x2loop[matrixrunghop,     "RUNGHOP",      PREFIX <> "-runghop.dat"];
  dooffdiagNEW[];
  dodiag[];
  dorecalcf[];

  (* Global spin-difference operator for two-channel problems *)
  If[channels == 2,
    dogeneralloop[matrixopqdiff, "QDIFF", PREFIX <> "-qdiff.dat", True];
    dogeneralloop[matrixopq1,    "Q1",    PREFIX <> "-q1.dat",    True];
    dogeneralloop[matrixopq2,    "Q2",    PREFIX <> "-q2.dat",    True];
    dogeneralloop[matrixopqtot,  "QTOT",  PREFIX <> "-qtot.dat",  True];
    dogeneralloop[matrixnimix,      "OFFDIAG_MIX",       PREFIX <> "-offdiag-mix.dat"];
    do2x2loop[matrixrunghop,     "RUNGHOP",      PREFIX <> "-runghop.dat"];
  ];

  recalcop[{1/2}, diffsD];
  recalcop[{0}, diffsS];
  recalcop[{1}, diffsT];
];
