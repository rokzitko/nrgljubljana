(* 
SYMTYPE=QSZ
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

debug=1;

PREFIX = "qsz-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = qszbasis[ Take[allops, channels] ];
mult[{__}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffq_, diffsz_}] := 
  "Invar(" <> ToString[diffq] <> ", " <> ToString[2diffsz] <> ");";

InvarQN[{q_, sz_}] := 
  "Invar(" <> ToString[q] <> ", " <> ToString[2sz] <> ");";

simpl2[expr_] := Simplify[expr, 2Sz \[Element] Integers];

simpl4[expr_] := FullSimplify[expr, 2Sz \[Element] Integers];
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

newstates1[qn:{q_, sz_}, state_] := 
  newstates2[qn, {-q, -sz}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, sz_}, d:{diffq_, diffsz_}, state_] := 
  nc[coef[1], dirpdt[{{Q+diffq, Sz+diffsz}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{q1_, sz1_}, diff1:{}, qn2:{q2_, sz2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)
irrule[opsz_] = {Sz -> Sz + opsz, Q -> Q+1};
ircg[opsz_] = 1;

outmake[{a_, b_}] := 
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-Sz)+ssz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q1_, sz1_}, diff1:{}, 
            qn2:{q2_, sz2_}, diff2:{}, 
            {m_}] := 1;

(* Left is (Q,Sz), right is (Q', Sz')=(Qp, Szp). *)
QNrecalcop[delta:{deltaq_, deltasz_}] := 
  (Qp = Q+deltaq; Szp = Sz+deltasz);

RULE1recalcop[{m_}] := {};

(* Just renamed according to QNrecalcop[] definitions! *)
RULE2recalcop[{m_}] := {Q -> Qp, Sz -> Szp};

CG1recalcop[{m_}] := 1;

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

matrixopnumberup[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opnumberup = number[ op[], UP];
  matrixop1[i1, i2, opnumberup]
];

matrixopnumberdown[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opnumberdown = number[ op[], DO];
  matrixop1[i1, i2, opnumberdown]
];

diffsS = { {{0, 0}, FNSINGLET<>".dat"} };

diffsD = {
  {{-1, 1/2}, FNDOUBLET<>"p.dat"},
  {{-1, -1/2}, FNDOUBLET<>"m.dat"}
};

diffsT = { 
  {{0, 0}, FNTRIPLET<>"s.dat"}, 
  {{0, 1}, FNTRIPLET<>"p.dat"},
  {{0, -1}, FNTRIPLET<>"m.dat"}
};  

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

matrixrunghopup[i1_, i2_] := Module[{op1, op2, op},
  op1 = ch2op[0];
  op2 = ch2op[1];
  op = hop[ op1[], op2[], UP];
  matrixop1[i1, i2, op]
];

matrixrunghopdown[i1_, i2_] := Module[{op1, op2, op},
  op1 = ch2op[0];
  op2 = ch2op[1];
  op = hop[ op1[], op2[], DO];
  matrixop1[i1, i2, op]
];

(* ch is the orbital/channel index of the operator on the previous site (f).
   op1 corresponds to the operator on the added site. *)

matrixnimix[i1_, i2_, ch_] := Module[{expr=0, op1},
  op1 = ch2op[1-ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op1[AN,UP]] ];
  expr += matrixel[i1, i2, nc[f[CR,ch,DO], op1[AN,DO]] ];
  res = simpl4 @ expr;
  Print[res];
  res
];

matrixnimixup[i1_, i2_, ch_] := Module[{expr=0, op1},
  op1 = ch2op[1-ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op1[AN,UP]] ];
  simpl4 @ expr
];

matrixnimixdown[i1_, i2_, ch_] := Module[{expr=0, op1},
  op1 = ch2op[1-ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,DO], op1[AN,DO]] ];
  simpl4 @ expr
];


matrixopn1u[i1_, i2_, ch_] := Module[{op, opspinz},
  If[ch == 0,
  op = ch2op[ch];
  opn = number[ op[], UP ];
  matrixop1[i1, i2, opn],
  0]
];

matrixopn1d[i1_, i2_, ch_] := Module[{op, opspinz},
  If[ch == 0,
  op = ch2op[ch];
  opn = number[ op[], DO ];
  matrixop1[i1, i2, opn],
  0]
];

matrixopq1u[i1_, i2_, ch_] := Module[{op, opspinz},
  If[ch == 0,
  op = ch2op[ch];
  opn = number[ op[], UP ]-1/2;
  matrixop1[i1, i2, opn],
  0]
];

matrixopq1d[i1_, i2_, ch_] := Module[{op, opspinz},
  If[ch == 0,
  op = ch2op[ch];
  opn = number[ op[], DO ]-1/2;
  matrixop1[i1, i2, opn],
  0]
];


doall[] := Module[{},
  donew[];

  dogeneralloop[matrixnimix,      "OFFDIAG_MIX",       PREFIX <> "-offdiag-mix.dat"];
  dogeneralloop[matrixnimixup,    "OFFDIAG_MIX_UP",    PREFIX <> "-offdiag-mix-UP.dat"];
  dogeneralloop[matrixnimixdown,  "OFFDIAG_MIX_DOWN",  PREFIX <> "-offdiag-mix-DOWN.dat"];

  do2x2loop[matrixrunghop,     "RUNGHOP",      PREFIX <> "-runghop.dat"];
  do2x2loop[matrixrunghopup,   "RUNGHOP_UP",   PREFIX <> "-runghop-UP.dat"];
  do2x2loop[matrixrunghopdown, "RUNGHOP_DOWN", PREFIX <> "-runghop-DOWN.dat"];

  dooffdiag[];
  dodiag[];

  dorecalcf[];

  (* Arguments:
  - matrixopspinz: function which calls matrixop1[] to determine
    matrix elements
  - "SPINZ": function name in the C++ code
  - PREFIX <> "-spinz.dat": filename  
  - upper=True: only upper diagonal blocks will be considered. This is the
    approach required if the operator in matrixopspinz[] is the full
    self-conjugate operator, otherwise we would do double counting.
  *)
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", True];

  (* 2020 *)
  dogeneralloop[matrixopn1u, "N1U", PREFIX <> "-n1u.dat", True];
  dogeneralloop[matrixopn1d, "N1D", PREFIX <> "-n1d.dat", True];
  dogeneralloop[matrixopq1u, "Q1U", PREFIX <> "-q1u.dat", True];
  dogeneralloop[matrixopq1d, "Q1D", PREFIX <> "-q1d.dat", True];
  If[channels >= 2,
(*   dogeneralloop[matrixopn2u, "N2U", PREFIX <> "-n2u.dat", True]; *)
(*   dogeneralloop[matrixopn2d, "N2D", PREFIX <> "-n2d.dat", True]; *)
  ];

  (* The following are used for spin-polarized Wilson chains, therefore we need to
     include the zero terms (INCLUDEZEROONDIAG). *)
  INCLUDEZEROONDIAG = True;
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];

  dogeneralloop[matrixniup,    "OFFDIAG_UP",    PREFIX <> "-offdiag-UP.dat"];
  dogeneralloop[matrixnidown,  "OFFDIAG_DOWN",  PREFIX <> "-offdiag-DOWN.dat"];

  recalcop[{0}, diffsS];
  recalcop[{1/2}, diffsD];
  recalcop[{1}, diffsT];
];
