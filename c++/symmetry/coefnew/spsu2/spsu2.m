(* 
SYMTYPE=SPSU2
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

(* CHANGE LOG
22.11.2008 - support for anomalous hopping
17.10.2010 - qdiff, qtot
           - fixed simplification to work with mma7
15.12.2010 - Q1, Q2
13.10.2015 - three channels
10.5.2023 - fixed signs in 2nd channel
13.5.2023 - support for complex valued Wilson chain coefficients
*)

PREFIX = "spsu2-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = sbasis[ Take[allops, channels] ];
mult[{s_}] := 2s+1;
basisprep[];

(* Problem dependent utility functions *)
fixstate[s_, sz_, state_] := fixspin[s, sz, state];

Invar[{diffs_}] := "Invar(" <> ToString[2diffs] <> ");";
InvarQN[{s_}] := "Invar(" <> ToString[2s+1] <> ");";

simplifynew[expr_] := FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simpl2[expr_] := FullSimplify[ expr /. rule ,S+Sz \[Element] Integers];

simpl4[expr_] := FullSimplify[expr, 2S \[Element] Integers && 
                                    2Sz \[Element] Integers &&
                                    S+Sz \[Element] Integers &&
                                    S>3];

(* NOTE: taking Sz->S is dangerous when debugging, but it improves
   the calculation time considerably! *)
simpl5[expr_] := simpl4[expr /. Sz->S]; (* Take maximal Sz!*)
(* simpl5 = simpl4; *)

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

newstates1[qn:{s_}, state_] := 
  Do[newstates2[qn, {diffs}, state], {diffs, -s, s}];

(* NOTE: CG[] is just an alias for ClebschGordan[] *)

(* Produce the expression for the new state *)

getexpr1[qn:{s_}, d:{diffs_}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}],
    dirpdt[{{S+diffs}, {Sz-sz}, r}, fixstate[s, sz, state]]], {sz, -s, s}];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{s1_}, diff1:{sz1_}, qn2:{s2_}, diff2:{sz2_}, 
  opf:f[izn_, _, szn_]] := Module[{cg},
    cg = CG[{s2, sz2}, {1/2, spinofop[izn,szn]}, {s1, sz1}];
    If[ch =!= 0, Print["cg ", opf]];
    cg
  ];

matrixopisospin[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  If[ch == 2, Ix = isospinx[c[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinx1[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx1[a[],0] ];
  If[ch == 1, Ix = isospinx1[b[],0] ];
  If[ch == 2, Ix = isospinx1[c[],0] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinx2[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx2[a[],0] ];
  If[ch == 1, Ix = isospinx2[b[],0] ];
  If[ch == 2, Ix = isospinx2[c[],0] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

(* GLOBAL OPERATORS: phase important! *)

matrixopisospinGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  If[ch == 2, Ix = isospinx[c[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinzGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinz[a[]] ];
  If[ch == 1, Ix = isospinz[b[]] ];
  If[ch == 2, Ix = isospinz[c[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinplusGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinplus[a[]] ];
  If[ch == 1, Ix = isospinplus[b[]] ];
  If[ch == 2, Ix = isospinplus[c[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinminusGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinminus[a[]] ];
  If[ch == 1, Ix = isospinminus[b[]] ];
  If[ch == 2, Ix = isospinminus[c[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopchargeGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = number[a[]]-1 ]; (* Average charge subtracted! *)
  If[ch == 1, Ix = number[b[]]-1 ];
  If[ch == 2, Ix = number[c[]]-1 ];
  matrixop1[i1, i2, Ix]
];

(* For ch=2 only *)
matrixopqdiffGL[i1_, i2_, ch_] := Module[{op, opqdiff},
  op = ch2op[ch];
  (* ch numbering actually starts with ch=0, but we keep the following
     line as it is in qs/qs.m *)
  opqdiff = number[ op[] ]  If[ch == 1, +1, -1, MyError[]];
  matrixop1[i1, i2, opqdiff]
];

matrixopQ1GL[i1_, i2_, ch_] := Module[{op, opq1},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq1 = number[ op[] ];
  If[ch == 0, matrixop1[i1, i2, opq1], 0]
];

matrixopQ2GL[i1_, i2_, ch_] := Module[{op, opq2},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq2 = number[ op[] ];
  If[ch == 1, matrixop1[i1, i2, opq2], 0]
];

(* Be careful about the signs! In initial.m, we use anomaloushop[f[n], f[n+1]]
= f[CR, n, UP] f[CR, n+1, DO] - f[CR, n, DO] f[CR, n+1, UP] + h.c.
= f[CR, UP] op[CR, DO] - f[CR, DO] op[CR, UP + h.c.
*)

matrixanomalous[i1_, i2_, ch_] := Module[{expr=0, op},
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op[CR,DO]]];
  expr -= matrixel[i1, i2, nc[f[CR,ch,DO], op[CR,UP]]];
  (* NOTE: no Hermitian conjugate part here! *)
  If[expr =!= 0,
    MyPrint["#### ", expr];
  ];
  simpl4 @ expr
];

(* Recall: AN=1, CR=0 *)

(* recalcf[] stuff *)
irrule[opsz_] = {S -> S + opsz, Sz -> Sz + opsz};
ircg[opsz_] = CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

outmake[{b_}] := "Invar(" <> toc[2(b-S)+ss1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{s1_}, diff1:{sz1_}, 
            qn2:{s2_}, diff2:{sz2_}, 
            {m_}] :=
  CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];
 
(* Unprimed is bra, primed is ket! *)
QNrecalcop[delta:{deltas_}] := (Sp = S+deltas);

(* OLD: THIS IS NOT GENERAL ENOUGH! WE NEED CONSIDER {m,m} AS WELL AS {m,-m}
QUANTUM NUMBER COMBINATIONS!! *)
RULE1recalcop[{m_}] := {Sz -> S}; (* bra! *)
RULE2recalcop[{m_}] := {S -> Sp, Sz -> S-m}; (* ket *)
CG1recalcop[{m_}] := simplifynew @ CG[{Sp, S-m}, {m, m}, {S, S}];
(* NOTE: usually, CG1 should always be one, since we are combining
|Sp,S-m> with |m,m> to compare with <S,S|. *)

diffsD = {
  {{1/2}, FNDOUBLET<>"p.dat"},
  {{-1/2}, FNDOUBLET<>"m.dat"}
};

diffsS = { {{0}, FNSINGLET<>".dat"} };

diffsT = { 
  {{0}, FNTRIPLET<>"s.dat"}, 
  {{1}, FNTRIPLET<>"p.dat"},
  {{-1}, FNTRIPLET<>"m.dat"}
};

doall[] := Module[{},
  donew[];

  (* NEW, 22.10.2009 *)
  (* upper=False! Used only for global operator recalculations! *)
  dogeneralloop[matrixopisospinGL,      "ISOSPINX", PREFIX <> "-Ixtot.dat", False];
  dogeneralloop[matrixopisospinzGL,     "ISOSPINZ", PREFIX <> "-Iztot.dat", False];
  dogeneralloop[matrixopisospinplusGL,  "ISOSPINP", PREFIX <> "-Iptot.dat", False];
  dogeneralloop[matrixopisospinminusGL, "ISOSPINM", PREFIX <> "-Imtot.dat", False];
  dogeneralloop[matrixopchargeGL,       "CHARGE",   PREFIX <> "-Qtot.dat",  False];

  If[channels == 2,
    dogeneralloop[matrixopqdiffGL, "QDIFF", PREFIX <> "-qdiff.dat", False];
    dogeneralloop[matrixopQ1GL, "Q1", PREFIX <> "-q1.dat", False];
    dogeneralloop[matrixopQ2GL, "Q2", PREFIX <> "-q2.dat", False];
  ];

  dooffdiag[];
  dorecalcf[];
  recalcop[{1/2}, diffsD];
  dodiag[];

  dogeneralloop[matrixopisospin,  "ISOSPINX",  PREFIX <> "-isospinx.dat", True];
(*  dogeneralloop[matrixopisospin,  "ISOSPINX",  PREFIX <> "-isospinx-v2.dat", False]; *)
(*  dogeneralloop[matrixopisospinx1,"ISOSPINX1", PREFIX <> "-isospinx1.dat", False]; *)
(*  dogeneralloop[matrixopisospinx2,"ISOSPINX2", PREFIX <> "-isospinx2.dat", False]; *)

  dogeneralloop[matrixanomalous, "ANOMALOUS", PREFIX <> "-anomalous.dat", False]; (* note False! *)

  recalcop[{0}, diffsS];
  recalcop[{1}, diffsT]; (* TO DO: LOOP ??? *)
];
