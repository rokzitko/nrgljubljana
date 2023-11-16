(*
SYMTYPE=SPSU2LR
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2015-2023
*)

(* CHANGE LOG
24.9.2015 - first version
*)

PREFIX = "spsu2lr-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
myops = Take[allops, channels];
makebasis[myops];
Print["myops=", myops, " BASIS=", BASIS];

basis = sbasis[myops];

Print["basis (no LR) = ", basis];

basis = transformtoLRvc[basis, myops, vacuum[] ];
basis = bzvc2bzop @ basis;

Print["basis (after LR) = ", basis];

mult[{s_, _}] := 2s+1;
basisprep[];

(* Problem dependent utility functions *)
fixstate[s_, sz_, state_] := fixspin[s, sz, state];

Invar[{diffs_, diffp_}] := "Invar(" <> ToString[2diffs] <> ", " <> ToString[diffp] <> ");";
InvarQN[{s_, p_}] := "Invar(" <> ToString[2s+1] <> ", " <> ToString[p] <> ");";

simplifynew[expr_] := FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simpl2[expr_] := FullSimplify[expr /. rule,S+Sz \[Element] Integers];

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

newstates1[qn:{s_, p_}, state_] := 
  Do[newstates2[qn, {diffs, p}, state], {diffs, -s, s}];

(* NOTE: CG[] is just an alias for ClebschGordan[] *)

(* Produce the expression for the new state *)

getexpr1[qn:{s_, p_}, d:{diffs_, diffp_}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}],
    dirpdt[{{S+diffs, P p}, {Sz-sz}, r}, fixstate[s, sz, state]]], {sz, -s, s}];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{s1_, p1_}, diff1:{sz1_}, qn2:{s2_, p2_}, diff2:{sz2_}, 
  opf:f[izn_, _, szn_]] := Module[{cg},
    cg = CG[{s2, sz2}, {1/2, spinofop[izn,szn]}, {s1, sz1}];
    If[ch =!= 0, Print["cg ", opf]];
    cg
  ];

matrixopisospin[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

(* Be careful about the signs! In initial.m, we use anomaloushop[f[n], f[n+1]]
= f[CR, n, UP] f[CR, n+1, DO] - f[CR, n, DO] f[CR, n+1, UP] + h.c.
= f[CR, UP] op[CR, DO] - f[CR, DO] op[CR, UP + h.c.
*)

matrixanomalous[i1_, i2_, ch_] := Module[{expr=0, op},
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op[CR,DO]]];
  expr -= matrixel[i1, i2, nc[f[CR,ch,DO], op[CR,UP]]];
  If[expr =!= 0,
    MyPrint["#### ", expr];
  ];
  simpl4 @ expr
];

(* Recall: AN=1, CR=0 *)

(* recalcf[] stuff *)
irrule[opsz_, opp_] = {S -> S + opsz, Sz -> Sz + opsz, P -> P opp};
ircg[opsz_, opp_] = CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

(* toc2 adds parenthesis *)
outmake[{b_,c_}] := "Invar" <> toc2[2(b-S)+ss1, c /. P -> p1];

diffqn[qn1:{s1_, p1_}, qn2:{s2_, p2_}] := Module[{diff},
  diff = {s2-s1, p2 p1}; (* CAREFUL HERE! *)
  diff = Simplify[diff, P^2 == 1];
  diff
];

CG2recalcop[qn1:{s1_, p1_}, diff1:{sz1_}, 
            qn2:{s2_, p2_}, diff2:{sz2_}, 
            {m_}] :=
  CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];
 
(* Unprimed is bra, primed is ket! *)
QNrecalcop[delta:{deltas_, deltap_}] := (Sp = S+deltas; Pp = P deltap);

(* OLD: THIS IS NOT GENERAL ENOUGH! WE NEED CONSIDER {m,m} AS WELL AS {m,-m}
QUANTUM NUMBER COMBINATIONS!! *)
RULE1recalcop[{m_}] := {Sz -> S}; (* bra! *)
RULE2recalcop[{m_}] := {S -> Sp, Sz -> S-m, P -> Pp}; (* ket *)
CG1recalcop[{m_}] := simplifynew @ CG[{Sp, S-m}, {m, m}, {S, S}];
(* NOTE: usually, CG1 should always be one, since we are combining
|Sp,S-m> with |m,m> to compare with <S,S|. *)

diffsD = {
  {{1/2, 1}, FNDOUBLET<>"p.dat"},
  {{-1/2, 1}, FNDOUBLET<>"m.dat"}
};

diffsDodd = {
  {{1/2, -1}, FNDOUBLET<>"diffp.dat"},
  {{-1/2, -1}, FNDOUBLET<>"diffm.dat"}
};

diffsS = { {{0, 1}, FNSINGLET<>".dat"} };

diffsT = {
  {{0, 1}, FNTRIPLET<>"s.dat"}, 
  {{1, 1}, FNTRIPLET<>"p.dat"},
  {{-1, 1}, FNTRIPLET<>"m.dat"}
};

dorecalcfLOOP["SPSU2LR", p_] := Module[{str, ch, op},
  str = If[p == 1, "", "diff"];
  For[ch = 0, ch <= channels-1, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2, p},
      FNSPINUP <> str <> tos[op]];
    makerecalcf[op[CR, DO], {-1/2, p},
      FNSPINDOWN <> str <> tos[op] ];
  ];
];

dorecalcfLOOP["SPSU2LR"] := Module[{},
  dorecalcfLOOP["SPSU2LR", 1];
  If[channels == 2,
    dorecalcfLOOP["SPSU2LR", -1]
  ];
];

doall[] := Module[{},
  donew[];
  dogeneralloop[matrixanomalous, "ANOMALOUS", PREFIX <> "-anomalous.dat"];
  dooffdiag[];
  dorecalcf["SPSU2LR"];
  dodiag[];
  dogeneralloop[matrixopisospin,"ISOSPINX", PREFIX <> "-isospinx.dat", True];

  recalcop[{0}, diffsS];
  recalcop[{1/2}, diffsD];
  recalcop[{1/2}, diffsDodd];
  recalcop[{1}, diffsT];
];
