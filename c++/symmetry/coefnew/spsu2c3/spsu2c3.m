(*
SYMTYPE=SPSU2C3
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2015-2023
$Id$
*)

(* CHANGE LOG *)
(* 2.11.2015 - first version *)
(* 5.11.2015 - modified in line with QSC3-standard *)

symtype = "SPSU2C3";
PREFIX = "spsu2c3";
fnprep[];

(* Make basis *)
myops = Take[allops, channels];
MyPrint["myops=", myops];

basis = Get["one_site_basis_spsu2c3_new_f.data"];
basis = basis //. {
 f[CR, 0, s_] :> a[CR, s],
 f[CR, 1, s_] :> b[CR, s],
 f[CR, 2, s_] :> c[CR, s]
 };

MyPrint["basis=", basis];

mult[{s_, p_}] := 2s+1;
basisprep[];

(* Problem dependent utility functions *)
fixstate[s_, sz_, state_] := fixspin[s, sz, state];

Invar[{diffs_, diffp_}] := "Invar" <> tos2[2diffs, diffp] <> ";"; 
InvarQN[{s_, p_}] := "Invar" <> tos2[2s+1, p] <> ";";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simpl2[expr_] := Simplify[expr /. rule, S+Sz \[Element] Integers && 2S \[Element] Integers];

simpl4[expr_] := Simplify[expr, 2S \[Element] Integers && S>2];
simpl5[expr_] := (simpl4[expr //. Sz->S] /. r:Root[__] :> N[r]) (* Take maximal Sz!*)

skpdtsimpl[expr_] := Simplify[expr, P \[Element] Integers];

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str, "S"->"S(ss)"] 
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
  Do[newstates2[qn, {diffs, Mod[-p, 3]}, state], {diffs, -s, s}];

(* Produce the expression for the new state *)

getexpr1[qn:{s_, p_}, d:{diffs_, _}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}],
    dirpdt[{{S+diffs, Mod[P-p,3]}, {Sz-sz}, r}, fixstate[s, sz, state]]], 
    {sz, -s, s}];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{s1_, p1_}, diff1:{sz1_}, 
           qn2:{s2_, p2_}, diff2:{sz2_}, 
  f[izn_, _, szn_]] := CG[{s2, sz2}, {1/2, spinofop[izn,szn]}, {s1, sz1}];

(* recalcf[] stuff; two arguments. called from makerecalcf[]. *)
irrule[opsz_, opp_] = {S -> S + opsz, Sz -> Sz + opsz, P -> Mod[P+opp,3]};
ircg[opsz_, opp_] = CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

outmake[{b_, c_}] :=
  "Invar" <> toc2[2(b-S)+ss1, c /. P->p1];

diffqn[qn1:{s1_, p1_}, qn2:{s2_, p2_}] := Module[{diff},
  diff = {s2-s1, Mod[p2-p1,3]}
];

CG2recalcop[qn1:{s1_, p1_}, diff1:{sz1_},
            qn2:{s2_, p2_}, diff2:{sz2_},
            {m_}] :=
  CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];
 
QNrecalcop[delta:{deltas_, deltap_}] := 
  (Sp = S+deltas; Pp = Mod[P+deltap,3]);

RULE1recalcop[{m_}] := {Sz -> S};
RULE2recalcop[{m_}] := {S -> Sp, Sz -> S-m, P -> Pp};

CG1recalcop[{m_}] := CG[{Sp, S-m}, {m, m}, {S, S}];

diffsD = {
  {{1/2, 0}, FNDOUBLET<>"p.dat"},
  {{-1/2, 0}, FNDOUBLET<>"m.dat"}
};

diffsS = { {{0, 0}, FNSINGLET<>".dat"} };

diffsT = {
  {{0, 0}, FNTRIPLET<>"s.dat"}, 
  {{1, 0}, FNTRIPLET<>"p.dat"},
  {{-1, 0}, FNTRIPLET<>"m.dat"}
};

dorecalcfLOOP["SPSU2C3", p_] := Module[{str, ch, op},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    MyPrint[1, "dorecalcfLOOP ch=", ch, " op=", op];
    makerecalcf[op[CR, UP], {1/2, p},  FNSPINUP   <> tos[p] <> "-" <> tos[op] ];
    makerecalcf[op[CR, DO], {-1/2, p}, FNSPINDOWN <> tos[p] <> "-" <> tos[op] ];
  ];
];

dorecalcfLOOP["SPSU2C3"] := Module[{},
 dorecalcfLOOP["SPSU2C3", 0];
(* dorecalcfLOOP["SPSU2C3", 1];
 dorecalcfLOOP["SPSU2C3", 2]; *)
];

SIMPLIFYNEW = True;
simplifynew[expr_] := Simplify[expr, 2S \[Element] Integers && S>=0];

matrixelsimpl[0] = 0;
matrixelsimpl[0.] = 0;
matrixelsimpl[0.+0.I] = 0;
matrixelsimpl[x_] := Chop[N[Simplify[x]]];

(* matrixop1 calls simpl2[] internally *)
matrixopdiagSUMCH[i_, i_] := matrixelsimpl @ simpl4[ matrixop1[i, i, number[a] + number[b] + number[c] ] ];

MATRIXNICHECKSYM[i1_, i2_] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1/2;

(* matrixel calls simpl5 internally. *)

matrixniNEWCR[i1_, i2_, ch_, spin_] := Module[{expr=0, o},
  If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
  o = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR, ch, spin], o[AN, spin]] ];
  expr = simpl4 @ expr;
  MyPrint["matrixNEWCR === ", expr];
  expr
  ];

(* h.c. of the above *)
matrixniNEWAN[i1_, i2_, ch_, spin_] := Module[{expr=0, o},
  If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
  o = ch2op[ch];
  expr += matrixel[i1, i2, nc[o[CR, spin], f[AN, ch, spin]] ];
  expr = simpl4 @ expr;
  MyPrint["matrixNEWAN === ", expr];
  expr
];

(* fixed *)
matrixopisospinx[i1_, i2_, op_] := Module[{Ix},
  Ix = simpl5 @ isospinx[op[]]; (* !! *)
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

(* sum over tz and sz: that's the final coefficient *)
matrixisospinxSUMCH[i1_, i2_] := Module[{expr},
  expr = matrixopisospinx[i1,i2,a] + matrixopisospinx[i1,i2,b] + matrixopisospinx[i1,i2,c];
  expr = simpl5[expr];
  MyPrint["matrixisospinxSUMCH === ", expr];
  expr
];

doisospinxSUMCH[] := Module[{},
  MyPrint["doisospinxSUMCH[]"];
  dogeneralloopSUMCH[matrixisospinxSUMCH,
  "ISOSPINX",
  PREFIX <> "-isospinx.dat", True];
];

doall[] := Module[{},
  donew[];

  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;

  dooffdiag[];
(*  dodiagSUMCH[];
  doisospinxSUMCH[]; *)

  dorecalcf["SPSU2C3"];

  (* doanomalousSUMCH[]; NOT YET *)

  recalcop[{1/2}, diffsD];
  recalcop[{0}, diffsS];
  recalcop[{1}, diffsT];
];
