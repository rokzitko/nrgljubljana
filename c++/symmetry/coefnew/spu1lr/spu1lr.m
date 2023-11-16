(*
SYMTYPE=SPU1LR
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2010-2023
*)

(* CHANGE LOG
11.5.2010 - first version
11.5.2023 - fixed signs for isospinx (2nd channel)
*)

PREFIX = "spu1lr-" <> ToString[channels] <> "ch";
fnprep[];

szbasisvc[l_List] := Module[{bvc},
  bvc = qszbasisvc[l];
  bvc = Map[{ {#[[1,2]]}, #[[2]] }&, bvc];
  bvc = mergebasis[bvc];
  bvc
];

szbasis[l_List] := bzvc2bzop @ szbasisvc[l];

(* Make basis *)
myops = Take[allops, channels];
makebasis[myops];
Print["myops=", myops, " BASIS=", BASIS];

basis = szbasis[myops];
basis = transformtoLRvc[basis, myops, vacuum[]];
basis = bzvc2bzop @ basis;

Print["basis=", basis];
mult[{__}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffsz_, diffp_}] := 
 "Invar(" <> ToString[2diffsz] <> ", " <> ToString[diffp] <> ");";

InvarQN[{sz_, p_}] := 
 "Invar(" <> ToString[2sz] <> ", " <> ToString[p] <> ");";

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

newstates1[qn:{sz_, p_}, state_] := 
  newstates2[qn, {-sz, p}, state];

(* NOTE: CG[] is just an alias for ClebschGordan[] *)

(* Produce the expression for the new state *)

getexpr1[qn:{sz_, p_}, d:{diffsz_, diffp_}, state_] := 
  nc[coef[1], dirpdt[{{Sz+diffsz, P p}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{sz1_, p1_}, diff1:{}, qn2:{sz2_, p2_}, diff2:{}, 
  opf:f[izn_, _, szn_]] := 1;

matrixopisospin[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

(* GLOBAL operators *)

matrixopisospinGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinzGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinz[a[]] ];
  If[ch == 1, Ix = isospinz[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinplusGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinplus[a[]] ];
  If[ch == 1, Ix = isospinplus[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinminusGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinminus[a[]] ];
  If[ch == 1, Ix = isospinminus[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopchargeGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = number[a[]]-1 ]; (* Average charge subtracted! *)
  If[ch == 1, Ix = number[b[]]-1 ];
  matrixop1[i1, i2, Ix]
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
irrule[opsz_, opp_] = {Sz -> Sz + opsz, P -> P opp};
ircg[opsz_,opp_] = 1;

outmake[{b_,c_}] := "Invar" <> toc2[2(b-Sz)+ssz1, c /. P -> p1];

diffqn[qn1:{sz1_, p1_}, qn2:{sz2_,p2_}] := Module[{diff},
 diff = {sz2-sz1, p1 p2};
 Simplify[diff, P^2 == 1]
];

CG2recalcop[qn1:{sz1_, p1_}, diff1:{}, 
            qn2:{sz2_, p2_}, diff2:{}, 
            {m_}] := 1;

(* Unprimed is bra, primed is ket! *)
QNrecalcop[delta:{deltasz_, deltap_}] := (Szp = Sz+deltasz; Pp = P deltap);

RULE1recalcop[{m_}] := {};
RULE2recalcop[{m_}] := { Sz -> Szp };

CG1recalcop[{m_}] := 1;

diffsS = { {{0, 1}, FNSINGLET<>".dat"} };

diffsD = {
  {{1/2, 1},  FNDOUBLET<>"p.dat"},
  {{-1/2, 1}, FNDOUBLET<>"m.dat"}
};

diffsT = { 
  {{0, 1}, FNTRIPLET<>"s.dat"}, 
  {{1, 1}, FNTRIPLET<>"p.dat"},
  {{-1, 1}, FNTRIPLET<>"m.dat"}
};

dorecalcfLOOP["SPU1LR", p_] := Module[{str, ch, op},
  str = If[p == 1, "", "diff"];
  For[ch = 0, ch <= channels-1, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2, p},
      FNSPINUP <> str <> tos[op]];
    makerecalcf[op[CR, DO], {-1/2, p},
      FNSPINDOWN <> str <> tos[op] ];
  ];
];

dorecalcfLOOP["SPU1LR"] := Module[{},  
  dorecalcfLOOP["SPU1LR", 1];
  If[channels == 2, 
    dorecalcfLOOP["SPU1LR", -1]
  ];
];

doall[] := Module[{},
  donew[];
  dooffdiag[];
  dorecalcf["SPU1LR"];

  (* NEW, 23.10.2009 *)
  (* upper=False! Used only for global operator recalculations! *)
  dogeneralloop[matrixopisospinGL,      "ISOSPINX", PREFIX <> "-Ixtot.dat", False];

  dogeneralloop[matrixopisospinzGL,     "ISOSPINZ", PREFIX <> "-Iztot.dat", False];
  dogeneralloop[matrixopisospinplusGL,  "ISOSPINP", PREFIX <> "-Iptot.dat", False];
  dogeneralloop[matrixopisospinminusGL, "ISOSPINM", PREFIX <> "-Imtot.dat", False];
  dogeneralloop[matrixopchargeGL,  "CHARGE",   PREFIX <> "-Qtot.dat",  False];

  (* NEW, 22.11.2008 *)
  dogeneralloop[matrixanomalous, "ANOMALOUS", PREFIX <> "-anomalous.dat"];

  recalcop[{1/2}, diffsD];
  AppendTo[defs, defstr["LENGTH_D_" <> ToString[channels] <> "CH",LEN]];

  dodiag[];
  dogeneralloop[matrixopisospin,"ISOSPINX", PREFIX <> "-isospinx.dat", True];

  recalcop[{0}, diffsS];
  recalcop[{1/2}, diffsDself];
  recalcop[{1}, diffsT]; (* TO DO: LOOP ??? *)
];
