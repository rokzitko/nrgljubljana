(*
SYMTYPE=SPU1
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

(* CHANGE LOG
4.12.2008 - first version
24.1.2011 - QDIFF, Q1, Q2, Q1UP, Q1DO, Q2UP, Q2DO
10.5.2023 - fix sign in isospix for 2nd channel
*)

PREFIX = "spu1-" <> ToString[channels] <> "ch";
fnprep[];

szbasisvc[l_List] := Module[{bvc},
  bvc = qszbasisvc[l];
  bvc = Map[{ {#[[1,2]]}, #[[2]] }&, bvc];
  bvc = mergebasis[bvc];
  bvc
];

szbasis[l_List] := bzvc2bzop @ szbasisvc[l];

(* Make basis *)
basis = szbasis[ Take[allops, channels] ];
mult[{__}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffsz_}] := "Invar(" <> ToString[2diffsz] <> ");";

InvarQN[{sz_}] := "Invar(" <> ToString[2sz] <> ");";

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

newstates1[qn:{sz_}, state_] := 
  newstates2[qn, {-sz}, state];

(* NOTE: CG[] is just an alias for ClebschGordan[] *)

(* Produce the expression for the new state *)

getexpr1[qn:{sz_}, d:{diffsz_}, state_] := 
  nc[coef[1], dirpdt[{{Sz+diffsz}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{sz1_}, diff1:{}, qn2:{sz2_}, diff2:{}, 
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

matrixopQ1UPGL[i1_, i2_, ch_] := Module[{op, opq1},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq1 = number[ op[], UP ];
  If[ch == 0, matrixop1[i1, i2, opq1], 0]
];

matrixopQ2UPGL[i1_, i2_, ch_] := Module[{op, opq2},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq2 = number[ op[], UP ];
  If[ch == 1, matrixop1[i1, i2, opq2], 0]
];

matrixopQ1DOGL[i1_, i2_, ch_] := Module[{op, opq1},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq1 = number[ op[], DO ];
  If[ch == 0, matrixop1[i1, i2, opq1], 0]
];

matrixopQ2DOGL[i1_, i2_, ch_] := Module[{op, opq2},
  op = ch2op[ch];
  (* ch numbering starts with ch=0 *)
  opq2 = number[ op[], DO ];
  If[ch == 1, matrixop1[i1, i2, opq2], 0]
];



(* Recall: AN=1, CR=0 *)

(* recalcf[] stuff *)
irrule[opsz_] = {Sz -> Sz + opsz};
ircg[opsz_] = 1;

outmake[{b_}] := "Invar(" <> toc[2(b-Sz)+ssz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{sz1_}, diff1:{}, 
            qn2:{sz2_}, diff2:{}, 
            {m_}] := 1;
 
(* Unprimed is bra, primed is ket! *)
QNrecalcop[delta:{deltasz_}] := (Szp = Sz+deltasz);

RULE1recalcop[{m_}] := {};
RULE2recalcop[{m_}] := { Sz -> Szp };

CG1recalcop[{m_}] := 1;

diffsS = { {{0}, FNSINGLET<>".dat"} };

diffsD = {
  {{1/2},  FNDOUBLET<>"p.dat"},
  {{-1/2}, FNDOUBLET<>"m.dat"}
};

diffsT = { 
  {{0}, FNTRIPLET<>"s.dat"}, 
  {{1}, FNTRIPLET<>"p.dat"},
  {{-1}, FNTRIPLET<>"m.dat"}
  };

(* NEW, 28.5.2012 *)
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

doall[] := Module[{},
  donew[];

  (* NEW, 23.10.2009 *)
  (* upper=False! Used only for global operator recalculations! *)
  dogeneralloop[matrixopisospinGL,      "ISOSPINX", PREFIX <> "-Ixtot.dat", False];
  dogeneralloop[matrixopisospinzGL,     "ISOSPINZ", PREFIX <> "-Iztot.dat", False];
  dogeneralloop[matrixopisospinplusGL,  "ISOSPINP", PREFIX <> "-Iptot.dat", False];
  dogeneralloop[matrixopisospinminusGL, "ISOSPINM", PREFIX <> "-Imtot.dat", False];
  dogeneralloop[matrixopchargeGL,  "CHARGE",   PREFIX <> "-Qtot.dat",  False];

  (* NEW, 22.11.2008 *)
  dogeneralloop[matrixanomalous, "ANOMALOUS", PREFIX <> "-anomalous.dat"];

  If[channels == 2,
    dogeneralloop[matrixopqdiffGL, "QDIFF", PREFIX <> "-qdiff.dat", False];
    dogeneralloop[matrixopQ1GL, "Q1", PREFIX <> "-q1.dat", False];
    dogeneralloop[matrixopQ2GL, "Q2", PREFIX <> "-q2.dat", False];
    dogeneralloop[matrixopQ1UPGL, "Q1UP", PREFIX <> "-q1up.dat", False];
    dogeneralloop[matrixopQ2UPGL, "Q2UP", PREFIX <> "-q2up.dat", False];
    dogeneralloop[matrixopQ1DOGL, "Q1DO", PREFIX <> "-q1do.dat", False];
    dogeneralloop[matrixopQ2DOGL, "Q2DO", PREFIX <> "-q2do.dat", False];
  ];

  dooffdiag[];
  dorecalcf[];
  dodiag[];

  recalcop[{0}, diffsS];
  recalcop[{1/2}, diffsD];
  recalcop[{1}, diffsT]; (* TO DO: LOOP ??? *)

  dogeneralloop[matrixopisospin,"ISOSPINX", PREFIX <> "-isospinx.dat", True];

  (* The following are used for spin-polarized Wilson chains, therefore we need to                                             
    include the zero terms (INCLUDEZEROONDIAG). *)                                                                            
  INCLUDEZEROONDIAG = True;                                                                                                    
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];                         
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];                         

  dogeneralloop[matrixniup,    "OFFDIAG_UP",    PREFIX <> "-offdiag-UP.dat"];                                                  
  dogeneralloop[matrixnidown,  "OFFDIAG_DOWN",  PREFIX <> "-offdiag-DOWN.dat"];   
];
