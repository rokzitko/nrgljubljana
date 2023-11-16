(*
SYMTYPE=U1
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

(* 
CHANGE LOG
----------
10.9.2012 - out-of-diagonal elements
*)

PREFIX = "u1-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = qbasis[ Take[allops, channels] ];
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

  CGmatrixelV2[qn1:{q1_}, diff1:{}, qn2:{q2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)
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

diffsD = {
  {{-1}, FNDOUBLET<>".dat"}
};

diffsS = { {{0}, FNSINGLET<>".dat"} };

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

matrixopspiny[i1_, i2_, ch_] := Module[{op, opspiny},
  op = ch2op[ch];
  opspiny = spiny[ op[] ];
  matrixop1[i1, i2, opspiny]
];

matrixopspinx[i1_, i2_, ch_] := Module[{op, opspinx},
  op = ch2op[ch];
  opspinx = spinx[ op[] ];
  matrixop1[i1, i2, opspinx]
];

matrixopnumberup[i1_, i2_, ch_] := Module[{op},
  op = ch2op[ch];
  opnumberup = number[ op[], UP];
  matrixop1[i1, i2, opnumberup]
];

matrixopnumberdown[i1_, i2_, ch_] := Module[{op},
  op = ch2op[ch];
  opnumberdown = number[ op[], DO];
  matrixop1[i1, i2, opnumberdown]
];

(* Added 10.9.2012 *)
matrixopnumberupdo[i1_, i2_, ch_] := Module[{op},
  op = ch2op[ch];
  opx = op[CR, UP] ~ nc ~ op[AN,DO];
  matrixop1[i1, i2, opx]
];

matrixopnumberdoup[i1_, i2_, ch_] := Module[{op},
  op = ch2op[ch];
  opx = op[CR, DO] ~ nc ~ op[AN,UP];
  matrixop1[i1, i2, opx]
  ];

matrixhopdoup[i1_, i2_, ch_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];  
  op = ch2op[ch];
  expr += matrixelV2[i1, i2, nc[f[CR,ch,DO], op[AN,UP]]];
  If[expr =!= 0,
  MyPrint["#### ", expr];
  ];
  simpl4 @ expr
];

matrixhopupdo[i1_, i2_, ch_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];  
  op = ch2op[ch];
  expr += matrixelV2[i1, i2, nc[f[CR,ch,UP], op[AN,DO]]];
  If[expr =!= 0,
  MyPrint["#### ", expr];
  ];
  simpl4 @ expr
];

doall[] := Module[{},
  donew[];

  dogeneralloop[matrixniNEW[#1,#2,#3,DO]&,  "OFFDIAG_DO",
    PREFIX <> "-offdiag-DO.dat", True]; (* True necessary to avoid double counting *)
  dogeneralloop[matrixniNEW[#1,#2,#3,UP]&,  "OFFDIAG_UP",
    PREFIX <> "-offdiag-UP.dat", True];

  (* Added 10.9.2012 *)
  dogeneralloop[matrixhopdoup, "OFFDIAG_DOUP", 
    PREFIX <> "-offdiag-DOUP.dat"];
  dogeneralloop[matrixhopupdo, "OFFDIAG_UPDO", 
    PREFIX <> "-offdiag-UPDO.dat"];

  (* For gluing together matrix parts (off-diagonal matrix blocks between
     different invariant subspaces). *)
  dooffdiagNEW[]; (* ***!!*** *)

  (* u1-?ch-diag.dat is not actually used, only the u1-?ch-diag-UP/DOWN/DOUP.dat part
     are needed. *)
  dodiag[];

  (* ATTENTION: upper=False here! These are not used to construct the Hamiltonian... *)
  dogeneralloop[matrixopspinx, "SPINX", PREFIX <> "-spinx.dat", False];
  dogeneralloop[matrixopspiny, "SPINY", PREFIX <> "-spiny.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];

  dorecalcf["U1"];

  (* Added 14.10.2009 *)
  INCLUDEZEROONDIAG = True;
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];

  (* Added 10.9.2012 *)
  (* upper=True because these are used to construct the Hamiltonian. DOUP does contribute in the upper part, but UPDO doesn't! *)
  (* We handle this by analogy with SYMTYPE=SPSU2: we just use one of them! The Hamiltonian is Hermitian anyway!! *)
  (* THIS ONE NOT USED: dogeneralloop[matrixopnumberupdo, "DIAG_UPDO", PREFIX <> "-diag-UPDO.dat", True, INCLUDEZEROONDIAG]; *)
  dogeneralloop[matrixopnumberdoup, "DIAG_DOUP", PREFIX <> "-diag-DOUP.dat", True, INCLUDEZEROONDIAG];

  checksignfn[{}] := True;
  recalcop[{}, diffsD];
  checksignfn[{}] := False;
  recalcop[{}, diffsS];
];
