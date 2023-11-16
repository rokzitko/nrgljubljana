(*
SYMTYPE=NONE
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "none-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = nonebasis[ Take[allops, channels] ];
mult[{}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{}] :=
  "Invar(0);"

InvarQN[{}] := 
  "Invar(0);";

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

newstates1[qn:{}, state_] := newstates2[qn, {}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{}, d:{}, state_] := 
  nc[coef[1], dirpdt[{{}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{}, diff1:{}, qn2:{}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)

irrule[___] = {};
ircg[___] = 1;

outmake[{}] :=
  "Invar(0)";

diffqn[qn1:{___}, qn2:{___}] := qn2-qn1;

CG2recalcop[qn1:{}, diff1:{},
            qn2:{}, diff2:{}, 
            {}] := 1;

QNrecalcop[delta:{}] := 1;

RULE1recalcop[{}] := {};
RULE2recalcop[{}] := {};

CG1recalcop[{}] := 1;

diffsS = { {{}, FNSINGLET<>".dat"} };
diffsD = { {{}, FNDOUBLET<>".dat"} };

matrixopisospinGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinx[a[]] ];
  If[ch == 1, Ix = isospinx[b[]] ];
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinyGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospiny[a[]] ];
  If[ch == 1, Ix = isospiny[b[]] ];
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

(* NEW, from SPU1, 24.6.2016 *)
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

  dogeneralloop[matrixopisospinGL,      "ISOSPINX", PREFIX <> "-Ixtot.dat", False];
  dogeneralloop[matrixopisospinyGL,     "ISOSPINY", PREFIX <> "-Iytot.dat", False];
  dogeneralloop[matrixopisospinzGL,     "ISOSPINZ", PREFIX <> "-Iztot.dat", False];
  dogeneralloop[matrixopisospinplusGL,  "ISOSPINP", PREFIX <> "-Iptot.dat", False];
  dogeneralloop[matrixopisospinminusGL, "ISOSPINM", PREFIX <> "-Imtot.dat", False];
  dogeneralloop[matrixopchargeGL,       "CHARGE",   PREFIX <> "-Qtot.dat",  False];

  dogeneralloop[matrixopspinx, "SPINX", PREFIX <> "-spinx.dat", False];
  dogeneralloop[matrixopspiny, "SPINY", PREFIX <> "-spiny.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];

  (* NOT SUFFICIENT! : dooffdiagNEW[]; *)

  dogeneralloop[matrixniNEWCR[#1,#2,#3,DO]&,  "OFFDIAG_CR_DO",
    PREFIX <> "-offdiag-CR-DO.dat", True];
  dogeneralloop[matrixniNEWCR[#1,#2,#3,UP]&,  "OFFDIAG_CR_UP",
    PREFIX <> "-offdiag-CR-UP.dat", True];
  dogeneralloop[matrixniNEWAN[#1,#2,#3,DO]&,  "OFFDIAG_AN_DO",
    PREFIX <> "-offdiag-AN-DO.dat", True];
  dogeneralloop[matrixniNEWAN[#1,#2,#3,UP]&,  "OFFDIAG_AN_UP",
    PREFIX <> "-offdiag-AN-UP.dat", True];

  dodiag[];
  dorecalcf["NONE"];

  checksignfn[{}] := False;
  recalcop[{}, diffsS];
  checksignfn[{}] := True;
  recalcop[{}, diffsD];

  INCLUDEZEROONDIAG = True;
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];
];
