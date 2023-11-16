(*
SYMTYPE=P
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2017-2023
*)

PREFIX = "p-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis = qbasis[ Take[allops, channels] ];

Print[basis];

basis = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, basis]; (* Charge parity *)
basis = mergebasis[basis];

Print[basis];

mult[_] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffp_}] :=
  "Invar(" <> ToString[diffp] <> ");";

InvarQN[{p_}] :=
  "Invar(" <> ToString[p] <> ");";

simpl2[expr_] := Simplify[expr];

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

newstates1[qn:{p_}, state_] := newstates2[qn, {p}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{p_}, d:{_}, state_] := 
  nc[coef[1], dirpdt[{{p P}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{p1_}, diff1:{}, qn2:{p2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)

irrule[opp_] = {P -> P opp};
ircg[opp_] = 1;

outmake[{p_}] := Module[{repl,str},
  repl = p /. P->p1;
  str = "Invar(" <> ToString[repl] <> ")";
  Print["outmake", {p}, {repl}, str];
  str
];

diffqn[qn1:{p1_}, qn2:{p2_}] := {Simplify[p1 p2, P^2 == 1]}; (* Must be a list! *)

CG2recalcop[qn1:{p1_}, diff1:{},
            qn2:{p2_}, diff2:{}, 
            {_}] := 1;
 
QNrecalcop[delta:{deltap_}] := Pp = P deltap;

RULE1recalcop[{p_}] := {};
RULE2recalcop[{p_}] := {P -> Pp};

CG1recalcop[{p_}] := 1;

diffsS = { {{1}, FNSINGLET<>".dat"} };
diffsD = { {{-1}, FNDOUBLET<>".dat"} }; (* changes parity: odd operator! *)

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
  If[ch == 1, Ix = isospinplus[b[]] ]; (* sign? *)
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

matrixopisospinminusGL[i1_, i2_, ch_] := Module[{Ix},
  If[ch == 0, Ix = isospinminus[a[]] ];
  If[ch == 1, Ix = isospinminus[b[]] ]; (* sign? *)
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

(* Based on SYMTYPE=NONE *)
dorecalcfLOOP["P"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    allcode = {};
    op = ch2op[ch];
    (* fermion operators change parity!! *)
    makerecalcf[op[CR, UP], {-1}, FNRECALCF <> "-" <> tos[op] <> "-CR-UP"];
    makerecalcf[op[CR, DO], {-1}, FNRECALCF <> "-" <> tos[op] <> "-CR-DO"];
    makerecalcf[op[AN, UP], {-1}, FNRECALCF <> "-" <> tos[op] <> "-AN-UP"];
    makerecalcf[op[AN, DO], {-1}, FNRECALCF <> "-" <> tos[op] <> "-AN-DO"];
  ];
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
  dorecalcf["P"];

  checksignfn[{}] := False;
  recalcop[{1}, diffsS];

  checksignfn[{}] := True;
  recalcop[{-1}, diffsD]; (* changes parity! *)

  INCLUDEZEROONDIAG = True;
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];
];
