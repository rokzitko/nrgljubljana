(*
SYMTYPE=PP
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2017-2024
*)

PREFIX = "pp-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basis1 = qbasis[ {a[]} ];
basis2 = qbasis[ {b[]} ];

basis1 = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, basis1]; (* Charge parity *)
basis1 = mergebasis[basis1];

basis2 = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, basis2]; (* Charge parity *)
basis2 = mergebasis[basis2];

Print[basis1];
Print[basis2];

basis = basistensorproduct[basis1, basis2, Function[{qn1,qn2}, {qn1[[1]],qn2[[1]]}]];

Print[basis];

mult[_] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffpa_, diffpb_}] :=
  "Invar(" <> ToString[diffpa] <> ", " <> ToString[diffpb] <> ");";

InvarQN[{pa_, pb_}] :=  
  "Invar(" <> ToString[pa] <> ", " <> ToString[pb] <> ");";

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

newstates1[qn:{pa_, pb_}, state_] := newstates2[qn, {pa,pb}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{pa_,pb_}, d:{_,_}, state_] := 
  nc[coef[1], dirpdt[{{pa Pa, pb Pb}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{pa1_, pb1_}, diff1:{}, qn2:{pa2_, pb2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)

irrule[oppa_, oppb_] = {Pa -> Pa oppa, Pb -> Pb oppb};
ircg[oppa_, oppb_] = 1;

outmake[{pa_,pb_}] := "Invar(" <> ToString[pa /. Pa->pa1] <> ", " <> ToString[pb /. Pb->pb1] <> ")"; (* XXX ??? check in SYMTYPE=P *)

diffqn[qn1:{pa1_,pb1_}, qn2:{pa2_,pb2_}] := {Simplify[pa1 pa2, Pa^2 == 1 && Pb^2 == 1], Simplify[pb1 pb2, Pa^2 == 1 && Pb^2 == 1]}; (* Must be a list! *)

CG2recalcop[qn1:{pa1_,pb1_}, diff1:{},
            qn2:{pa2_,pb2_}, diff2:{}, 
            {_}] := 1;
 
QNrecalcop[delta:{deltapa_, deltapb_}] := (Pap = Pa deltapa; Pbp = Pb deltapb;);

RULE1recalcop[{pa_,pb_}] := {};
RULE2recalcop[{pa_,pb_}] := {Pa -> Pap, Pb -> Pbp};

CG1recalcop[{pa_, pb_}] := 1;

diffsS = { {{1, 1}, FNSINGLET<>".dat"} };
diffsDa = { {{-1, 1}, FNDOUBLET<>"a.dat"} }; (* changes parity: odd operator! *)
diffsDb = { {{1, -1}, FNDOUBLET<>"b.dat"} }; (* changes parity: odd operator! *)

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
dorecalcfLOOP["PP"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    allcode = {};
    op = ch2op[ch];
    (* fermion operators change parity!! *)
    If[ ch == 0, qn = {-1, 1}];
    If[ ch == 1, qn = {1, -1}];
    makerecalcf[op[CR, UP], qn, FNRECALCF <> "-" <> tos[op] <> "-CR-UP"];
    makerecalcf[op[CR, DO], qn, FNRECALCF <> "-" <> tos[op] <> "-CR-DO"];
    makerecalcf[op[AN, UP], qn, FNRECALCF <> "-" <> tos[op] <> "-AN-UP"];
    makerecalcf[op[AN, DO], qn, FNRECALCF <> "-" <> tos[op] <> "-AN-DO"];
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
  dorecalcf["PP"];

  checksignfn[{}] := False;
  recalcop[{1, 1}, diffsS];

  checksignfn[{}] := True;
  recalcop[{-1, 1}, diffsDa]; (* changes parity! *)
  recalcop[{1, -1}, diffsDb]; (* changes parity! *)

  INCLUDEZEROONDIAG = True;
  dogeneralloop[matrixopnumberup,   "DIAG_UP",   PREFIX <> "-diag-UP.dat",   True, INCLUDEZEROONDIAG];
  dogeneralloop[matrixopnumberdown, "DIAG_DOWN", PREFIX <> "-diag-DOWN.dat", True, INCLUDEZEROONDIAG];
];
