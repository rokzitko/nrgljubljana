If[!ValueQ[SYMTYPE],  SYMTYPE = "QS"];

(* Implemented symmetry types are:
 "QS", ==>                 U(1)_charge x SU(2)_spin
 "DBLQS", ==>              U(1)_charge1 x U(1)_charge2 x SU(2)_spin [TO DO!!]
 "QST", ==>                U(1)_charge x SU(2)_spin x SO(3)_orbital [three orbitals]
 "QSTZ", ==>               U(1)_charge x SU(2)_spin x U(1)_orbital [three orbitals]
 "QSZTZ", ==>              U(1)_charge x U(1)_spin x U(1)_orbital [three orbitals]
 "QJ", ==>                 U(1)_charge x SU(2)_total_momentum [three orbitals]
 "ISO" and "ISO2", ==>     SU(2)_isospin x SU(2)_spin
 "QSZ", ==>                U(1)_charge x U(1)_spin
 "DBLQSZ", ==>             U(1)_charge1 x U(1)_charge2 x U(1)_spin
 "ISOSZ", ==>              SU(2)_isospin x U(1)_spin
 "ISOSZLR", ==>            SU(2)_isospin x U(1)_spin x Z_2
 "ISOLR" and "ISO2LR", ==> SU(2)_isospin x SU(2)_spin x Z_2
 "QSLR", ==>               U(1)_charge x SU(2)_spin x Z_2
 "QSC3", ==>               U(1)_charge x SU(2)_spin x Z_3 [3 channels]
 "QSZLR", ==>              U(1)_charge x U(1)_spin x Z_2
 "SU2", ==>                SU(2)_charge
 "DBLSU2", ==>             SU(2)_charge1 x SU(2)_charge2
 "DBLISOSZ", ==>           SU(2)_charge1 x SU(2)_charge2 x U(1)_spin
 "U1", ==>                 U(1)_charge
 "SPSU2", ==>              SU(2)_spin
 "SPU1", ==>               U(1)_spin
 "SPU1LR", ==>             U(1)_spin x Z_2
 "SPSU2LR", ==>            SU(2)_spin x Z_2
 "SPSU2C3", ==>            SU(2)_spin x Z_3 [3 channels]
 "SPSU2T", ==>             SU(2)_spin x SO(3)_orbital [three orbitals]
 "SL", ==>                 spinless fermions [U(1)_charge symmetry]
 "SL3", ==>                spinless fermions [threefold U(1)_charge]
 "P", ==>                  Z_2 fermion parity,
 "PP", ==>                 (Z_2)^2 fermion parity,
 "NONE", ==>               no symmetry
*)

If[SYMTYPE == "runtime",
  If[!paramexists["symtype"],
    MyError["Cannot determine the symmetry type."];
  ];
  SYMTYPE = param["symtype"];
];

knownsymtypes =
{"QS", "QST", "QSTZ", "QSZTZ", "QJ", "ISO", "ISO2", "QSZ", "ISOLR", "ISO2LR", "QSLR", "QSC3",
"QSZLR", "DBLQSZ", "DBLSU2", "DBLISOSZ",
"SU2", "U1", "SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2C3", "SPSU2T",
"SL", "SL3", "P", "PP", "NONE", "ISOSZ", "ISOSZLR"};
If[!(MemberQ[knownsymtypes, SYMTYPE]),
  MyError["Unknown SYMTYPE."];
];

(* isLR[] returns True if the problem has Z_2 parity symmetry. For symmetry
types in this list, we create a parity-adapted set of basis states as
step number 5 in the basis-state-generation part of the code. *)
lrsymtypes = {"QSLR", "QSZLR", "ISOLR", "ISO2LR", "ISOSZLR", "SPU1LR", "SPSU2LR"};
isLR[] := MemberQ[lrsymtypes, SYMTYPE];

qstypes = {"QS", "QSLR"};
qsztypes = {"QSZ", "QSZLR"};
isotypes = {"ISO", "ISO2", "ISOLR", "ISO2LR"};
isosztypes = { "ISOSZ", "ISOSZLR" };
sctypes = {"SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2T", "P", "PP", "NONE", "SPSU2C3"};
orbtypes = {"QST", "SPSU2T", "QSTZ", "QSZTZ"}; (* quantum number T *)
isQS[] := MemberQ[qstypes, SYMTYPE];
isQSZ[] := MemberQ[qsztypes, SYMTYPE];
isISO[] := MemberQ[isotypes, SYMTYPE];
isISOSZ[] := MemberQ[isosztypes, SYMTYPE];
isSC[] := MemberQ[sctypes, SYMTYPE];
isORB[] := MemberQ[orbtypes, SYMTYPE];

isSU2[] := (SYMTYPE === "SU2");
isDBLSU2[] := (SYMTYPE === "DBLSU2");
isDBLISOSZ[] := (SYMTYPE === "DBLISOSZ");
isDBLQSZ[] := (SYMTYPE == "DBLQSZ");
isU1[] := (SYMTYPE === "U1");
isSPSU2[] := (SYMTYPE === "SPSU2");
isSPU1[] := (SYMTYPE === "SPU1");
isSPU1LR[] := (SYMTYPE === "SPU1LR");
isSPSU2LR[] := (SYMTYPE === "SPSU2LR");
isSL[] := (SYMTYPE === "SL");
isSL3[] := (SYMTYPE === "SL3");
isP[] := (SYMTYPE === "P");
isPP[] := (SYMTYPE === "PP");
isNONE[] := (SYMTYPE === "NONE");
isQST[] := (SYMTYPE === "QST");
isQSTZ[] := (SYMTYPE === "QSTZ");
isQSZTZ[] := (SYMTYPE === "QSZTZ");
isQJ[] := (SYMTYPE === "QJ");
isSPSU2T[] := (SYMTYPE === "SPSU2T");
isQSC3[] := (SYMTYPE === "QSC3");
isSPSU2C3[] := (SYMTYPE === "SPSU2C3");
(* TODO: remove the above! *)
is[sym_] := SYMTYPE === sym;
