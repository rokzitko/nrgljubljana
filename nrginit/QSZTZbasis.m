(* Full occupation-vector basis labelled by (Q, 1+2 Sz, Tz).
   Q is relative to half filling; Tz = Sum[t_i n_i] is not shifted or doubled.
   Loading this file outside SYMTYPE="QSZTZ" only defines the constructor. *)
qsztzbasisvc::metadata = "Supply exactly one integer-charge rule for each distinct orbital in the operator list.";
qsztzbasisvc::spin = "QSZTZ requires spin-1/2 fermion operators. Unsupported operators: ``.";
qsztzbasisvc::dimension = "The QSZTZ basis does not span the fermionic Hilbert space.";

qsztzbasisvc[ops_List, charges_] := Module[{keys, bad, qsz, mask, result},
  If[!MatchQ[charges, {___Rule}],
    Message[qsztzbasisvc::metadata]; Return[$Failed]
  ];
  keys = First /@ charges;
  If[!DuplicateFreeQ[ops] || !DuplicateFreeQ[keys] ||
      Sort[keys] =!= Sort[ops] || !VectorQ[Last /@ charges, IntegerQ],
    Message[qsztzbasisvc::metadata]; Return[$Failed]
  ];
  bad = qsbasisBadOps[ops]; (* qszbasisvc itself is not limited to spin-1/2. *)
  If[bad =!= {}, Message[qsztzbasisvc::spin, bad]; Return[$Failed]];
  If[ops === {}, makebasis[{}]; Return[{{{0, 1, 0}, {vacuum[]}}}]];

  qsz = qszbasisvc[ops]; (* Sets BASIS in the order used by the occupation vectors. *)
  mask = (BASIS /. head_[CR, indices___, spin_] :> head[indices]) /. charges;
  If[Length[mask] =!= 2 Length[ops] || !VectorQ[mask, IntegerQ],
    Message[qsztzbasisvc::metadata]; Return[$Failed]
  ];
  result = mergebasis[Flatten[Table[
    {{sector[[1, 1]], 1 + 2 sector[[1, 2]], mask . (List @@ v)}, {v}},
    {sector, qsz}, {v, sector[[2]]}], 1]];
  If[Total[Length /@ result[[All, 2]]] =!= 4^Length[ops],
    Message[qsztzbasisvc::dimension]; Return[$Failed]
  ];
  result
];

(* The writer packs f[+1], f[0], f[-1] into f 0 0. The Hamiltonian must
   already use these labels; only the generic bath basis inventory changes. *)
If[SYMTYPE === "QSZTZ", Module[{},
  If[!TrueQ[CHANNELS == 3],
    MyError["QSZTZ requires exactly three channels."]; Return[$Failed]
  ];
  If[!TrueQ[Ninit == 0],
    MyError["QSZTZ supports only Ninit=0."]; Return[$Failed]
  ];
  If[option["READBASIS"],
    MyError["QSZTZ does not support READBASIS; generate the basis instead."]; Return[$Failed]
  ];
  If[paramdefaultbool["spinless", False],
    MyError["QSZTZ requires spin-1/2 fermions; spinless=true is unsupported."]; Return[$Failed]
  ];
  If[{MAKESPINKET, MAKEORBKET, MAKEPHONON, MAKEPHONONTENSOR, MAKEFLOQUET} =!=
      {Null, Null, Null, Null, Null},
    MyError["QSZTZ does not support auxiliary spin, orbital, phonon, or Floquet kets."];
    Return[$Failed]
  ];
  If[!MatchQ[qsztzImpurityCharges, {___Rule}] ||
      Sort[First /@ qsztzImpurityCharges] =!= Sort[basopsdot],
    MyError["QSZTZ requires qsztzImpurityCharges rules for exactly basopsdot ({} for a bath-only model)."];
    Return[$Failed]
  ];
  basopsch = {f[1], f[0], f[-1]};
  basisops = basops = Sort[Join[basopsdot, basopsch]];
  NROPS = Length[basisops];
  bvc = qsztzbasisvc[basisops,
    Join[qsztzImpurityCharges, {f[1] -> 1, f[0] -> 0, f[-1] -> -1}]];
  If[bvc === $Failed, MyError["QSZTZ basis construction failed."]; Return[$Failed]];
  vak = vacuum[];
  bz = bzvc2bzop[bvc];
  degnr[{q_Integer, ssz_Integer, tz_Integer, rest___}] := 1;
]];
