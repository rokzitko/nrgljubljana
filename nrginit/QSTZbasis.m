(* Highest-weight (Q, 2S+1, Tz) basis for spin-1/2 fermion orbitals.
   Q is measured from half filling; Tz = Sum[t_i n_i] is not shifted or
   doubled. Charges are explicit rules for whole orbitals, including indices.
   Loading this file outside SYMTYPE="QSTZ" only defines the constructor. *)
qstzbasisvc::metadata = "Supply exactly one integer-charge rule for each distinct orbital in the operator list.";
qstzbasisvc::tz = "A spin-adapted state does not have a unique integer Tz.";
qstzbasisvc::dimension = "The QSTZ basis does not span the fermionic Hilbert space.";

qstzbasisvc[ops_List, charges_] := Module[{keys, mask, qs, lift, result},
  If[!MatchQ[charges, {___Rule}],
    Message[qstzbasisvc::metadata]; Return[$Failed]
  ];
  keys = First /@ charges;
  If[!DuplicateFreeQ[ops] || !DuplicateFreeQ[keys] ||
      Sort[keys] =!= Sort[ops] || !VectorQ[Last /@ charges, IntegerQ],
    Message[qstzbasisvc::metadata]; Return[$Failed]
  ];
  qs = qsbasisvc[ops]; (* Also rejects nonfermionic and non-spin-1/2 orbitals. *)
  If[qs === $Failed, Return[$Failed]];
  mask = (makebasis[ops] /. head_[CR, indices___, spin_] :> head[indices]) /. charges;
  lift[qn_, state_] := Module[{tz},
    tz = Union[Cases[state, v_vc :> mask . (List @@ v), {0, Infinity}]];
    If[!MatchQ[tz, {_Integer}], Return[$Failed]];
    {{qn[[1]], 2 qn[[2]] + 1, First[tz]}, {state}}
  ];
  result = Flatten[Table[lift[sector[[1]], state],
    {sector, qs}, {state, sector[[2]]}], 1];
  If[!FreeQ[result, $Failed], Message[qstzbasisvc::tz]; Return[$Failed]];
  result = mergebasis[result];
  If[Total[(#[[1, 2]] Length[#[[2]]]) & /@ result] =!= 4^Length[ops],
    Message[qstzbasisvc::dimension]; Return[$Failed]
  ];
  result
];

(* The initializer's QSTZ writer packs f[+1], f[0], f[-1] into f 0 0.
   Models must already use these labels in H; only the generic basis operator
   list from initial.m is replaced here. Bath-only models declare {} metadata. *)
If[SYMTYPE === "QSTZ", Module[{},
  If[!TrueQ[CHANNELS == 3],
    MyError["QSTZ requires exactly three channels."]; Return[$Failed]
  ];
  If[!TrueQ[Ninit == 0],
    MyError["QSTZ supports only Ninit=0."]; Return[$Failed]
  ];
  If[option["READBASIS"],
    MyError["QSTZ does not support READBASIS; generate the basis instead."]; Return[$Failed]
  ];
  If[paramdefaultbool["spinless", False],
    MyError["QSTZ requires spin-1/2 fermions; spinless=true is unsupported."]; Return[$Failed]
  ];
  If[{MAKESPINKET, MAKEORBKET, MAKEPHONON, MAKEPHONONTENSOR, MAKEFLOQUET} =!=
      {Null, Null, Null, Null, Null},
    MyError["QSTZ does not support auxiliary spin, orbital, phonon, or Floquet kets."];
    Return[$Failed]
  ];
  If[!MatchQ[qstzImpurityCharges, {___Rule}] ||
      Sort[First /@ qstzImpurityCharges] =!= Sort[basopsdot],
    MyError["QSTZ requires qstzImpurityCharges rules for exactly basopsdot ({} for a bath-only model)."];
    Return[$Failed]
  ];
  basopsch = {f[1], f[0], f[-1]};
  basisops = basops = Sort[Join[basopsdot, basopsch]];
  NROPS = Length[basisops];
  bvc = qstzbasisvc[basisops,
    Join[qstzImpurityCharges, {f[1] -> 1, f[0] -> 0, f[-1] -> -1}]];
  If[bvc === $Failed, MyError["QSTZ basis construction failed."]; Return[$Failed]];
  vak = vacuum[];
  bz = bzvc2bzop[bvc];
  degnr[{q_Integer, ss_Integer, tz_Integer, rest___}] := ss;
]];
