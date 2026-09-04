Module[{t = {}, floquetKets},
  floquetKets = Cases[bvc, ket[a___] :> {a}, Infinity];
  If[floquetKets === {} ||
     !And @@ Map[Length[#] == 3 && Take[#, 2] === {0, 0} &&
                   MemberQ[Range[-ncut, ncut], Last[#]]&, floquetKets],
    MyError["Floquet coordinate is not the final basis coordinate."]
  ];
  t = Join[t, mtSingletOp["m", opm]];
  t = Join[t, mtSingletOp["m^2", opm2]];
  t
]
