(*
   Single-frequency Floquet helpers for NRG Ljubljana.

   The Floquet index is the final ket/bra coordinate and runs from
   -Ncut through +Ncut. Shift operators have hard boundaries.
*)

floquetbasis::usage =
  "floquetbasis[Ncut] returns ket[-Ncut] through ket[Ncut].";
floquetid::usage =
  "floquetid[Ncut] returns the identity in the truncated Floquet space.";
floquetnumber::usage =
  "floquetnumber[Ncut] returns the Floquet-index operator.";
floquetshift::usage =
  "floquetshift[k,Ncut] shifts the final Floquet coordinate by k without wrapping.";
floquetplus::usage =
  "floquetplus[Ncut] is floquetshift[1,Ncut].";
floquetminus::usage =
  "floquetminus[Ncut] is floquetshift[-1,Ncut].";
floquetcos::usage =
  "floquetcos[Ncut] is floquetplus[Ncut]+floquetminus[Ncut].";
floquetlift::usage =
  "floquetlift[op,k,Ncut] tensors op with a shift of the final Floquet coordinate.";
transformtoFL::usage =
  "transformtoFL[basis,{Ncut}] appends a Floquet coordinate to basis.";

floquetbasis[Ncut_Integer] /; Ncut >= 0 :=
  Table[ket[i], {i, -Ncut, Ncut}];
floquetbasis[{Ncut_Integer}] /; Ncut >= 0 := floquetbasis[Ncut];

floquetcomplete[head_, args_List, coords_List] := Module[{padded},
  If[Length[args] > Length[coords],
    MyError["Floquet operator has more ket/bra coordinates than the basis."];
    Return[$Failed]
  ];
  padded = Join[ConstantArray[Null, Length[coords]-Length[args]], args];
  head @@ MapThread[If[#1 === Null, #2, #1]&, {padded, coords}]
];

floquetresolve[op_, coords_List] := op /. {
  ket[a___] :> floquetcomplete[ket, {a}, coords],
  bra[a___] :> floquetcomplete[bra, {a}, coords]
};

floquetembedded[x_?NumberQ] := x;
conj[floquetembedded[op_]] := floquetembedded[conj[op]];
ap[floquetembedded[op_], state:vc[prefix___, ket[coords___]]] :=
  ap[floquetresolve[op, {coords}], state];
ap[floquetembedded[_], _vc] := (
  MyError["A Floquet operator requires a trailing Floquet ket coordinate."];
  0
);
nc[a___, floquetembedded[op_], b___, state:ket[coords___]] /;
    FreeQ[{b}, _floquetembedded | _ket | _bra] :=
  nc[a, floquetresolve[op, {coords}], b, state];

floquetshiftbare[k_Integer, Ncut_Integer] /; Ncut >= 0 :=
  Module[{first, last},
    first = Max[-Ncut, -Ncut-k];
    last = Min[Ncut, Ncut-k];
    If[first > last,
      0,
      Sum[nc[ket[i+k], bra[i]], {i, first, last}]
    ]
  ];

floquetshift[k_Integer, Ncut_Integer] /; Ncut >= 0 :=
  floquetembedded[floquetshiftbare[k, Ncut]];
floquetid[Ncut_Integer] /; Ncut >= 0 := floquetshift[0, Ncut];
floquetnumber[Ncut_Integer] /; Ncut >= 0 :=
  floquetembedded[Sum[i nc[ket[i], bra[i]], {i, -Ncut, Ncut}]];
floquetplus[Ncut_Integer] /; Ncut >= 0 := floquetshift[1, Ncut];
floquetminus[Ncut_Integer] /; Ncut >= 0 := floquetshift[-1, Ncut];
floquetcos[Ncut_Integer] /; Ncut >= 0 :=
  floquetplus[Ncut] + floquetminus[Ncut];

floquetlift[op_, k_Integer, Ncut_Integer] /; Ncut >= 0 :=
  floquetembedded[ketbratensorproduct[op, floquetshiftbare[k, Ncut]]];

transformtoFL[basis_List, {Ncut_Integer}] /; Ncut >= 0 :=
  Module[{vecs},
    vecs = Map[
      Function[floquetstate,
        applybasis[basis,
          Function[state, ketbratensorproduct[state, floquetstate]]]],
      floquetbasis[Ncut]];
    mergebasis @ Flatten[vecs, 1]
  ];
transformtoFL[basis_List, Ncut_Integer] /; Ncut >= 0 :=
  transformtoFL[basis, {Ncut}];

floquetvalidatebasis[basis_List, Ncut_Integer] /; Ncut >= 0 :=
  Module[{states, terms, coords, widths, width},
    If[!TrueQ[bzQ[basis]] || basis === {} ||
        !And @@ Map[Last[#] =!= {}&, basis] || !FreeQ[basis, Null],
      MyError["Floquet basis is not a canonical occupation-vector basis."];
      Return[$Failed]
    ];
    states = Flatten[basis[[All, 2]], 1];
    terms = Flatten[Map[sum2list, Expand /@ states], 1];
    If[terms === {} || MemberQ[terms, 0] ||
        !And @@ Map[
          Count[#, HoldPattern[ket[___]], {0, Infinity}] == 1&, terms],
      MyError["Every Floquet basis state must contain one nonzero auxiliary ket."];
      Return[$Failed]
    ];
    coords = Cases[terms, ket[a___] :> {a}, Infinity];
    widths = DeleteDuplicates[Length /@ coords];
    If[coords === {} || Length[widths] =!= 1,
      MyError["Floquet basis states must have a uniform trailing mode coordinate."];
      Return[$Failed]
    ];
    width = First[widths];
    If[!And @@ Map[
        IntegerQ[Last[#]] && -Ncut <= Last[#] <= Ncut&, coords],
      MyError["Floquet basis states have an invalid final mode coordinate."];
      Return[$Failed]
    ];
    width
  ];

floquetpreparebasis[basis_, vacuum_] := Module[{states, opbasis},
  If[!ListQ[basis] || !TrueQ[bzQ[basis]] || basis === {} ||
      !And @@ Map[Last[#] =!= {}&, basis] || !FreeQ[basis, Null],
    MyError["hook_basis returned a noncanonical occupation-vector basis."];
    Return[$Failed]
  ];
  states = Flatten[basis[[All, 2]], 1];
  If[states === {} || MemberQ[states, 0],
    MyError["hook_basis returned a noncanonical occupation-vector basis."];
    Return[$Failed]
  ];
  opbasis = bzvc2bzop[basis];
  If[!FreeQ[opbasis, _vc | _vc2ops | _ap | _DirectedInfinity] ||
      bzop2bzvc[opbasis, vacuum] =!= basis,
    MyError["hook_basis returned a noncanonical occupation-vector basis."];
    Return[$Failed]
  ];
  opbasis
];

(* Exact parsing for the public options=Ncut=... compatibility syntax. *)
floquetoptionvalue[keyword_String, options_List] :=
  Module[{prefix, matches},
    prefix = keyword <> "=";
    matches = Select[options,
      StringLength[#] >= StringLength[prefix] &&
      StringTake[#, StringLength[prefix]] === prefix&];
    If[Length[matches] =!= 1,
      MyError["Expected exactly one option ", prefix, "VALUE."];
      Return[$Failed]
    ];
    StringDrop[First[matches], StringLength[prefix]]
  ];

parsefloquetcutoff[options_List] := Module[{raw, chars},
  raw = floquetoptionvalue["Ncut", options];
  If[raw === $Failed, Return[$Failed]];
  chars = Characters[raw];
  If[chars === {} || !(And @@ Map[DigitQ, chars]),
    MyError["Ncut must be a non-negative decimal integer, got: ", raw];
    Return[$Failed]
  ];
  FromDigits[Map[ToExpression, chars]]
];

floquetdropsign[chars_List] :=
  If[chars =!= {} && MemberQ[{"+", "-"}, First[chars]], Rest[chars], chars];

floquetunsigneddecimalq[chars_List] :=
  chars =!= {} && Count[chars, "."] <= 1 &&
  Count[chars, _?DigitQ] >= 1 &&
  And @@ Map[(DigitQ[#] || # === ".")&, chars];

floquetsignedintegerq[chars_List] := Module[{digits},
  digits = floquetdropsign[chars];
  digits =!= {} && And @@ Map[DigitQ, digits]
];

floquetnumericstringq[value_String] :=
  Module[{chars, exponentPositions, mantissa, exponent},
    chars = Characters[stripws[value]];
    If[chars === {}, Return[False]];
    exponentPositions = Flatten @ Position[chars, "e" | "E"];
    If[Length[exponentPositions] > 1, Return[False]];
    If[exponentPositions === {},
      Return[floquetunsigneddecimalq[floquetdropsign[chars]]]
    ];
    mantissa = Take[chars, exponentPositions[[1]]-1];
    exponent = Drop[chars, exponentPositions[[1]]];
    floquetunsigneddecimalq[floquetdropsign[mantissa]] &&
      floquetsignedintegerq[exponent]
  ];

floquetexactdecimal[raw_String] :=
  Module[{chars, sign = 1, exponentPositions, exponent = 0,
      exponentChars, mantissa, dotPositions, fractionalDigits = 0,
      digits, significand, power, order},
    chars = Characters[stripws[raw]];
    If[First[chars] === "+", chars = Rest[chars]];
    If[First[chars] === "-", sign = -1; chars = Rest[chars]];
    exponentPositions = Flatten @ Position[chars, "e" | "E"];
    If[exponentPositions =!= {},
      exponentChars = Drop[chars, First[exponentPositions]];
      mantissa = Take[chars, First[exponentPositions]-1];
      If[First[exponentChars] === "+", exponentChars = Rest[exponentChars]];
      If[First[exponentChars] === "-",
        exponentChars = Rest[exponentChars];
        exponent = -FromDigits[Map[ToExpression, exponentChars]],
        exponent = FromDigits[Map[ToExpression, exponentChars]]
      ],
      mantissa = chars
    ];
    dotPositions = Flatten @ Position[mantissa, "."];
    If[dotPositions =!= {},
      fractionalDigits = Length[mantissa]-First[dotPositions]
    ];
    digits = DeleteCases[mantissa, "."];
    significand = FromDigits[Map[ToExpression, digits]];
    If[significand == 0, Return[0]];
    power = exponent-fractionalDigits;
    order = IntegerLength[significand]+power-1;
    If[order > 308 || order < -324, Return[$Failed]];
    sign If[power >= 0,
      significand 10^power,
      significand/10^-power
    ]
  ];

floquetroundeven[value_] := Module[{lower, remainder},
  lower = Floor[value];
  remainder = value-lower;
  Which[
    remainder < 1/2, lower,
    remainder > 1/2, lower+1,
    EvenQ[lower], lower,
    True, lower+1
  ]
];

floquetmachinedouble[value_] :=
  Module[{binaryExponent, significand, scale},
    If[!TrueQ[value > 0], Return[$Failed]];
    If[value < 2^-1022,
      significand = floquetroundeven[value 2^1074];
      If[significand == 0, Return[$Failed]];
      scale = Quiet[N[2^-1074, MachinePrecision], General::munfl];
      Return[N[significand, MachinePrecision] scale]
    ];
    binaryExponent = IntegerLength[Numerator[value], 2] -
      IntegerLength[Denominator[value], 2];
    If[value < 2^binaryExponent, binaryExponent--];
    significand = floquetroundeven[value 2^(52-binaryExponent)];
    If[significand == 2^53,
      significand = 2^52;
      binaryExponent++
    ];
    If[binaryExponent > 1023, Return[$Failed]];
    scale = Quiet[N[2^(binaryExponent-52), MachinePrecision],
      General::munfl];
    N[significand, MachinePrecision] scale
  ];

parsefloquetomega[raw_String] := Module[{exactValue, value},
  If[!floquetnumericstringq[raw],
    MyError["Omega must be a finite positive double-precision number, got: ", raw];
    Return[$Failed]
  ];
  exactValue = floquetexactdecimal[raw];
  If[exactValue === $Failed,
    MyError["Omega must be a finite positive double-precision number, got: ", raw];
    Return[$Failed]
  ];
  value = floquetmachinedouble[exactValue];
  If[!MachineNumberQ[value] || !TrueQ[0. < value <= $MaxMachineNumber],
    MyError["Omega must be a finite positive double-precision number, got: ", raw];
    Return[$Failed]
  ];
  value
];

validatefloquetmodedeclarations[blocks_List] := Module[{headers},
  headers = Select[Cases[blocks, {header_String} :> header],
    StringLength[#] == 2 && StringTake[#, -1] === "m" &&
      MemberQ[{"s", "p", "g", "d", "t", "o", "q"},
        StringTake[#, 1]]&];
  If[headers =!= {"sm"},
    MyError["Floquet input must contain exactly one singlet operator named m."];
    Return[$Failed]
  ];
  True
];
