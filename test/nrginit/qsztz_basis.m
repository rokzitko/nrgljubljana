$HistoryLength = 0;
$IterationLimit = Infinity;
package = Environment["NRG_INIT_DIR"];
Quiet[Get[FileNameJoin[{package, "sneg.m"}]], General::shdw];
Get[FileNameJoin[{package, "QSTZbasis.m"}]];
Get[FileNameJoin[{package, "QSZTZbasis.m"}]];
check[label_, ok_] := If[!TrueQ[ok], Print["FAIL ", label]; Exit[1]];

snegfermionoperators[f, impurity];
bath = {f[1], f[0], f[-1]};
checkBasis[ops_, tz_] := Module[{rules, b, states, qstz, expanded, labels},
  rules = Reverse[Thread[ops -> tz]];
  b = qsztzbasisvc[ops, rules];
  states = Join @@ b[[All, 2]];
  check["each Fock state appears once", Sort[states] === Sort[vc @@@ Tuples[{0, 1}, 2 Length[ops]]]];
  labels[v_vc] := With[{n = Partition[List @@ v, 2]},
    {Total[Flatten[n]] - Length[ops], 1 + Total[n . {1, -1}], tz . (Total /@ n)}];
  check["Q, 1+2Sz, unshifted Tz of every state",
    And @@ (Union[labels /@ Last[#]] === {First[#]} & /@ b)];

  (* Each QSTZ highest weight represents all 2S+1 spin projections. *)
  qstz = qstzbasisvc[ops, rules];
  expanded = Flatten[Table[With[{q = sector[[1, 1]], ss = sector[[1, 2]], tz0 = sector[[1, 3]]},
    Table[{q, ssz, tz0}, {ssz, 2 - ss, ss, 2}, {Length[Last[sector]]}]], {sector, qstz}], 2];
  check["QSZTZ sector dimensions vs expanded QSTZ",
    Sort[Tally[expanded]] === Sort[{First[#], Length[Last[#]]} & /@ b]];
  True
];
check["bath", checkBasis[bath, {1, 0, -1}]];
check["impurity + bath", checkBasis[Join[{impurity[], impurity[2], impurity["edge", 4]}, bath], {0, -1, 2, 1, 0, -1}]];
check["vacuum", qsztzbasisvc[{}, {}] === {{{0, 1, 0}, {vc[]}}}];
snegspinlessfermionoperators[spinless];
check["spinless rejected", Quiet[qsztzbasisvc[{spinless[]}, {spinless[] -> 0}]] === $Failed];
Do[check["invalid charge metadata", Quiet[qsztzbasisvc[bath, rules]] === $Failed],
  {rules, {{}, {f[1] -> 1, f[1] -> 0, f[-1] -> -1}, Thread[bath -> {1/2, 0, -1}]}}];

(* Exercise the nonempty model metadata as well as the standalone constructor. *)
SYMTYPE = "QSZTZ";
CHANNELS = 3; Ninit = 0;
option[_] := False;
paramdefaultbool[_, default_] := default;
MAKESPINKET = MAKEORBKET = MAKEPHONON = MAKEPHONONTENSOR = MAKEFLOQUET = Null;
basopsdot = {impurity["edge", 4]};
qsztzImpurityCharges = {impurity["edge", 4] -> 2};
Get[FileNameJoin[{package, "QSZTZbasis.m"}]];
expectedOps = Sort[Join[basopsdot, bath]];
check["adapter charges and multiplicity", basopsch === bath && basisops === expectedOps &&
  bvc === qsztzbasisvc[expectedOps, Join[qsztzImpurityCharges, Thread[bath -> {1, 0, -1}]]] &&
  And @@ (degnr[First[#]] === 1 & /@ bvc)];

(* The cache path in basis.m precedes module dispatch: it must not be touched. *)
option[name_] := name === "READBASIS";
silentGet[___] := (Print["FAIL touched basis cache"]; Exit[1]);
MyError[msg_String] := Throw[msg, "InitializerError"];
error = Catch[Get[FileNameJoin[{package, "basis.m"}]]; "no error", "InitializerError"];
check["READBASIS guard", StringContainsQ[error, "READBASIS"]];
Print["PASS QSZTZ basis"];
Exit[0];
