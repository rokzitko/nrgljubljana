Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];

loadWilson[{{"tri", "rkpw"}, {"tridiag_method", "rkpw"}, {"disccheck", "true"}}];
analytic = Table[(1 + 1/2)/2 2^(-n/2) (1 - 2^(-n-1))/Sqrt[(1 - 2^(-2n-1)) (1 - 2^(-2n-3))], {n, 0, 4}];
check["initializer flat-band analytic oracle", close[Flatten[xitable[1]], analytic] && close[Flatten[zetatable[1]], ConstantArray[0, 5]]];
check["exact particle-hole symmetry in initializer tables", Flatten[zetatable[1]] === ConstantArray[0., 5]];
check["upstream precision remains separate", PREC == 30 && Precision[du[1][0, 0]] > 20 && Precision[de[1, 0]] > 20];
check["machine coefficients eagerly populated", And @@ (MachineNumberQ /@ Flatten[{Table[xi[1][n], {n, 0, 4}], Table[dzeta[1][n], {n, 0, 4}]}])];
check["no higher vectors or squared hoppings", !NumericQ[du[1][1, 0]] && !NumericQ[dv[1][1, 0]] && SubValues[xi2] == {}];
Block[{rkpwScalar, xi, dzeta},
  rkpwScalar[input_, n_] := (observedPoles = input; ConstantArray[0., {2, n}]);
  dothelanczosrkpw[];
];
check["initializer shell interleaving", observedPoles === Flatten[Table[{{de[1, m], du[1][0, m]}, {-deminus[1, m], dv[1][0, m]}}, {m, 0, mMAX}], 1]];
fullCoefficients = {Flatten[zetatable[1]], Flatten[xitable[1]]};
Ninit = 2;
loadWilson[{{"tri", "cpp"}, {"tridiag_method", "rkpw"}, {"disccheck", "true"}}];
check["cpp rkpw Ninit seed", RKPW && DISCNMAX == 2 && TRUEDISCNMAX == 4 && close[{Flatten[zetatable[1]], Flatten[xitable[1]]}, Take[#, 3] & /@ fullCoefficients]];
loadWilson[{{"tri", "none"}, {"tridiag_method", "rkpw"}}];
check["none rkpw Ninit seed", RKPW && DISCNMAX == 2 && close[{Flatten[zetatable[1]], Flatten[xitable[1]]}, Take[#, 3] & /@ fullCoefficients]];
Block[{hookfile},
  hookfile["hook_post_lanczosinit"] := (ClearAll[deminus]; deminus[a_, m_] := de[a, m]/2);
  loadWilson[{{"tri", "rkpw"}}];
  check["asymmetric initializer onsite is not removed", Abs[dzeta[1][0]] > 0.01];
];
Block[{hookfile},
  nearAmplitude = 1. + 4 2.^-52;
  nearChain = rkpwScalar[{{1., 1.}, {-1., nearAmplitude}}, 1];
  hookfile["hook_post_lanczosinit"] := (ClearAll[de, deminus, du, dv];
    de[_, _] = deminus[_, _] = 1.;
    du[_][0, m_] := If[m == 0, 1., 0.];
    dv[_][0, m_] := If[m == 0, nearAmplitude, 0.]);
  loadWilson[{{"tri", "rkpw"}, {"nrxi", "0"}}];
  check["near-machine asymmetry is not thresholded",
    nearChain[[1, 1]] != 0. && Abs[nearChain[[1, 1]]] < 10^-14 && dzeta[1][0] === nearChain[[1, 1]]];
];
COEFCHANNELS = CHANNELS = 2;
loadWilson[{{"tri", "rkpw"}}];
check["independent scalar channels", close[xitable[1], xitable[2]] && MachineNumberQ[xi[2][4]] && Precision[dv[2][0, 0]] > 20];
COEFCHANNELS = CHANNELS = 1;
expectError["unknown tri", loadWilson[{{"tri", "typo"}}], "Unknown tri backend: typo"];
expectError["unknown method with full tri", loadWilson[{{"tri", "rkpw"}, {"tridiag_method", "typo"}}], "Unknown tridiag_method backend: typo"];
expectError["unknown method with cpp", loadWilson[{{"tri", "cpp"}, {"tridiag_method", "typo"}}], "Unknown tridiag_method backend: typo"];
POL2x2 = True;
expectError["no matrix extension", loadWilson[{{"tri", "rkpw"}}], "scalar chains"];
POL2x2 = False;
RUNGS = True;
expectError["no rung extension", loadWilson[{{"tri", "rkpw"}}], "scalar chains"];
RUNGS = False;
expectError["no matrix interface", loadWilson[{{"tri", "rkpw"}, {"wilsonchain", "matrix"}}], "scalar chains"];
Block[{BAND = "nambu"}, expectError["no Nambu band", loadWilson[{{"tri", "rkpw"}}], "scalar chains"]];
Block[{isSC, SYMTYPE, COEFCHANNELS, CHANNELS, BAND},
  isSC[] := True; BAND = "flat"; COEFCHANNELS = CHANNELS = 1;
  Do[
    SYMTYPE = symmetry;
    Do[
      loadWilson[{{"tri", "rkpw"}, {"tridiag_method", "rkpw"}, {"bcsgap", gap}}];
      check["prescribed scalar pairing " <> symmetry <> " " <> gap,
        close[Flatten[xitable[1]], analytic] && close[Flatten[scdeltatable[1]], ConstantArray[ToExpression[gap], 5]] &&
        Flatten[sckappatable[1]] == ConstantArray[0, 5]],
      {gap, {"0", "0.125"}}],
    {symmetry, {"SPSU2", "SPU1", "SPU1LR", "P", "PP", "NONE"}}];
  SYMTYPE = "SPSU2";
  Do[
    COEFCHANNELS = CHANNELS = channels;
    pairs = Join[{{"tri", "rkpw"}}, Table[{"bcsgap" <> ToString[a], "!" <> ToString[a/8, InputForm]}, {a, channels}]];
    loadWilson[pairs];
    check["channel-complete prescribed pairing " <> ToString[channels],
      And @@ Table[close[Flatten[scdeltatable[a]], ConstantArray[a/8, 5]], {a, channels}]],
    {channels, {2, 3}}];
  expectError["incomplete prescribed pairing", loadWilson[{{"tri", "rkpw"}, {"bcsgap1", "0"}, {"bcsgap2", "0"}}], "channel-complete"];
  COEFCHANNELS = CHANNELS = 1;
  expectError["missing prescribed pairing", loadWilson[{{"tri", "rkpw"}}], "explicit finite constant"];
  Do[expectError["invalid prescribed pairing " <> gap, loadWilson[{{"tri", "rkpw"}, {"bcsgap", gap}}], "explicit finite constant"],
    {gap, {"!Infinity", "!Indeterminate", "!I", "!omega"}}];
  Do[expectError["no SC runtime handoff " <> mode,
    loadWilson[{{"tri", mode}, {"tridiag_method", "rkpw"}, {"bcsgap", "0"}}], "pairing-table handoff"], {mode, {"cpp", "none"}}];
  BAND = "cosine";
  expectError["no non-flat SC extension", loadWilson[{{"tri", "rkpw"}, {"bcsgap", "0"}}], "flat-band"];
  BAND = "flat"; SYMTYPE = "SPSU2T";
  expectError["no unaudited SC extension", loadWilson[{{"tri", "rkpw"}, {"bcsgap", "0"}}], "audited scalar"];
  SYMTYPE = "SPSU2";
  Block[{hookfile},
    hookfile["hook_bcs"] := (gapdefined = True);
    expectError["no custom pairing tables", loadWilson[{{"tri", "rkpw"}, {"bcsgap", "0"}}], "hook-defined pairing"];
  ];
];

Do[expectError["invalid RKPW bandrescale " <> value,
  loadWilson[{{"tri", "rkpw"}, {"bandrescale", value}}], "bandrescale must be a finite positive machine real"],
  {value, {"0", "-1", "!Infinity", "!Indeterminate", "!I", "!10^400", "!10^-400"}}];
expectError["scaled hopping underflow", loadWilson[{{"tri", "rkpw"}, {"bandrescale", "!2^-1074"}, {"gap", "1"}}], "final scaled coefficient"];
expectError["cpp seed scaling checked", loadWilson[{{"tri", "cpp"}, {"tridiag_method", "rkpw"}, {"bandrescale", "!2^-1074"}, {"gap", "1"}}], "final scaled coefficient"];
loadWilson[{{"tri", "rkpw"}, {"bandrescale", "!10^-310"}, {"gap", "1"}}];
scaledTables = {Flatten[zetatable[1]], Flatten[xitable[1]]};
rawTables = {Table[zeta[1][n], {n, 0, 4}], Table[xi[1][n], {n, 0, 4}]};
check["representable subnormal output tables", And @@ (MachineNumberQ /@ Flatten[scaledTables]) &&
  Min[Flatten[xitable[1]]] > 0. && close[SetPrecision[scaledTables, Infinity]/SetPrecision[rkpwBandscale, Infinity], rawTables, 10^-12]];
expectError["final onsite includes shift0", loadWilson[{{"tri", "rkpw"}, {"bandrescale", "!2^1023"}, {"shift0", "4"}}], "final scaled coefficient"];
expectError["final onsite includes gap", loadWilson[{{"tri", "rkpw"}, {"bandrescale", "!2^1023"}, {"gap", "4"}}], "final scaled coefficient"];
expectError["nonfinite onsite rejected", loadWilson[{{"tri", "rkpw"}, {"shift0", "!Infinity"}}], "final scaled coefficient"];
Block[{POLARIZED = True, CHANNELS = 1, COEFCHANNELS = 2, isQSZ, isU1, isSPU1, isP, isPP, isNONE},
  isQSZ[] = isU1[] = isSPU1[] = isPP[] = isNONE[] = False; isP[] = True;
  expectError["final onsite includes globalh", loadWilson[{{"tri", "rkpw"}, {"bandrescale", "!2^1023"}, {"globalh", "8"}}], "final scaled coefficient"];
];
(* A one-pole star has an exact terminal hopping and zero raw onsite. *)
Block[{hookfile},
  hookfile["hook_post_lanczosinit"] := (ClearAll[du, dv];
    du[_][0, m_] := If[m == 0, 1, 0]; dv[_][0, _] = 0; de[_, 0] = 0);
  loadWilson[{{"tri", "rkpw"}, {"nrxi", "0"}, {"bandrescale", "!2^-1074"}, {"shift0", "1"}}];
  check["least subnormal onsite and exact terminal zero", SetPrecision[zetatable[1][[1, 1]], Infinity] == 2^-1074 && xitable[1] === {{0.}}];
  expectError["physical onsite shift underflow", loadWilson[{{"tri", "rkpw"}, {"nrxi", "0"}, {"bandrescale", "!10^-310"}, {"shift0", "!10^-20"}}], "final scaled coefficient"];
];
(* Independent 100-digit vector Lanczos oracle, used only for small stars. *)
reference[poles_, count_] := Module[{e, v, previous, residual, alpha, beta = 0, aa = {}, bb = {}, j},
  {e, v} = Transpose[N[poles, 100]];
  v = v/Sqrt[v.v]; previous = 0 v;
  Do[
    alpha = v.(e v);
    AppendTo[aa, alpha];
    If[j == Length[e], AppendTo[bb, 0]; Break[]];
    residual = e v - alpha v - beta previous;
    beta = Sqrt[residual.residual];
    AppendTo[bb, beta]; previous = v; v = residual/beta,
    {j, count}];
  {aa, bb}
];
star = {{9/10, 1/3}, {-7/10, 1/2}, {2/5, 2/3}, {-1/5, 1/5}, {1/20, 1/7}, {-1/30, 3/7}};
chain = rkpwScalar[star, 6];
check["asymmetric independent oracle", close[chain, reference[star, 6]]];
check["support exhaustion terminal zero", Last[chain[[2]]] === 0. && Min[Most[chain[[2]]]] > 0.];
Do[check["prefix incorporates all poles " <> ToString[n], close[rkpwScalar[star, n], Take[#, n] & /@ chain]], {n, 1, 5}];
jacobi = DiagonalMatrix[chain[[1]]] + DiagonalMatrix[Most[chain[[2]]], 1] + DiagonalMatrix[Most[chain[[2]]], -1];
{eigenvalues, eigenvectors} = Eigensystem[jacobi];
recovered = SortBy[Transpose[{eigenvalues, eigenvectors[[All, 1]]^2}], First];
expected = SortBy[Transpose[{star[[All, 1]], star[[All, 2]]^2/Total[star[[All, 2]]^2]}], First];
check["asymmetric spectral nodes and weights", close[recovered, expected]];

check["one pole", rkpwScalar[{{2, 3}}, 1] === {{2.}, {0.}}];
check["zero-energy pole", rkpwScalar[{{0, 1}}, 1] === {{0.}, {0.}}];
check["zeros stripped and exact duplicates merged", close[rkpwScalar[{{2, 0}, {1, 3}, {-2, 4}, {1, 4}, {7, 0}}, 2], reference[{{1, 5}, {-2, 4}}, 2]]];
check["merge after machine conversion", rkpwScalar[{{1, 1}, {1 + 10^-20, 1}}, 1] === {{1.}, {0.}}];
adjacent = rkpwScalar[{{1., 1.}, {1. + 2.^-52, 1.}}, 2];
check["adjacent doubles not merged", adjacent[[2, 1]] > 0. && adjacent[[2, 2]] === 0.];
expectError["empty input", rkpwScalar[{}, 1], "nonempty"];
expectError["zero support", rkpwScalar[{{0, 0}, {1, 0}}, 1], "support is 0"];
expectError["count exceeds effective support", rkpwScalar[{{1, 1}, {1, 2}, {2, 0}}, 2], "support is 1"];
expectError["zero count", rkpwScalar[star, 0], "positive integer"];
expectError["fractional count", rkpwScalar[star, 3/2], "positive integer"];
expectError["negative amplitude", rkpwScalar[{{1, -1}}, 1], "nonnegative"];
expectError["complex energy", rkpwScalar[{{I, 1}}, 1], "machine real"];
expectError["nonfinite input", rkpwScalar[{{Infinity, 1}}, 1], "machine real"];
expectError["energy underflow", rkpwScalar[{{10^-400, 1}}, 1], "rounded to zero"];
expectError["amplitude underflow", rkpwScalar[{{1, 10^-400}}, 1], "rounded to zero"];
expectError["energy overflow", rkpwScalar[{{10^400, 1}}, 1], "machine real"];
check["common scaling prevents merged amplitude overflow", rkpwScalar[{{1, 15 10^307}, {1, 15 10^307}}, 1] === {{1.}, {0.}}];
check["common scaling prevents norm overflow", close[rkpwScalar[{{1, 15 10^307}, {-1, 15 10^307}}, 2], {{0, 0}, {1, 0}}]];
expectError["common scaling must not erase a nonzero amplitude", rkpwScalar[{{1, 2^1023}, {-1, 2^-1074}}, 2], "during common power-of-two scaling"];
expectError["recurrence overflow fails without fallback", rkpwScalar[{{10^308, 1}, {-10^308, 1}}, 2], "numerical breakdown"];
expectError["unrepresentable hopping fails", rkpwScalar[{{0, 1}, {10^-200, 10^-200}}, 2], "underflowed hopping"];

(* Tiny amplitudes and energies must not be squared in the reconstruction. *)
tiny = rkpwScalar[{{10^-200, 10^-200}, {-10^-200, 10^-200}}, 2];
check["tiny unsquared hopping", Abs[tiny[[2, 1]]/1.*^-200 - 1] < 10^-14];
tinyMerged = rkpwScalar[{{1, 3 10^-200}, {-2, 5 10^-200}, {1, 4 10^-200}}, 2];
check["safe hypot in duplicate merging", close[tinyMerged, reference[{{1, 5}, {-2, 5}}, 2]]];
subnormal = rkpwScalar[{{10^-310, 1}, {-10^-310, 1}}, 2];
check["representable subnormal hopping", MachineNumberQ[subnormal[[2, 1]]] && Abs[subnormal[[2, 1]]/1.*^-310 - 1] < 10^-12];
subnormalAmplitudes = rkpwScalar[{{1, 10^-310}, {-1, 10^-310}}, 2];
check["subnormal amplitude ratios", close[subnormalAmplitudes, {{0, 0}, {1, 0}}, 10^-12]];
leastSubnormal = Quiet[N[2^-1074, MachinePrecision], General::munfl];
check["least subnormal equal amplitudes", close[rkpwScalar[{{1., leastSubnormal}, {-1., leastSubnormal}}, 2], {{0, 0}, {1, 0}}]];
(* A later rotation can have subnormal operands despite initial normalization.
   Its direction must not be computed from the already rounded physical norm. *)
Do[
  mixedAmplitude = Quiet[N[multiple 2^-1074, MachinePrecision], General::munfl];
  mixedChain = rkpwScalar[{{0., 1.}, {1., mixedAmplitude}, {-1., mixedAmplitude}}, 3];
  roundedFirstHop = Quiet[N[N[Sqrt[2] multiple 2^-1074, 80], MachinePrecision], General::munfl];
  check["later subnormal rotation " <> ToString[multiple],
    Max[Abs[mixedChain[[1]]]] < 10^-14 && Abs[mixedChain[[2, 2]] - 1.] < 10^-14 && mixedChain[[2, 3]] === 0.];
  check["rounded subnormal first hopping " <> ToString[multiple],
    Abs[SetPrecision[mixedChain[[2, 1]], Infinity]/SetPrecision[roundedFirstHop, Infinity] - 1] < 1/10];
  check["mixed extreme chain stays machine " <> ToString[multiple], And @@ (MachineNumberQ /@ Flatten[mixedChain])],
  {multiple, {1, 2, 3, 4, 7, 16, 255, 1024}}];
scaleStar = {{1, 1}, {-1, 2}, {1/4, 3}};
scaleReference = reference[scaleStar, 3];
Do[
  scaledStar = Quiet[N[Map[{#[[1]], #[[2]] multiple 2^exponent} &, scaleStar], MachinePrecision], General::munfl];
  check["common amplitude scale invariance " <> ToString[{multiple, exponent}], close[rkpwScalar[scaledStar, 3], scaleReference]],
  {multiple, {1, 2, 3}}, {exponent, {-1074, -1073, -1022, -500, 0, 500, 1020}}];
repeatedTiny = Quiet[N[{{1, 2^-1074}, {-1, 2 2^-1074}, {1, 3 2^-1074}}, MachinePrecision], General::munfl];
check["normalize before merging subnormal duplicates", close[rkpwScalar[repeatedTiny, 2], reference[{{1, Sqrt[10]}, {-1, 2}}, 2]]];

(* Wilson's flat-band analytic chain, Lambda=4, z=1, Yoshida mesh.
   Shell ordering matters: alternating signs from the outermost inward.
   900 hops reach ~10^-271; squared hops underflow much earlier. *)
shells = 510; count = 900;
flatStar = Flatten[Table[{{(1 + 1/4)/2 4^-m, Sqrt[(1 - 1/4) 4^-m/2]},
  {-(1 + 1/4)/2 4^-m, Sqrt[(1 - 1/4) 4^-m/2]}}, {m, 0, shells}], 1];
flatChain = rkpwScalar[N[flatStar, 40], count];
analytic = N[Table[(1 + 1/4)/2 4^(-n/2) (1 - 4^(-n-1))/Sqrt[(1 - 4^(-2n-1)) (1 - 4^(-2n-3))], {n, 0, count-1}], 40];
relativeError = Max[Abs[SetPrecision[flatChain[[2]], 40]/analytic - 1]];
onsiteError = Max[Abs[SetPrecision[flatChain[[1]], 40]/analytic]];
Print["Flat analytic maximum relative hopping error: ", relativeError, "; scaled onsite error: ", onsiteError];
check["long flat analytic tail", relativeError < 10^-11 && onsiteError < 10^-11 && Min[flatChain[[2]]] > 0.];
check["no arbitrary precision fallback", And @@ (MachineNumberQ /@ Flatten[flatChain])];
check["short flat prefix uses all shells", close[rkpwScalar[N[flatStar, 40], 3], Take[#, 3] & /@ flatChain]];

Print["rkpw failures: ", failures];
If[failures == 0, True, $Failed]
