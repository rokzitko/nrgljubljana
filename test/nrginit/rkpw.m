(* Standalone numerical and wilson.m dispatch tests: no SNEG or model generation. *)
DEBUG = 0;
Get[FileNameJoin[{sourceDir, "nrginit", "misc.m"}]];
Get[FileNameJoin[{sourceDir, "nrginit", "initialparse.m"}]];
MyError[args__] := Throw[StringJoin[ToString /@ {args}], "initializer-error"];
hook[_] := Null;
hookfile[_] := Null;
option[_] := False;
isSC[] := False;
SYMTYPE = "QS";
COEFCHANNELS = CHANNELS = 1;
POLARIZED = POL2x2 = RUNGS = False;
BAND = "flat";
DY = True; DC = DZ = False;
lambda = 2.; z = 1.; bandrescale = 1;
Ninit = 0;

failures = 0;
check[label_, condition_] := If[!TrueQ[condition], Print["FAILED: ", label]; failures++];
SetAttributes[expectError, HoldRest];
expectError[label_, expression_, fragment_] := Module[{result},
  result = Catch[expression, "initializer-error"];
  check[label, StringQ[result] && StringContainsQ[result, fragment]];
];
close[x_, y_, tolerance_:10^-13] := Max[Abs[Flatten[x - y]]] < tolerance;

loadWilson[pairs_] := Module[{},
  ClearAll[data, de, deminus, df, dfminus, eps, thetaCh, demem, deminusmem, diagA, zeta,
    du0, dv0, uvrescalefactor, xitable, zetatable, eptable, emtable, u0ptable, u0mtable, i, m];
  listkeywords["param"] = First /@ pairs;
  listkeywords["dmft"] = {};
  Scan[(data["param"][#[[1]]] = #[[2]]) &, pairs];
  bandrescale = paramdefaultnum["bandrescale", 1];
  Nmax = 4;
  Get[FileNameJoin[{sourceDir, "nrginit", "wilson.m"}]];
];

loadWilson[{}];
check["legacy defaults", TRI == "old" && TRIDIAGMETHOD == "lanczos" && PREC == 1000 && !RKPW];
oldCoefficients = {Flatten[zetatable[1]], Flatten[xitable[1]]};
loadWilson[{{"tri", "rkpw"}, {"disccheck", "true"}}];
check["rkpw agrees with short high-precision legacy chain", close[{Flatten[zetatable[1]], Flatten[xitable[1]]}, oldCoefficients]];
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
loadWilson[{{"tri", "cpp"}}];
check["cpp default seed unchanged", !RKPW && dothelanczos === dothelanczosold && PREC == 30 && !MachineNumberQ[xi[1][0]]];
loadWilson[{{"tri", "none"}}];
check["none default seed unchanged", !RKPW && dothelanczos === dothelanczosold];
COEFCHANNELS = CHANNELS = 2;
loadWilson[{{"tri", "rkpw"}}];
check["independent scalar channels", close[xitable[1], xitable[2]] && MachineNumberQ[xi[2][4]] && Precision[dv[2][0, 0]] > 20];
COEFCHANNELS = CHANNELS = 1;
expectError["unknown tri", loadWilson[{{"tri", "typo"}}], "Unknown tri backend: typo"];
expectError["unknown method even with old tri", loadWilson[{{"tridiag_method", "typo"}}], "Unknown tridiag_method backend: typo"];
expectError["unknown method with cpp", loadWilson[{{"tri", "cpp"}, {"tridiag_method", "typo"}}], "Unknown tridiag_method backend: typo"];
POL2x2 = True;
expectError["no matrix extension", loadWilson[{{"tri", "rkpw"}}], "scalar normal-state"];
POL2x2 = False;
isSC[] := True;
expectError["no superconducting extension", loadWilson[{{"tri", "cpp"}, {"tridiag_method", "rkpw"}}], "scalar normal-state"];
isSC[] := False;

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
loadWilson[{{"bandrescale", "0"}}];
check["legacy bandrescale behavior unchanged", !RKPW && Max[Abs[Flatten[xitable[1]]]] == 0];

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
