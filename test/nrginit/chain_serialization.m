(* Full model generation and actual makedata/maketable serialization.
   Each invocation uses a fresh kernel and a private working directory. *)
sourceDir = Environment["CHAIN_SERIALIZATION_SOURCE"];
backend = Environment["CHAIN_TEST_BACKEND"];
mode = Environment["CHAIN_SERIALIZATION_MODE"];
exponent = ToExpression[Environment["CHAIN_SERIALIZATION_EXPONENT"]];
If[!MemberQ[{"legacy", "rkpw"}, backend] ||
    !MemberQ[{"full", "cpp", "none"}, mode] || !MemberQ[{20, 305}, exponent],
  Print["Invalid chain serialization test selection"]; Exit[1]];
SetDirectory[Environment["CHAIN_SERIALIZATION_WORK"]];

tri = If[mode == "full", If[backend == "rkpw", "rkpw", "old"], mode];
(* Full tri selection takes precedence over the runtime method selector. *)
method = If[mode == "full", If[backend == "rkpw", "lanczos", "rkpw"],
  If[backend == "rkpw", "rkpw", "lanczos"]];
clipOption = If[exponent == 20, "CHOP", "EPSCLIP"];
threshold = If[exponent == 20, 10^-12, 10^-300];
parameterText = StringRiffle[{
  "[param]", "symtype=QS", "model=CLEAN", "Ninit=1", "Nmax=3", "nrxi=3",
  "discretization=Y", "Lambda=2", "prec=80", "mmadebug=0",
  "tri=" <> tri, "tridiag_method=" <> method,
  "bandrescale=1e-" <> ToString[exponent], "shift0=0.125",
  "mMAX=" <> If[exponent == 20, "80", "998"]}, "\n"] <> "\n";
Export["param", parameterText, "Text"];

(* Use only source-tree packages, not a developer's installed initializer. *)
PACKAGEPATH = {FileNameJoin[{sourceDir, "nrginit"}]};
SYMTYPE = "runtime";
Quiet[Get["sneg.m", Path -> PACKAGEPATH], General::shdw];
Get["initial.m", Path -> PACKAGEPATH];
failures = 0;
check[label_, condition_] := If[!TrueQ[condition], Print["FAILED: ", label]; failures++];
check["selected reconstruction", RKPW === (backend == "rkpw")];

(* Parse only the numeric suffix blocks of the emitted file. Use arbitrary
   precision for comparisons so squaring or subtracting tiny doubles cannot
   make a broken serialization check pass through underflow. *)
readBlock[file_, marker_] := Module[{lines, start, tokens, count, arrays = {}},
  lines = StringTrim /@ Import[file, "Lines"];
  start = FirstPosition[lines, marker];
  If[MissingQ[start], Return[{}]];
  tokens = StringSplit[StringRiffle[Drop[lines, First[start]], " "]];
  Do[
    count = ToExpression[First[tokens]] + 1;
    AppendTo[arrays, ToExpression[StringReplace[#, {"e" -> "*^", "E" -> "*^"}]] & /@
      Take[Rest[tokens], count]];
    tokens = Drop[tokens, count + 1],
    {If[marker == "z", 2, 4]}];
  arrays
];
sameNumbers[actual_, expected_] := Dimensions[actual] === Dimensions[expected] &&
  And @@ MapThread[If[#2 == 0, #1 == 0,
    #1 != 0 && Abs[SetPrecision[#1, 50]/SetPrecision[#2, 50] - 1] < 10^-14] &,
    {Flatten[actual], Flatten[expected]}];

chain = {Flatten[xitable[1]], Flatten[zetatable[1]]};
star = Flatten /@ {eptable[1], emtable[1], u0ptable[1], u0mtable[1]};
check["tiny representable physical hopping", 0 < chain[[1, 1]] < threshold];
check["tiny physical onsite shift", 0 < chain[[2, 1]] < threshold &&
  Abs[SetPrecision[chain[[2, 1]], 50]/(10^-exponent/8) - 1] < 10^-14];
check["tiny representable star pole", 0 < star[[1, -1]] < threshold &&
  MachineNumberQ[N[star[[1, -1]]]] && N[star[[1, -1]]] != 0.];

makedata["unclipped"];
check["unclipped data written", FileExistsQ["unclipped"]];
check["tiny physical seed energy", 0 < Abs[GSenergy] < threshold &&
  !StringContainsQ[Import["unclipped", "Text"], "\ne\n0\n"]];
If[mode != "none",
  check["serialized physical chain preserved", sameNumbers[readBlock["unclipped", "z"], chain]],
  check["none omits chain", readBlock["unclipped", "z"] === {}]];
If[mode == "cpp",
  check["serialized C++ star preserved", sameNumbers[readBlock["unclipped", "T"], star]]];

If[backend == "rkpw",
  (* Exercise the public serialization path, including options added after
     initialization. A rejected write must not truncate a previous data file. *)
  Do[
    target = "rejected-" <> StringRiffle[options, "-"];
    result = Block[{loptions = options, MyError, perturbhamiltonian},
      MyError[args__] := Throw[StringJoin[ToString /@ {args}], "serialization-error"];
      perturbhamiltonian[] := Throw["table assembly reached", "serialization-error"];
      Catch[makedata[target], "serialization-error"]
    ];
    check["reject " <> target, StringQ[result] &&
      StringContainsQ[result, "rkpw serialization is incompatible with CHOP/EPSCLIP"]];
    check["no rejected output " <> target, !FileExistsQ[target]];
    Export[target, "previous data\n", "Text"];
    previousHash = FileHash[target];
    Block[{loptions = options, MyError},
      MyError[args__] := Throw["rejected", "serialization-error"];
      Catch[makedata[target], "serialization-error"]
    ];
    check["existing output unchanged " <> target, FileHash[target] === previousHash],
    {options, {{"CHOP"}, {"EPSCLIP"}, {"CHOP", "EPSCLIP"}}}],

  Block[{loptions = {clipOption}}, makedata["clipped"]];
  check["legacy seed energy clipped", StringContainsQ[Import["clipped", "Text"], "\ne\n0\n"]];
  If[mode != "none",
    clippedChain = readBlock["clipped", "z"];
    check["legacy physical hopping clipped", clippedChain[[1, 1]] == 0];
    check["legacy physical onsite clipped", clippedChain[[2, 1]] == 0]];
  If[mode == "cpp",
    clippedStar = readBlock["clipped", "T"];
    expectedStar = star /. value_Real /; Abs[value] < threshold -> 0;
    check["legacy C++ star clipping unchanged", sameNumbers[clippedStar, expectedStar]]];
];

Print["chain serialization ", backend, " ", mode, " ", clipOption, " failures: ", failures];
If[failures == 0, True, $Failed]
