Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];
Get[FileNameJoin[{sourceDir, "nrginit", "sym.m"}]];
backend = Environment["CHAIN_TEST_BACKEND"];
If[!MemberQ[{"legacy", "rkpw"}, backend], Print["Missing explicit CHAIN_TEST_BACKEND"]; Exit[1]];
method = If[backend == "legacy", "lanczos", "rkpw"];
fullTri = If[backend == "legacy", "old", "rkpw"];

(* Piecewise-constant one-sided densities with analytic geometric shell sums.
   Mock only file import: the actual asymode interpolants, integrals, star
   normalization and selected reconstruction all execute. *)
Block[{BAND = "asymode", DY = False, DC = False, DZ = True, COEFCHANNELS = 2, CHANNELS = 2,
    MyImport, intrho, intrhoneg, theta0Ch},
  positive = {1, 3}; negative = {2, 1};
  energyFactor = 1/(2 Log[2]);
  Do[With[{file = "Delta.dat" <> If[ch == 1, "", ToString[ch]], channel = ch},
    MyImport[file, "Table"] :=
      Join[Table[{bandrescale x, negative[[channel]]}, {x, {-1, -1/2, -1/4}}],
           Table[{bandrescale x, positive[[channel]]}, {x, {1/4, 1/2, 1}}]]
  ], {ch, 2}];
  Do[With[{file = "../" <> table <> If[ch == 1, "", ToString[ch]]},
    MyImport[file, "Table"] = {{1, N[energyFactor, 80]}, {30, N[energyFactor, 80]}}
  ], {ch, 2}, {table, {"FSOL.dat", "FSOLNEG.dat"}}];
  Do[
    loadWilson[{{"tri", fullTri}, {"tridiag_method", method}, {"hardgap", "false"},
      {"mMAX", "80"}, {"nrxi", "4"}, {"bandrescale", ToString[width]}}];
    Do[
      check["asymmetric theta " <> ToString[{width, ch}],
        Abs[thetaCh[ch] - width (positive[[ch]] + negative[[ch]])] < 10^-20];
      Do[
        shellWidth = 2^(-m-1);
        check["positive shell weight " <> ToString[{width, ch, m}],
          Abs[df[ch, m]/(width positive[[ch]] shellWidth) - 1] < 10^-20];
        check["negative shell weight " <> ToString[{width, ch, m}],
          Abs[dfminus[ch, m]/(width negative[[ch]] shellWidth) - 1] < 10^-20],
        {m, {0, 1, 5, 30, 60}}];
      mass = 1 - 2^(-mMAX-1);
      mean = (positive[[ch]] - negative[[ch]])/(positive[[ch]] + negative[[ch]]) *
        energyFactor (1/2) (1 - (1/2)^(2 (mMAX+1)))/((1 - 1/4) mass);
      secondMoment = energyFactor^2 (1/2) (1 - (1/2)^(3 (mMAX+1)))/((1 - 1/8) mass);
      check["normalized asymmetric mean " <> ToString[{width, ch}], Abs[dzeta[ch][0] - mean] < 10^-13];
      check["asymmetric first hopping " <> ToString[{width, ch}], Abs[xi[ch][0] - Sqrt[secondMoment - mean^2]] < 10^-13];
      check["physical onsite scale " <> ToString[{width, ch}], Abs[zetatable[ch][[1, 1]] - width mean] < 10^-13];
      check["normalized star " <> ToString[{width, ch}],
        Abs[Sum[du[ch][0, m]^2 + dv[ch][0, m]^2, {m, 0, mMAX}] - 1] < 10^-20],
      {ch, 2}],
    {width, {1, 2}}];
];

(* Stop before any reconstruction when checking guard activation. *)
dispatch[pairs_] := Block[{hookfile},
  hookfile["hook_pre_lanczosinit"] := Throw["ready", "dispatch"];
  Catch[loadWilson[Join[{{"tri", "cpp"}, {"tridiag_method", method}}, pairs]], "dispatch"]
];
Do[
  Block[{Ninit = seed},
    Do[expectError["deferred gap " <> ToString[{seed, width, value}],
      dispatch[{{"gap", value}, {"bandrescale", ToString[width]}}], "all-site onsite corrections: gap"],
      {width, {1, 2}}, {value, {"0.125", "-0.125", "1e-300", "1e-999", "!0", "!1-1", "!Infinity", "!I", "nonsense"}}];
    Do[check["zero deferred correction " <> value, dispatch[{{"gap", value}}] === "ready"],
      {value, {"0", "-0", "+0.0", ".0", "0.", "0e-999", "00.00E+999"}}];
    check["absent deferred correction", dispatch[{}] === "ready"];
    loadWilson[{{"tri", "cpp"}, {"tridiag_method", method}, {"bandrescale", "2"}, {"shift0", "0.125"}}];
    check["seed-only shift0 " <> ToString[seed], Abs[zetatable[1][[1, 1]] - 1/4] < 10^-13];
    If[seed > 0, check["shift0 does not affect later sites", Abs[zetatable[1][[2, 1]]] < 10^-13]];
  ], {seed, {0, 1}}];

Block[{POLARIZED = True, SYMTYPE, CHANNELS = 1, COEFCHANNELS = 2},
  Do[
    SYMTYPE = symmetry;
    expectError["active deferred globalh " <> symmetry, dispatch[{{"globalh", "0.125"}}], "all-site onsite corrections: globalh"];
    If[method == "rkpw",
      expectError["zero globalh retains SC handoff guard " <> symmetry, dispatch[{{"globalh", "0"}}], "pairing-table handoff"],
      check["zero globalh " <> symmetry, dispatch[{{"globalh", "0"}}] === "ready"]],
    {symmetry, {"SPU1", "P", "PP", "NONE"}}];
  SYMTYPE = "QSZ";
  check["inactive globalh is not rejected", dispatch[{{"globalh", "0.125"}}] === "ready"];
  loadWilson[{{"tri", "cpp"}, {"tridiag_method", method}, {"globalB", "0.5"}, {"bandrescale", "2"}}];
  check["globalB stays in seed", Abs[zetatable[1][[1, 1]] - 1/2] < 10^-13 && Abs[zetatable[2][[1, 1]] + 1/2] < 10^-13];
  Block[{BAND = "flat_with_bulk_field", DY = False, DZ = True},
    check["star-encoded bulkh is allowed", dispatch[{{"bulkh", "0.2"}}] === "ready"]];
];
Block[{POLARIZED = False, SYMTYPE = "SPU1"},
  If[method == "rkpw",
    expectError["inactive globalh retains SC guard", dispatch[{{"globalh", "0.125"}}], "pairing-table handoff"],
    check["unpolarized globalh inactive", dispatch[{{"globalh", "0.125"}}] === "ready"]];
];
loadWilson[{{"tri", fullTri}, {"tridiag_method", method}, {"gap", "0.125"}, {"bandrescale", "2"}}];
check["full initializer retains all-site gap", close[Flatten[zetatable[1]], Table[-(-1)^n/4, {n, 0, 4}]]];
loadWilson[{{"tri", "none"}, {"tridiag_method", method}, {"gap", "0.125"}, {"bandrescale", "2"}}];
check["external handoff remains unchanged", Abs[zetatable[1][[1, 1]] + 1/4] < 10^-13];

Print["mapping ", backend, " failures: ", failures];
If[failures == 0, True, $Failed]
