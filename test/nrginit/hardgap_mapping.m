Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];
backend = Environment["CHAIN_TEST_BACKEND"];
If[!MemberQ[{"legacy", "rkpw"}, backend], Print["Missing explicit CHAIN_TEST_BACKEND"]; Exit[1]];
method = If[backend == "legacy", "lanczos", "rkpw"];
fullTri = If[backend == "legacy", "old", "rkpw"];
settings = {{"tri", fullTri}, {"tridiag_method", method}, {"hardgap", "true"},
  {"boundary", "0.25"}, {"mMAX", "4"}, {"nrxi", "4"}, {"xmax", "99"}};
flatE[x_] := 1/4 + 3/4 If[x < 2, 2-x+(1-2^(1-x))/Log[2], 2^(2-x)/(2 Log[2])];
table = N[Table[{x, flatE[x]/2^(2-x)}, {x, 1, 11/2, 1/8}], 80];
Block[{BAND = "asymode", DY = False, DC = False, DZ = True, z = 0.5,
    MyImport, intrho, intrhoneg, theta0Ch, negativeTable},
  MyImport["Delta.dat", "Table"] := Join[
    Table[{bandrescale x, 1/(2 bandrescale)}, {x, {-1, -1/2, -1/4}}],
    Table[{bandrescale x, 1/(2 bandrescale)}, {x, {1/4, 1/2, 1}}]];
  MyImport["../FSOL.dat", "Table"] := table;
  negativeTable := table;
  MyImport["../FSOLNEG.dat", "Table"] := negativeTable;
  Do[
    loadWilson[Append[settings, {"bandrescale", ToString[width]}]];
    theta = 3/4 (1-2^(-1/2-4));
    check["retained hard-gap theta", Abs[thetaCh[1] - theta] < 10^-20];
    Do[
      hi = If[m == 0, 1, 1/4 + 3/4 2^(1/2-m)];
      lo = 1/4 + 3/4 2^(-1/2-m);
      check["hard-gap shell mass", Abs[df[1, m]/((hi-lo)/2)-1] < 10^-20];
      check["hard-gap representative", Abs[de[1, m]-flatE[3/2+m]]/(hi-lo) < 2 10^-14 && lo < de[1, m] < hi];
      check["hard-gap normalized amplitude", Abs[du[1][0, m]^2-((hi-lo)/2)/theta] < 10^-20],
      {m, 0, 4}];
    moment2 = Sum[flatE[3/2+m]^2 (If[m == 0, 1, 1/4+3/4 2^(1/2-m)]-(1/4+3/4 2^(-1/2-m))), {m, 0, 4}]/theta;
    check["physical first hopping", Abs[xitable[1][[1, 1]]/width-Sqrt[moment2]] < 10^-13],
    {width, {1, 2}}];
  Block[{z = 0.37},
    loadWilson[settings];
    interpolation = Interpolation[Map[{#[[1]], #[[2]] 2^(2-#[[1]])} &, table], InterpolationOrder -> 1];
    check["interpolate energies for off-grid twist", Abs[de[1, 4]-interpolation[1+4+z]] < 10^-20]];
  Block[{table = Most[table]},
    expectError["actual table extent, not declared xmax", loadWilson[settings], "does not cover"]];
  Block[{negativeTable = Most[table]},
    expectError["negative table coverage checked independently", loadWilson[settings], "FSOLNEG.dat"]];
  Block[{table = Reverse[table]}, expectError["unordered table", loadWilson[settings], "increasing abscissas"]];
  Block[{table = Map[{#[[1]], 0.01} &, table]}, expectError["wrong gap energies", loadWilson[settings], "outside its hard-gap shell"]];
  expectError["unsupported adaptive gap", loadWilson[Append[settings, {"adapt", "true"}]], "adapt=false"];
];
Block[{BAND = "flat"}, expectError["unsupported initializer gap band", loadWilson[settings], "band=asymode"]];
Do[expectError["invalid boundary " <> value,
  loadWilson[{{"tri", fullTri}, {"hardgap", "true"}, {"boundary", value}}], "boundary must be"],
  {value, {"-0.1", "1", "!Infinity", "!I"}}];
loadWilson[{{"tri", fullTri}, {"hardgap", "true"}, {"boundary", "0"}}];
check["zero boundary remains ungapped", Abs[xi[1][0] - (3/4)/Sqrt[7/4]] < 10^-13];
Print["hardgap ", backend, " failures: ", failures];
If[failures == 0, True, $Failed]
