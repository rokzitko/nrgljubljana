Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];

(* Inspect dispatch only: even a legacy default must not execute a producer
   when this neutral test runs with TEST_CHAIN_LEGACY=OFF. *)
hookfile["hook_pre_lanczosinit"] := Throw[
  {TRI, TRIDIAGMETHOD, PREC, RKPW, DISCNMAX, dothelanczos}, "dispatch-only"];
Do[
  observed = Catch[loadWilson[pairs], "dispatch-only"];
  reported = {{"tri", TRI}, {"tridiag_method", TRIDIAGMETHOD}};
  Print["Reported initializer selection: ", reported];
  explicit = Catch[loadWilson[reported], "dispatch-only"];
  check["omitted selection equals explicit reported default", ListQ[observed] && observed === explicit];
  check["no coefficient tables generated", !ValueQ[xitable[1]] && !ValueQ[zetatable[1]]],
  {pairs, {{}, {{"tri", "cpp"}}, {{"tri", "none"}}}}];

Do[
  method = {"lanczos", "rkpw"}[[methodIndex]];
  pairs = {{"tri", spec[[1]]}, {"tridiag_method", method}};
  producer = spec[[methodIndex + 2]];
  expected = {spec[[2]], producer === dothelanczosrkpw, producer};
  label = spec[[1]] <> "/" <> method;
  observed = Catch[loadWilson[pairs], "dispatch-only"];
  check[label <> " precision and producer", observed[[{3, 4, 6}]] === expected];
  observed = Catch[loadWilson[Append[pairs, {"prec", "73"}]], "dispatch-only"];
  check[label <> " explicit precision", observed[[{3, 4, 6}]] === ReplacePart[expected, 1 -> 73]];
  Block[{option},
    option["GENERATE_TEMPLATE"] = True;
    observed = Catch[loadWilson[pairs], "dispatch-only"];
    check[label <> " template override", observed[[{3, 4, 6}]] === ReplacePart[expected, 3 -> None]];
  ],
  {spec, {
    {"old", 1000, dothelanczosold, dothelanczosold},
    {"rkpw", 30, dothelanczosrkpw, dothelanczosrkpw},
    {"cpp", 30, dothelanczosold, dothelanczosrkpw},
    {"none", 30, dothelanczosold, dothelanczosrkpw},
    {"manual", 30, loaddiscretizationtables, loaddiscretizationtables},
    {"manual_nambu", 30, loaddiscretizationtables, loaddiscretizationtables},
    {"manual_nambu_new", 30, loaddiscretizationtables, loaddiscretizationtables}}},
  {methodIndex, 2}];

Print["chain default failures: ", failures];
If[failures == 0, True, $Failed]
