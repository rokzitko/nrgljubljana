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

Print["chain default failures: ", failures];
If[failures == 0, True, $Failed]
