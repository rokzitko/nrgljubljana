Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];

loadWilson[{{"tri", "old"}, {"tridiag_method", "lanczos"}}];
check["explicit legacy selection", TRI == "old" && TRIDIAGMETHOD == "lanczos" && PREC == 1000 && !RKPW];
analytic = Table[(1 + 1/2)/2 2^(-n/2) (1 - 2^(-n-1))/Sqrt[(1 - 2^(-2n-1)) (1 - 2^(-2n-3))], {n, 0, 4}];
check["legacy flat-band analytic oracle", close[Flatten[xitable[1]], analytic] && close[Flatten[zetatable[1]], ConstantArray[0, 5]]];
Do[
  loadWilson[{{"tri", mode}, {"tridiag_method", "lanczos"}}];
  check[mode <> " explicit legacy seed", !RKPW && dothelanczos === dothelanczosold && PREC == 30 && !MachineNumberQ[xi[1][0]]],
  {mode, {"cpp", "none"}}];
loadWilson[{{"tri", "old"}, {"tridiag_method", "lanczos"}, {"bandrescale", "0"}}];
check["legacy bandrescale behavior unchanged", !RKPW && Max[Abs[Flatten[xitable[1]]]] == 0];

Print["legacy chain failures: ", failures];
If[failures == 0, True, $Failed]
