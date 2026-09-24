Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];

loadWilson[{{"tri", "old"}, {"tridiag_method", "lanczos"}}];
oldCoefficients = {Flatten[zetatable[1]], Flatten[xitable[1]]};
loadWilson[{{"tri", "rkpw"}, {"tridiag_method", "rkpw"}}];
check["rkpw agrees with short high-precision legacy chain", close[{Flatten[zetatable[1]], Flatten[xitable[1]]}, oldCoefficients]];

Print["chain crosscheck failures: ", failures];
If[failures == 0, True, $Failed]
