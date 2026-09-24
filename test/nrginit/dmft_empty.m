Get[FileNameJoin[{sourceDir, "test", "nrginit", "chain_test_setup.m"}]];
Get[FileNameJoin[{sourceDir, "nrginit", "sym.m"}]];
backend = Environment["CHAIN_TEST_BACKEND"];
If[!MemberQ[{"legacy", "rkpw"}, backend], Print["Missing explicit CHAIN_TEST_BACKEND"]; Exit[1]];
method = If[backend == "legacy", "lanczos", "rkpw"];
fullTri = If[backend == "legacy", "old", "rkpw"];

(* Supply the [dmft] block when loadWilson dispatches the real constructor.
   No star coefficients or integration routines are substituted. *)
loadDMFT[density_, pairs_:{}] := Block[{BAND = "dmft", loadmodule},
  loadmodule["dmft.m", True] := (
    listkeywords["dmft"] = {"gamma"};
    data["dmft"]["gamma"] = density;
    Get[FileNameJoin[{sourceDir, "nrginit", "dmft.m"}]];
  );
  loadWilson[Join[{{"tri", fullTri}, {"tridiag_method", method}, {"mMAX", "5"},
    {"nrxi", If[backend == "rkpw", "3", "1"]}}, pairs]];
];

(* Four nonempty shells, asymmetric in both energy and mass. All other shells
   are truly empty, including internal gaps and the low-energy tail. *)
densityExpression = Piecewise[{{2 eps, 1/2 < eps <= 1}, {3 eps, 1/16 < eps <= 1/8},
  {-4 eps, -1/2 <= eps < -1/4}, {-5 eps, -1/4 <= eps < -1/8}}];
positiveMass = {3/4, 0, 0, 9/512, 0, 0};
negativeMass = {0, 3/8, 15/128, 0, 0, 0};
theta = 645/512;
weights = {3/4, 9/512, 3/8, 15/128}/theta;

Do[
  DY = scheme == "Y"; DC = scheme == "C"; DZ = False;
  density = ToString[scale densityExpression, InputForm];
  loadDMFT[density];
  label = backend <> " " <> scheme <> " scale=" <> ToString[N[scale], InputForm];
  check[label <> " explicit selectors", TRI == fullTri && TRIDIAGMETHOD == method];
  check[label <> " theta", Abs[thetaCh[1]/(scale theta) - 1] < 10^-20];
  check[label <> " analytic shell masses", close[
    {Table[df[1, m], {m, 0, mMAX}], Table[dfminus[1, m], {m, 0, mMAX}]}/scale,
    {positiveMass, negativeMass}, 10^-20]];
  check[label <> " fixed star table layout",
    Dimensions /@ {eptable[1], emtable[1], u0ptable[1], u0mtable[1]} === ConstantArray[{6, 1}, 4]];
  Do[
    masses = If[side == 1, positiveMass, negativeMass];
    Do[
      mass = If[side == 1, df[1, m], dfminus[1, m]];
      energy = If[side == 1, de[1, m], deminus[1, m]];
      amplitude = If[side == 1, du[1][0, m], dv[1][0, m]];
      check[label <> " finite in-shell energy " <> ToString[{side, m}],
        NumberQ[energy] && Element[energy, Reals] && km[m+1] <= energy <= km[m]];
      If[masses[[m+1]] == 0,
        check[label <> " exact zero mass/amplitude " <> ToString[{side, m}],
          NumberQ[mass] && Sign[mass] === 0 && NumberQ[amplitude] && Sign[amplitude] === 0];
        check[label <> " inert midpoint " <> ToString[{side, m}],
          Abs[energy/((km[m+1] + km[m])/2) - 1] < 10^-20],
        check[label <> " positive mass/amplitude " <> ToString[{side, m}], mass > 0 && amplitude > 0]
      ], {m, 0, mMAX}], {side, 2}];

  energies = If[DY, 7/9, 3/4] {1, 1/8, -1/2, -1/4};
  check[label <> " nonempty representative energies",
    close[{de[1, 0], de[1, 3], -deminus[1, 1], -deminus[1, 2]}, energies, 10^-20]];
  check[label <> " normalized star", Abs[Sum[du[1][0, m]^2 + dv[1][0, m]^2, {m, 0, mMAX}] - 1] < 10^-20];
  Do[
    expectedMoment = weights . energies^order;
    starMoment = Sum[de[1, m]^order du[1][0, m]^2 + (-deminus[1, m])^order dv[1][0, m]^2, {m, 0, mMAX}];
    check[label <> " star moment " <> ToString[order], Abs[starMoment - expectedMoment] < 10^-20],
    {order, 0, 7}];
  mean = weights . energies;
  check[label <> " first onsite and hopping",
    close[{dzeta[1][0], xi[1][0]}, {mean, Sqrt[weights . energies^2 - mean^2]}]];
  onsite = Flatten[zetatable[1]]; hopping = Flatten[xitable[1]];
  jacobi = DiagonalMatrix[onsite] + DiagonalMatrix[Most[hopping], 1] + DiagonalMatrix[Most[hopping], -1];
  Do[check[label <> " chain moment " <> ToString[order],
    Abs[MatrixPower[jacobi, order][[1, 1]] - weights . energies^order] < 10^-13],
    {order, 0, 2 Length[onsite] - 1}];
  If[backend == "rkpw",
    check[label <> " finite support terminal hopping", Length[onsite] == 4 && Last[hopping] === 0. && Min[Most[hopping]] > 0];
    expectError[label <> " support count", loadDMFT[density, {{"nrxi", "4"}}], "support is 4"],
    (* Legacy has no exact-support termination contract: request only a prefix. *)
    check[label <> " legacy short prefix", Length[onsite] == 2 && Min[hopping] > 0]
  ],
  {scheme, {"Y", "C"}}, {scale, {1, 10^-80}}];

expectError["negative shell mass remains invalid", loadDMFT[ToString[-densityExpression, InputForm]], "Coefficient is negative"];
expectError["non-real density remains invalid", loadDMFT["I"], "Check definition for gamma"];
expectError["nonnumeric density remains invalid", loadDMFT["unknownDensity"], "Check definition for gamma"];

Print["dmft empty ", backend, " failures: ", failures];
If[failures == 0, True, $Failed]
