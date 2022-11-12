Module[{t},
  t = {};
  t = Join[t, mtDoubletOp["A_1", d[], 1 ]];
  t = Join[t, mtDoubletOp["A_2", a[], 1 ]];
  MPVCFAST = False;
  t = Join[t, mtDoubletOp["self_1", selfopd, 1 ]];
  t = Join[t, mtDoubletOp["self_2", selfopa, 1 ]];
  t = Join[t, mtDoubletOp["selfP_1", selfopdP, 1 ]];
  t = Join[t, mtDoubletOp["selfP_2", selfopaP, 1 ]];
  t = Join[t, mtSingletOp["SigmaH1u", SigmaHd[UP] ]];
  t = Join[t, mtSingletOp["SigmaH1d", SigmaHd[DO] ]];
  t = Join[t, mtSingletOp["SigmaH2u", SigmaHa[UP] ]];
  t = Join[t, mtSingletOp["SigmaH2d", SigmaHa[DO] ]];
  MPVCFAST = True;
  texportable = t;
];
texportable
