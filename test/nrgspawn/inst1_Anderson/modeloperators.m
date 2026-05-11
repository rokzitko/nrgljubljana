Module[{t},
  t = {};
  MPVCFAST = False;
  t = Join[t, mtSingletOp["SigmaHd",   SigmaHdAvg ]];
  t = Join[t, mtSingletOp["SigmaHd-u", SigmaHd[UP] ]];
  t = Join[t, mtSingletOp["SigmaHd-d", SigmaHd[DO] ]];
  MPVCFAST = True;
  texportable = t;
];
texportable
