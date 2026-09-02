MyPrint["modelops begin"];
Module[{t},
  t = {};
  MPVCFAST = False;
  t = Join[t, mtSingletOp["NPH1",   phononnumber[1, cutoffs] ]];
  t = Join[t, mtSingletOp["DISPL1", phononx[1, cutoffs]      ]];
  t = Join[t, mtSingletOp["NPH2",   phononnumber[2, cutoffs] ]];
  t = Join[t, mtSingletOp["DISPL2", phononx[2, cutoffs]      ]];
  MPVCFAST = True;
  texportable = t;
];
MyPrint["modelops end"];
texportable
