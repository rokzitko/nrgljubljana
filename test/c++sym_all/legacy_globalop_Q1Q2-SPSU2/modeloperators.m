Module[{t},
  t = {};

  t = Join[t, mtGlobalOp["Qdiff", number[f[0]] - number[f[1]] ]];
  t = Join[t, mtGlobalOp["Qtot",  number[f[0]] + number[f[1]] + number[d[]] - 3]];
  t = Join[t, mtGlobalOp["Q1",  number[f[0]] ]];
  t = Join[t, mtGlobalOp["Q2",  number[f[1]] ]];

  texportable = t;
];
  
texportable  
