Module[{t},
  t = {};

  t = Join[t, mtGlobalOp["SZtot", spinz[f[0]] + spinketbraZ[SPIN]]];
  t = Join[t, mtGlobalOp["SYtot", spiny[f[0]] + spinketbraY[SPIN]]];
  t = Join[t, mtGlobalOp["SXtot", spinx[f[0]] + spinketbraX[SPIN]]];

  texportable = t;
];
  
texportable  
