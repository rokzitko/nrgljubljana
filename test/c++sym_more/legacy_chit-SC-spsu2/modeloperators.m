Module[{t},
  t = {};

  t = Join[t, mtGlobalOp["Qtot",    number[f[0]] +   number[d[]] - 2 ]];
  t = Join[t, mtGlobalOp["Iztot", isospinz[f[0], 0] + isospinz[d[], -1] ]]; (* !!! *)
  t = Join[t, mtGlobalOp["Ixtot", isospinx[f[0], 0] + isospinx[d[], -1] ]]; (* !!! *)

  texportable = t;
];
  
texportable  
