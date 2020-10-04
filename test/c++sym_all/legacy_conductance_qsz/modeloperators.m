Module[{t},
  t = {};

  t = Join[t, mtGlobalOp["SZtot", spinz[f[0]] + spinz[d[]] ] ];
  t = Join[t, mtGlobalOp["Qtot",  (number[f[0]]-1/2) + (number[d[]]-1) ] ];

  texportable = t;
];
  
texportable  
