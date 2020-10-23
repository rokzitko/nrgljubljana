Module[{t},
 t = {};
 (* \vc{S} \cdot \vc{n} with \vc{n}=(1,1,1)/sqrt{3} *)
 t = Join[t, mtSingletOp["SXYZd", (spinx[d[]] + spiny[d[]] + spinz[d[]])/Sqrt[3] ]];
 
 texportable = t;
];

texportable
