def2ch[0];

SPIN = ToExpression @ param["spin", "extra"];

Module[{sx, sy, sz, ox, oy, oz, ss},
 sx = spinketbraX[SPIN];
 sy = spinketbraY[SPIN];
 sz = spinketbraZ[SPIN];

 ox1 = nc[ sx, spinx[ f[0] ] ];
 oy1 = nc[ sy, spiny[ f[0] ] ];
 oz1 = nc[ sz, spinz[ f[0] ] ];

 ox2 = nc[ sx, spinx[ f[1] ] ];
 oy2 = nc[ sy, spiny[ f[1] ] ];
 oz2 = nc[ sz, spinz[ f[1] ] ];

 ss1 = Expand[ox1 + oy1 + oz1];
 ss2 = Expand[ox2 + oy2 + oz2];
 Hk = Jkondo1 ss1 + Jkondo2 ss2;
];

H = H0 + Hk;

MAKESPINKET = SPIN;
