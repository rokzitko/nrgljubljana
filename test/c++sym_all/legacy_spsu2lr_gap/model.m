def2ch[0];

If[!paramexists["spin", "extra"],
  MyError["Define the spin of the impurity!"];
];
SPIN = ToExpression @ param["spin", "extra"];
MyPrint["SPIN=", SPIN];

Module[{sz, sp, sm, oz1, op1, om1, oz2, op2, om2, ss1, ss2},
 sx = spinketbraX[SPIN];
 sy = spinketbraY[SPIN];
 sz = spinketbraZ[SPIN];
 sp = spinketbraP[SPIN];
 sm = spinketbraM[SPIN];

 oz1 = nc[sz, 1/2(number[f[0], UP] - number[f[0], DO])];
 op1 = nc[sp, f[CR, 0, DO], f[AN, 0, UP]];
 om1 = nc[sm, f[CR, 0, UP], f[AN, 0, DO]];

 oz2 = nc[sz, 1/2(number[f[1], UP] - number[f[1], DO])];
 op2 = nc[sp, f[CR, 1, DO], f[AN, 1, UP]];
 om2 = nc[sm, f[CR, 1, UP], f[AN, 1, DO]];

 ss1 = oz1 + 1/2 (op1 + om1) // Expand;
 ss2 = oz2 + 1/2 (op2 + om2) // Expand;

 Hc = Jkondo1 ss1 + Jkondo2 ss2;
];
 

H = H0 + Hc;

MAKESPINKET = SPIN; (* See QSZaddspinket[] function in SNEG *)

If[isLR[],
  lrchain = {f[0], f[1]};
];
  
