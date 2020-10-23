params = Join[params, {Bx -> extraBx, By -> extraBy, Bz -> extraBz, g -> extrag, kxz -> extrakxz, kyz -> extrakyz}];

def1ch[0];

If[!paramexists["spin", "extra"],
  MyError["Define the spin of the impurity!"];
];
SPIN = ToExpression @ param["spin", "extra"];
MyPrint["SPIN=", SPIN];

Module[{sz, sp, sm, sx, sy, oz, op, om, ss},
  sz = spinketbraZ[SPIN];
  sp = spinketbraP[SPIN];
  sm = spinketbraM[SPIN];
  sx = spinketbraX[SPIN];
  sy = spinketbraY[SPIN];

  oz = nc[ sz, spinz[ f[0] ] ];
  op = nc[ sp, spinminus[ f[0] ] ];
  om = nc[ sm, spinplus[ f[0] ] ];
  ss = oz + 1/2 (op + om) // Expand;
  MyVPrint[2, "KONDO ss=", ss];

  sxz = Simplify[nc[sx, sz] + nc[sz, sx]];
  syz = Simplify[nc[sy, sz] + nc[sz, sy]];

  Print2["sxz=", sxz];
  Print2["syz=", syz];

  (* Warning: sy is imaginary! *)
  Hc = Jkondo ss + g ( Bx sx + By sy + Bz sz ) + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]) +
       kxz sxz + kyz syz;
];

H = H0 + Hc;

MAKESPINKET = SPIN; (* See QSZaddspinket[] function in SNEG *)
