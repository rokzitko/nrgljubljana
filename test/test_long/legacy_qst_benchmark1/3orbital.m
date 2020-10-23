 Module[{SPIN, ORB, id1, id2, sx, sy, sz, tx, ty, tz, Sx, Sy, Sz, Tx, Ty, Tz, ss, tt},

  def3ch[0]; (* 3 channels, no fermionic impurities *)
  
  H0 = H0 /. {f[t_, 0, s_] :> f[t, -1, s],
              f[t_, 1, s_] :> f[t,  0, s],
              f[t_, 2, s_] :> f[t, +1, s]};

  SPIN = ToExpression @ param["spin", "extra"]; (* Should be SPIN=1 *)
  ORB  = ToExpression @ param["orb",  "extra"]; (* Should be ORB=1 *)
  MyPrint["SPIN=", SPIN, " ORB=", ORB];
  
  (* Identity operators *)
  id1 = spinketbraI[SPIN, {1,0}];
  id2 = spinketbraI[ORB,  {0,1}];
  
  (* Impurity spin operators *)
  sx = spinketbraX[SPIN, {1,0}] ~ nc ~ id2;
  sy = spinketbraY[SPIN, {1,0}] ~ nc ~ id2;
  sz = spinketbraZ[SPIN, {1,0}] ~ nc ~ id2;

  (* Impurity orbital operators *)  
  tx = spinketbraX[ORB, {0,1}] ~ nc ~ id1;
  ty = spinketbraY[ORB, {0,1}] ~ nc ~ id1;
  tz = spinketbraZ[ORB, {0,1}] ~ nc ~ id1;  
  
  (* Band spin operators *)
  Sx = spinx[f[-1]] + spinx[f[0]] + spinx[f[+1]];
  Sy = spiny[f[-1]] + spiny[f[0]] + spiny[f[+1]];
  Sz = spinz[f[-1]] + spinz[f[0]] + spinz[f[+1]];
  
  crops[sigma_] = { f[CR,1,sigma], f[CR,0,sigma], f[CR,-1,sigma] };
  anops[sigma_] = { f[AN,1,sigma], f[AN,0,sigma], f[AN,-1,sigma] };

  (* Band orbital operators *)
  Tx = Sum[VMV[crops[sigma], spinmatrixX[1], anops[sigma]], {sigma, 0, 1}];
  Ty = Sum[VMV[crops[sigma], spinmatrixY[1], anops[sigma]], {sigma, 0, 1}];
  Tz = Sum[VMV[crops[sigma], spinmatrixZ[1], anops[sigma]], {sigma, 0, 1}];

  (* Exchange operators *)
  ss = nc[sx, Sx] + nc[sy, Sy] + nc[sz, Sz];
  tt = nc[tx, Tx] + nc[ty, Ty] + nc[tz, Tz];
  
  ss = Expand[ss];
  tt = Expand[tt];
  
  (* Coupling Hamiltonian *)
  Hc = Jkondo1 ss + Jkondo2 tt;
  
  H = H0 + Hc;

  Print[H];

  MAKESPINKET = SPIN;
  MAKEORBKET = ORB;
];

