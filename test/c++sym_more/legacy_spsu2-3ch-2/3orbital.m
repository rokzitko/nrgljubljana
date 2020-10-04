(* V2 *)

(* If[MODEL == "3ORBITAL" && VARIANT == "SIMPLE1", *)

 Module[{SPIN, ORB, id1, id2, sx, sy, sz, tx, ty, tz, Sx, Sy, Sz, Tx, Ty, Tz, ss, tt},

  def3ch[0]; (* 3 channels, no fermionic impurities *)

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
  Sx = spinx[f[0]] + spinx[f[1]] + spinx[f[2]];
  Sy = spiny[f[0]] + spiny[f[1]] + spiny[f[2]];
  Sz = spinz[f[0]] + spinz[f[1]] + spinz[f[2]];
  
  crops[sigma_] = { f[CR,0,sigma], f[CR,1,sigma], f[CR,2,sigma] };
  anops[sigma_] = { f[AN,0,sigma], f[AN,1,sigma], f[AN,2,sigma] };

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

(*   H = H /. {bra[i_, _] :> bra[i], ket[i_, _] :> ket[i]}; *)

  Print[H];

  MAKESPINKET = SPIN;
  MAKEORBKET = ORB;
  (* Note: need to add the orbital part manually! *)

 ];

(* ]; *)
