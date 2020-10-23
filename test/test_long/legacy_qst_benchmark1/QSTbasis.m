(* 
 ch=1 <-> Tz=+1
 ch=2 <-> Tz=0
 ch=3 <-> Tz=-1
*)

(* f indexed as +1, 0, -1, in that order *)
Module[{},
 basisops = {f[+1], f[0], f[-1]};
 makebasis[basisops];
 MyPrint["basisops override -> ", basisops];
 vak = vacuum[];
 
 degnr[{q_Integer, ss_Integer, t_Integer, i___}] := ss (2t+1);

 bz = Get["one_site_basis.dat"];
 bz = Map[ {{#[[1,1]]-3, 2*#[[1,2]]+1, #[[1,3]]}, #[[2]]}&, bz];
 bz = bz /. d[CR,t_,s_]:>f[CR,t,s];
 MyPrint["bz=", bz];
 bvc = bzop2bzvc[bz, vak];
 MyPrint["bvc=", bvc];
];
