(*
  greens.m - Green's funciton on the Wilson chain
  rok.zitko@ijs.si, May 2008
*)

(* Lambda: NRG discretization parameter *)
If[!ValueQ[Lambda], Lambda = 2.];

(* Dimension of the sparse matrix *)
If[!ValueQ[M], M = 40];

MyImport[fn_, type___] := Module[{l},
  l = Import[fn, type];
  If[l === $Failed,
    Print["Import failed ", fn];
    Exit[];
  ];
  l
];

xilist = Import["xi", "Table"];
xi[n_Integer] := SetPrecision[xilist[[n+1, 1]], PREC] /; n >= 0;

(* Print[xilist]; *)

(* Get the Green's function on the first site *)
green[] := Module[{s},
  Print["green[]"];
  s = Table[
    Switch[i - j,  1, xi[j-1] Lambda^((M - i)/2),
                  -1, xi[i-1] Lambda^((M - j)/2),
                   _, 0 ],
                   {i, 1, M}, {j, 1, M}];

  v1 = gammaPol (1/omegaN) { Table[If[i==1, 1, 0], {i, 1, M}] };
  m0 = {{ 0 }};                 

  Print[v1];
  s2 = ArrayFlatten[{{ m0, v1 }, {Transpose[v1], s}}];
  
  {vals, vecs} = Transpose[Sort[Transpose[Eigensystem[s2]]]];
  site = vecs[[All,1]]^2;
  Print["sum first=", Plus @@ site];
  ( (omega-vals)^-1 ) . site
];

(* Extract the scattering phase using the Green function *)
(* 'grf' must be defined! *)
phifromgf[x_?NumberQ] := Module[{kappa, delta},
  kappa = (2/Pi) Lambda^((M - 1)/2) grf /. omega -> x;
  delta = -Arg[I + kappa];
  val = Floor[delta/Pi];
  (* Result in units of Pi *)
  1/Pi (delta - Pi val)
];

phifromgf[l_List] := Module[{l2, ndx, phis},
  l2 = Abs[l]; (* absolute values *)
  ndx = Position[l2, Min[l2]] [[1,1]]; (* Find the index *)
  phis = Map[phifromgf, l]; (* Calculate all phase shifts *)
  phis[[ndx]]
];
