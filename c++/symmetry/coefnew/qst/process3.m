fnmma2 = fn <> ".mma2";
l=Get[fnmma2];
Print[fn];
Print[Length[l]];

(* Output *)

toc[koef_] := ToString[CForm[koef]];
tos[koef_] := ToString[koef];

outmake[{a_, b_, c_}] :=
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-S)+ss1] <> ", " <> toc[c /. T -> t1] <> ")";

INVALIDVAL = 666;

string[t_] := Module[{x = t},
  x = x /. Indeterminate -> INVALIDVAL;
  x = toc[x];
  x = StringReplace[x, "Sqrt("->"sqrt("];
  x
];

proc[{i1_,ip_,IN1_,INp_,t_}] := "{ " <> tos[i1] <> ", " <> tos[ip] <> ", " <>
          outmake[IN1] <> ", " <> outmake[INp] <> ", " <> 
          string[t] <> " },";

(* Bug trap *)
proc[{i1_, ip_, IN1_, INp_, arg___}] := Module[{}, Print["ERROR! ", {i1,ip}, EXPR[arg]]; Exit[1] ];

ll = Map[proc, l];

Print[Export[fn, ll]];
