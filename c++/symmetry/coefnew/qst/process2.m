(* V1 *)

fntmp = fn <> ".tmp";
fnmma2 = fn <> ".mma2";
l=ReadList[fntmp];
Print[fn];
Print[Length[l]];

SetOptions[$Output, PageWidth -> 240];

(* simpl[] applied only on the factor *)
one[{a__, x_}] := Module[{},
  Print[a];
  l = Simplify[x, S > 5 && T > 5];
  Print["step1=", l];
  l = FullSimplify[x, S > 5 && T > 5, TimeConstraint -> 1];
  Print["step2=", l];
  {a, l}
];

l = Map[one, l];

Print[Put[l, fnmma2]];
