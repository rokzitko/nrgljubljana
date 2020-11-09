(* *)

PACKAGEPATH = {".", ".."};
Print["Dir: ", Directory[]];
Print["Path: ", PACKAGEPATH];
SYMTYPE = "runtime";

(* Load sneg first! *)
Get["sneg.m", Path->PACKAGEPATH];

res=Get["initialparse.m", Path->PACKAGEPATH];
If[res == $Failed,
 Print["Can't load initialparse.m"];
 Exit[1];
];

MyPrint["Parsing parameters"];
res = parse["param"];
If[res == $Failed,
 Print["Failed parsing parameter file."];
 Exit[1];
];
PARSED=True;

Get["initial.m", Path->PACKAGEPATH];

(* *)

ctrok = 0;
ctrfailed = 0;

test[a_, b_] := Module[{},
  If[a === b,
    ctrok++ ; Print["OK: a=", a, "   b=", b],
    ctrfailed++ ; Print["FAILED: a=", a, "   b=", b]
  ];
];
 
report[] := Module[{},
  Print["nr ok: ", ctrok];
  Print["nr failed: ", ctrfailed];
];

test[bazavc[{-2,1}], {vc[0,0,0,0]}]; (* no electrons *)
test[bazavc[{-1,2}], {vc[0,0,1,0],vc[1,0,0,0]}]; (* one spin-up electron *)
test[getdiffvc[{vc[0,0,1,0],vc[1,0,0,0]}], {vc[0,0,1,0],vc[1,0,0,0]}];
test[getdiffvc[{vc[0,0,1,0]+vc[1,0,0,0]}], {vc[0,0,1,0],vc[1,0,0,0]}];
test[getdiffvc[{{vc[0,0,1,0]/Sqrt[2]},{2vc[1,0,0,0]}}], {vc[0,0,1,0],vc[1,0,0,0]}];
test[transformrule[{vc[0,0,1,0],vc[1,0,0,0]}], a={vc[0, 0, 1, 0] -> {1, 0}, vc[1, 0, 0, 0] -> {0, 1}} ];
test[bazavcdiffvc[{-1,2}], a={vc[0, 0, 1, 0], vc[1, 0, 0, 0]} ];
test[bazavctransf[{-1,2}], {{1, 0}, {0, 1}} ];
test[StringJoinSep[" ", {"a", "b", "c"}], "a b c"];
mat = op2matrix[d[CR,UP], {-1,2}, {-2,1 }];
test[mat,  {{0}, {1}}];
test[Dimensions[mat], {2,1}];
test[Dimensions[mat], {Length[bazavc[{-1,2}]], Length[bazavc[{-2,1}]]}];
test[op2matrix[f[CR,0,UP], {-1,2}, {-2,1 }], {{1}, {0}}];


report[]
