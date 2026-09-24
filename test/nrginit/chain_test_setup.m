(* Standalone Wilson-chain tests: no SNEG or model generation. *)
DEBUG = 0;
Get[FileNameJoin[{sourceDir, "nrginit", "misc.m"}]];
Get[FileNameJoin[{sourceDir, "nrginit", "initialparse.m"}]];
MyError[args__] := Throw[StringJoin[ToString /@ {args}], "initializer-error"];
hook[_] := Null;
hookfile[_] := Null;
option[_] := False;
isSC[] := False;
SYMTYPE = "QS";
COEFCHANNELS = CHANNELS = 1;
POLARIZED = POL2x2 = RUNGS = False;
BAND = "flat";
DY = True; DC = DZ = False;
lambda = 2.; z = 1.; bandrescale = 1;
Ninit = 0;

failures = 0;
check[label_, condition_] := If[!TrueQ[condition], Print["FAILED: ", label]; failures++];
SetAttributes[expectError, HoldRest];
expectError[label_, expression_, fragment_] := Module[{result},
  result = Catch[expression, "initializer-error"];
  check[label, StringQ[result] && StringContainsQ[result, fragment]];
];
close[x_, y_, tolerance_:10^-13] := Max[Abs[Flatten[x - y]]] < tolerance;

loadWilson[pairs_] := Module[{},
  ClearAll[data, de, deminus, df, dfminus, eps, thetaCh, demem, deminusmem, diagA, zeta,
    du0, dv0, uvrescalefactor, xitable, zetatable, eptable, emtable, u0ptable, u0mtable,
    scdelta, sckappa, scdeltatable, sckappatable, i, m];
  listkeywords["param"] = First /@ pairs;
  listkeywords["dmft"] = {};
  Scan[(data["param"][#[[1]]] = #[[2]]) &, pairs];
  bandrescale = paramdefaultnum["bandrescale", 1];
  Nmax = 4;
  Get[FileNameJoin[{sourceDir, "nrginit", "wilson.m"}]];
];
