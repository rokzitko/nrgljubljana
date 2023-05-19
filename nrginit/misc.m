(** Misc functions **)

(* Generalization of Print[] which wraps all strings in StandardForm[]
format directives to remove the quotation marks. *)
Print2[l__] := Print @@ ({l} /. x_String -> StandardForm[x]);

(* Print a warning/error message at specified verbosity level *)
MyPrint[msg__] := If[DEBUG >= 1, Print2[msg]];
MyPrintForm[form_, msg___] := If[DEBUG >= 1,
  Print2 @ ToString @ StringForm[form, Sequence @@ (InputForm /@ {msg})]
];
MyVPrint[verbosity_, msg__] := If[DEBUG >= verbosity, Print2[msg]];

(* Prints an error message and exits *)
If[!ValueQ[ExitOnError], ExitOnError = True];
MyError[msg__] := (Print2[msg]; Print2["Aborting.\n"];
  If[ExitOnError, Exit[1] ]);

(* Prints a warning message -- does not exit *)
MyWarning[msg__] := Print2[msg];

SetAttributes[MyAssert, HoldFirst];
MyAssert[expr_] := If[Evaluate[expr] =!= True, MyError["Assertion failed: ", Hold[expr] ]];

(* Convert to a C-language compatible string, truncate to 10 significant digits *)
cstr10[x_] := ToString[N[x, 10], CForm];

c10[x_] := CForm[N[x,10]]; (* used in dmft.m *)

(* Timing code. Used for profiling and benchmarking. *)
ClearAll[timingdata, time0];
timingdata[_] = 0;
timestart[name_] := (time0[name] = AbsoluteTime[]);
timeadd[name_] := Module[{time1, timediff},
  time1 = AbsoluteTime[];
  timediff = time1-time0[name];
  timingdata[name] += timediff;
];
timereport[] := Module[{t},
  MyPrint["Timing report"];
  t = Map[{#[[1,1,1]], #[[2]]}&, DownValues[timingdata]];
  Scan[MyPrint, t];
];

(* Load a module, Get::noopen message is suppressed. *)
silentGet[x__] := Module[{ret},
  Off[Get::noopen];
  ret = Get[x];
  On[Get::noopen];
  ret
];

(* Load an external module (another .m) file. *)
loadmodule[filename_String, exitonfailure_:True] := Module[{ret},
  MyPrint["Loading module ", filename];
  ret = silentGet[filename, Path -> modulespath];
  If[ret === $Failed,
    MyWarning["Can't load " <> filename <> ". " <>
      If[exitonfailure, "Aborting.", "Continuing."]];
    If[exitonfailure == True, Exit[1]];
  ];
  ret (* The result is returned! *)
];

(* External hook: if the variable 'var' ends with a ".m" suffix,
   the corresponding external module is called. *)
hook[var_] := If[StringTake[var,-2] == ".m", loadmodule[var, True]];

(* External hook with a named external module *)
hookfile[keyword_] := Module[{filename},
  filename = paramdefault[keyword, ""];
  If[filename =!= "",
    loadmodule[filename, True]
  ];
];