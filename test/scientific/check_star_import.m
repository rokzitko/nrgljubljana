(* Exercise the real hook on adjacent doubles, without changing the physical inputs. *)
Block[{setpr, de, deminus, df, dfminus, thetaCh, theta0Ch, starRead,
       starEp, starEm, starUp, starUm, starTheta, mMAX},
  Module[{hook, work, stream, names, expected, actual},
    hook = FileNameJoin[{DirectoryName[$InputFileName], "import_star.m"}];
    work = Directory[];
    CreateDirectory["import-check"];
    SetDirectory["import-check"];
    names = {"de_pos.dat", "de_neg.dat", "du_pos.dat", "du_neg.dat", "theta.dat"};
    Do[
      stream = OpenWrite[name];
      WriteString[stream, If[name == "theta.dat", "0.2827433388230814\n",
        "0.5\n0.5000000000000001\n"]];
      Close[stream], {name, names}];
    expected = SetPrecision[Flatten[Import[#, "Table"]], Infinity] & /@ names;
    setpr[x_] := SetPrecision[x, 80];
    Get[hook];
    actual = SetPrecision[{starEp, starEm, starUp, starUm, {starTheta}}, Infinity];
    If[!TrueQ[actual == expected] || !TrueQ[starEp[[1]] < starEp[[2]]],
      Print["SCIENTIFIC_STAR_IMPORT_CHANGED_DOUBLES"];
      SetDirectory[work]; Exit[1]];
    SetDirectory[work];
    WriteString[$Output[[1]], "SCIENTIFIC_STAR_IMPORT_EXACT\n"];
  ]
];
