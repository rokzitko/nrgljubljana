(*
   Parameter file parsing code 
   Rok Zitko, rok.zitko@ijs.si, 2007
*)

(* Helper functions *)
parseisgroup[s_String] := StringTake[s, 1] == "[" && StringTake[s, -1] == "]";
parseisgroup[""] = False;
parsestrip[s_String] := Fold[StringDrop, s, {1, -1}];

(* Strip leading and trailing whitespace in Mma 5.0 compatible way *)
stripleading[s_String] := stripleading[StringDrop[s, 1]] /; 
  StringTake[s, 1] == " ";
stripleading[s_String] := s;
stripleading[""] = "";
striptrailing[s_String] := striptrailing[StringDrop[s, -1]] /; 
  StringTake[s, -1] == " ";
striptrailing[s_String] := s;
striptrailing[""] = "";
stripws[s_String] := stripleading[striptrailing[s]];
stripws[""] := ""; (* req'd in mma6 *)

(* Parameters starting with ! sign are evaluated as Mathematica
expressions! *)
klicaj[str_String] := StringLength[str] >= 1 && StringTake[str, 1] == "!";
evalstr[str_String] /; klicaj[str] := ToString[ToExpression[stripws @ StringDrop[str, 1]]];
evalstr[str_String] := str;

parse::nogroup = "Line `1` not contained in a group.";

parse[filename_]:=Module[{l, len, group, i, line},
  l = Import[filename, "Lines"];
  If[l == $Failed, Return[$Failed]];

  ClearAll[data];

  (* Strip leading/trailing whitespace *)
  l = Map[stripws, l];

  (* Drop blank lines and comments *)
  l = Select[l, # != "" &];
  l = Select[l, StringTake[#, 1] != "#" &];

  len = Length[l];
  group = Null;
  For[i = 1, i <= len, i++,
    line = l[[i]];
    If[parseisgroup[line],
      (* Header line *)
      group = parsestrip[line];
      listdata[group] = {}, (* !! *)
    (* elese *)
      (* Key-value line *)
      If[group == Null, Message[parse::nogroup, line]; Return[$Failed]];
      equalpos = StringPosition[line, "="][[1,1]];
      key = StringTake[line, {1, equalpos-1}];
      key = stripws[key];
      value = StringTake[line, {equalpos+1, StringLength[line]}];
      value = stripws[value];
      data[group][key] = value;
      AppendTo[listdata[group], {key,value}];
    ];
  ];
];

(* Get parameter with key 'key' in the group 'group'. *)
(* Parameters starting with ! sign are evaluated as Mathematica
   expressions! *)
param[key_, group_:"param"] := evalstr[ data[group][key] ];

(* For numerical quantities ImportString[] is used to handle C-like numeric
expression such as 1e-3. Note: only the first number is returned! *)
importnum[str_String] := First @ ImportString[ str, "List" ];
paramnum[key_, group_:"param"] := importnum @ param[key, group];

(* For logic quantities we allow several ways of specifying True. *)
parambool[key_, group_:"param"] := 
  MemberQ[{"Yes", "yes", "True", "true", "TRUE", "1"}, param[key, group] ];

paramexists[key_, group_:"param"] := ValueQ[ data[group][key] ];

paramdefault[key_, default_, group_:"param"] :=
  If[paramexists[key, group], param[key, group], default];

paramdefaultnum[key_, default_, group_:"param"] :=
  If[paramexists[key, group], paramnum[key, group], default];

paramdefaultbool[key_, default_, group_:"param"] :=
  If[paramexists[key, group], parambool[key, group], default];
