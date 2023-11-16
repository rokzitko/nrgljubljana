(*
SYMTYPE=QJ
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2016-2023
*)

(* CHANGE LOG *)
(* 4.3.2016 - first version *)
(* 14.3.2016 - new fixstate_qj.dat file *)

symtype = "QJ";
PREFIX = "qj";
fnprep[];

SIMPLIFYNEW = True;
simplifynew[x_] := Simplify[x /. Piecewise[{{val_, cond_}}, 0] -> val, J >= 0];

(* Make basis *)
basis = Get["one_site_qj_basis.data"];

mult[{q_, j_}] := 2j+1;
basisprep[];

checkpredznakone[expr_] := Print[{expr, predznak[expr]}];
checkpredznaksub[{qn_List, states_List}] := Map[checkpredznakone, states]
checkpredznak[basis_List] := Map[checkpredznaksub, basis];

checkpredznak[basis];

Get["fixstate_qj.data"]; (* TRICK! *)
Get["dj_qj.data"]; (* TRICK! *)

Invar[{diffq_, diffj_}] :=
  "Invar(" <> ToString[diffq] <> ", " <> ToString[2diffj] <> ");";

InvarQN[{q_, j_}] :=
  "Invar(" <> ToString[q] <> ", " <> ToString[2j+1] <> ");";

(* For expressions with J and Jz. *)
simpl2[expr_] := Simplify[expr, J+Jz \[Element] Integers &&
                                2J  \[Element] Integers && J >= 5 &&
                                2Jz \[Element] Integers];

(* For expressions with J only. *)
simpl4[expr_] := Simplify[expr, 2J \[Element] Integers && J >= 5];
simpl5[expr_] := simpl4[expr /. Jz->J]; (* Take maximal Jz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"J"->"J(jj)"] 
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"J"->"J(jjp)"] 
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"J"->"J(jj1)"] 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

newstates1[qn:{q_, j_}, state_] := 
  Do[newstates2[qn, {-q, diffj}, state], {diffj, -j, j}];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, j_}, d:{diffq_, diffj_}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{J+diffj, Jz-jz}, {j, jz}, {J, Jz}],
    dirpdt[{{Q+diffq, J+diffj}, {Jz-jz}, r}, xyz[j, jz, state]]], 
    {jz, -j, j}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely 
   symmetry reasons. *)
MATRIXNICHECKSYM[i1_, i2_] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1 &&
                              Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] > 3/2;

(* d - old site *)
(* f - new site *)

matrixniNEWCR[i1_, i2_, j_, jz_] := Module[{expr=0},
 If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
 expr += matrixel[i1, i2, nc[f[CR, j, jz], dj[AN, j, jz]] ];
 expr = simpl4 @ expr;
 expr
];

(* h.c. of the above *) (* TA NE PRISPEVA! *)
matrixniNEWAN[i1_, i2_, j_, jz_] := Module[{expr=0},
 If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
 expr += matrixel[i1, i2, nc[dj[CR, j, jz], f[AN, j, jz]] ];
 expr = simpl4 @ expr;
 expr
];

(* sum: hop one way and back *)
matrixniSUMCH[i1_, i2_, j_, jz_] := Module[{},
 MyPrint["  hop ", {j,jz}];
 matrixniNEWCR[i1,i2,j,jz] + matrixniNEWAN[i1,i2,j,jz]
];

(* sum over jz: that's the final coefficient *)
matrixniQJ[i1_, i2_, j_] := Module[{expr},
 MyPrint["matrixniQJ ",{i1,i2,j}];
 expr = Sum[matrixniSUMCH[i1,i2,j,jz], {jz, -j, j}];
 expr = simpl4[expr];
 MyPrint["RESULT === ", expr];
 expr
];

doloopQJ[function_, prefix_, filename_] := Module[
  {code, j, in1, in2, mx, koef, string},
  MyPrint["@@@@@ doloopQJ[", function, " ", prefix, " ", filename, "]"];
  filenametmp = filename <> ".tmp";
  Put["", filenametmp];
  code = {};
  For[j = 1/2, j <= 3/2, j += 1,
   MyPrint["***** doloopQJ j=", j];
   For[in1 = 1, in1 <= nrstates, in1++, (* in1 <= nrstates *)
      For[in2 = in1, in2 <= nrstates, in2++, (* in2 <= nrstates *)
          MyPrint["##### doloopQJ ##### ", in1, " ", in2];
          mx = function[in1, in2, j];
          If[mx =!= 0,
             koef = koefstr1[mx];
             nrch = If[j == 1/2, 0, 1];
             string = prefix <> "(" <>
               tos[in1] <> ", " <> tos[in2] <> ", " <>
               tos[nrch] <> ", " <> koef <> ");";
             line = {in1, in2, j, mx};
             Print[""];
             Print[line];
             Print[string];
             AppendTo[code,string];
             PutAppend[line, filenametmp];
             ];
          ];
      ];
  ];    
  store[code, filename];
];

dooffdiagQJ[] := Module[{},
    MyPrint["!!!!! dooffdiagQJ[]"];
    doloopQJ[matrixniQJ,
             "OFFDIAG",
             PREFIX <> "-offdiag.dat"];
];

CGmatrixel[qn1:{q1_, j1_}, diff1:{jz1_},
           qn2:{q2_, j2_}, diff2:{jz2_},
           f[CR, j_, jzn_]] /; (j == 1/2 || j == 3/2) && (-j <= jzn <= j) := 
             CG[{j2, jz2}, {j, jzn}, {j1, jz1}];

outmake[{a_, b_}] :=
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-J)+jj1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q1_, j1_}, diff1:{jz1_},
            qn2:{q2_, j2_}, diff2:{jz2_},
            {m_}] :=
  CG[{j2, jz2}, {m, jz1-jz2}, {j1, jz1}];
 
QNrecalcop[delta:{deltaq_, deltaj_}] :=
  (Qp = Q+deltaq; Jp = J+deltaj);

RULE1recalcop[{m_}] := {Jz -> J};
RULE2recalcop[{m_}] := {Q -> Qp, J -> Jp, Jz -> J-m};

CG1recalcop[{m_}] := CG[{Jp, J-m}, {m, m}, {J, J}];

diffsD = {
  {{-1, 1/2}, FNDOUBLET<>"p.dat"},
  {{-1, -1/2}, FNDOUBLET<>"m.dat"}
};

diffsS = { {{0, 0}, FNSINGLET<>".dat"} };

diffsT = {
  {{0, 0}, FNTRIPLET<>"s.dat"},
  {{0, 1}, FNTRIPLET<>"p.dat"},
  {{0, -1}, FNTRIPLET<>"m.dat"}
  };

FNQUAD = PREFIX <> "-quad";

diffsQ = {
  {{-1, 3/2},  FNQUAD<>"1.dat"},
  {{-1, 1/2},  FNQUAD<>"2.dat"},
  {{-1, -1/2}, FNQUAD<>"3.dat"},
  {{-1, -3/2}, FNQUAD<>"4.dat"}
};

(* recalcf[] stuff *)
(* OLD: irrule[opj_, opjz_] = {J -> J + opj, Jz -> Jz + opjz, Q -> Q+1}; *)
ircg[opj_, opjz_] = CG[{J, Jz}, {opj, opjz}, {J+opj, Jz+opjz}];
irrule[opj_, opjz_] = {J -> J + opjz, Jz -> Jz + opjz, Q -> Q+1}; (* ?? *)
ircg[opj_, opjz_] = CG[{J, Jz}, {opj, opjz}, {J+opjz, Jz+opjz}]; (* ?? *)

(* (J, Jz): kvantni stevili za ket. *)
(* (opj, opjz): d[CR, opj, opjz] -> operator f na dodanem mestu *)
(* (J+opjz, Jz+opjz) ?? *)
(* J+opjz -> ker zares spremenimo J za toliko! I1 = Invar(qp+1, jjp+NNN), kjer 
   je NNN povezan z opjz, ne pa z j! *)
(* opj=1/2 vs. opj=3/2 pa je potem upostevan z ircg[] *)

dorecalcfLOOP["QJ", nr_] := Module[{},
 MyPrint["dorecalcfLOOP[QJ]"];

 (* d is the new site. We are thus computing the irreducible matrix elements
 <i||d||k> of the newly added site. We don't need the matrix elements of
 the previous site for this purpose!! *)

If[nr == 1,
makerecalcf[dj[CR, 1/2, 1/2], {1/2, 1/2}, (* J, Jz of dj[CR, 1/2, 1/2] *)
"qj-spin_j1_2-jz1_2"];
];

If[nr == 2,
makerecalcf[dj[CR, 1/2, -1/2], {1/2, -1/2},
"qj-spin_j1_2-jz-1_2"];
];

If[nr ==3,
makerecalcf[dj[CR, 3/2, 3/2], {3/2, 3/2},
"qj-spin_j3_2-jz3_2"];
];

If[nr == 4,
makerecalcf[dj[CR, 3/2, 1/2], {3/2, 1/2},
"qj-spin_j3_2-jz1_2"];
];

If[nr == 5,
makerecalcf[dj[CR, 3/2, -1/2], {3/2, -1/2},
"qj-spin_j3_2-jz-1_2"];
];

If[nr == 6,
makerecalcf[dj[CR, 3/2, -3/2], {3/2, -3/2},
"qj-spin_j3_2-jz-3_2"];
];
];

matrixopdiagSUMCH[i_, i_] := matrixop1[i, i, Sum[ number[d[tz]], {tz, -1, 1}] ];
(* dodiagSUMCH[] in nrgcoef.m *)

doall[] := Module[{},
  (* donew[]; *)
  loadnew[];

  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;

  recalcop[{0}, diffsS];
  recalcop[{1/2}, diffsD];
  recalcop[{1}, diffsT];
  recalcop[{3/2}, diffsQ];

  dodiagSUMCH[];
  dooffdiagQJ[];
  dorecalcf["QJ", 1];
];
