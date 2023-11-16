(*
SYMTYPE=SPSU2T
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2015-2023
*)

(* CHANGE LOG *)
(* 22.10.2015 - first version *)
(* 25.10.2015 - fixes merged *)
(* 26.10.2015 - matrixanomalous[] fixed operator names *)
(* 26.10.2015 - idem for matrixopisospinx[] *)

symtype = "SPSU2T";
PREFIX = "spsu2t";
fnprep[];

(* Make basis *)
basis = Get["one_site_basis.data"];

basis = Map[{ Drop[First[#],1], Last[#] }&, basis]; (* Drop charge quantum number *)
basis = mergebasis[basis];

mult[{s_, t_}] := (2s+1)(2t+1);
basisprep[];

(* Problem dependent utility functions *)
<<"SSzTTz-mod.m";

fixstate[s_, sz_, t_, tz_, state_] := fixspinspin[state, s, sz, t, tz];

(* For -In.dat file: differences *)
Invar[{diffs_, difft_}] :=
  "Invar(" <> ToString[2diffs] <> ", " <> ToString[difft] <> ");";

(* For QN.dat file *)
InvarQN[{s_, t_}] :=
  "Invar(" <> ToString[2s+1] <> ", " <> ToString[t] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

(* For expressions with S, Sz, T, Tz. *)
(* Used in getexpr1[] when constructing the expressions for the NEW states *)
simpl2[expr_] := Simplify[expr /. rule, S+Sz \[Element] Integers &&
                                2S  \[Element] Integers && 
                                S >= 4 &&                 (* !!!*)
                                2Sz \[Element] Integers && 
                                T \[Element] Integers &&
                                Tz \[Element] Integers && 
                                T >= 4                     (* !!! *)
];

(* For expressions with S and T only. *)
(* Called from matrixniNEWCR, matrixni[]. (Applied on the final result in all routines, thus must
simplify expressions as far as possible.) *)
simpl4[expr_] := Simplify[expr, 2S \[Element] Integers && S >= 4 &&
                                T \[Element] Integers && T >= 4]; (* !!! *)

(* Called from matrixel[], calculation of offdiagonal matrix elements. *)
(* Note that simpl4 is also called! *)
(* simpl4 is used in more general routines once simpl5[] had already done
   the simplifications arising from Sz->Z, Tz->T. *)
simpl5[expr_] := simpl4[expr //. {Sz->S, Tz->T}]; (* Take maximal Sz and Tz!*)

(* Called from dogeneralloop[] *)
koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss)"] 
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"S"->"S(ssp)"] 
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss1)"] 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

newstates1[qn:{s_, t_}, state_] := 
  Do[newstates2[qn, {diffs, difft}, state], {diffs, -s, s}, {difft, -t, t}];

(* Produce the expression for the new state *)

getexpr1[qn:{s_, t_}, d:{diffs_, difft_}, state_] := 
 Sum[ nc[coef @ simpl2[ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}]
                        CG[{T+difft, Tz-tz}, {t, tz}, {T, Tz}] ],
    dirpdt[{{S+diffs, T+difft}, {Sz-sz,Tz-tz}, r}, fixstate[s, sz, t, tz, state]]], 
    {sz, -s, s}, {tz, -t, t}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely
   symmetry reasons. Called from matrixniNEWCR[]. *)
MATRIXNICHECKSYM[i1_, i2_] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1/2 ||
                              Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] > 1;

(* d - old site *)
(* f - new site *)

matrixniX[i1_, i2_, tz_, spin_] := Module[{op,expr=0},
  If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
  op=nc[f[CR, tz, spin], d[AN, tz, spin]];
  Print["X i1=", i1, " i2=", i2, " op=", op];
  expr += matrixel[i1, i2, op];
  expr = Simplify[simpl4 @ expr];
  MyPrint["X === ", expr];
  expr
];

(* sum over tz and sz: that's the final coefficient *)
matrixniSUMCH[i1_, i2_] := Module[{expr},
  expr = Sum[matrixniX[i1,i2,tz,UP] + matrixniX[i1,i2,tz,DO], {tz, -1, 1}];
  expr = Simplify[simpl4[expr]];
  MyPrint["matrixniSUMCH xy === ", expr];
  expr
];

(* Matrix element for the HAMILTONIAN. *)
(* Called from matrixel[] when building the offdiag parts using
   matrixniNEW* defined just above. *)
CGmatrixel[qn1:{s1_, t1_}, diff1:{sz1_, tz1_},
           qn2:{s2_, t2_}, diff2:{sz2_, tz2_},
           f[CR, tzn_, szn_]] := CG[{s2, sz2}, {1/2, If[szn == UP, 1/2, -1/2]}, {s1, sz1}] *
                                 CG[{t2, tz2}, {1, tzn}, {t1, tz1}];

CGmatrixel[qn1:{s1_, t1_}, diff1:{sz1_, tz1_},
           qn2:{s2_, t2_}, diff2:{sz2_, tz2_},
           f[AN, tzn_, szn_]] := (MyPrint["oops"]; Exit[1]) (* Bug trap *)

(* recalcf[] stuff *)
(* applied to bra in ireducel[] *)
irrule[opsz_, optz_] = 
 {S -> S + opsz, Sz -> Sz + opsz, 
  T -> T + optz, Tz -> Tz + optz};
ircg[opsz_, optz_] := CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}] *
                      CG[{T, Tz}, {1, optz}, {T+optz, Tz+optz}];

outmake[{b_, c_}] := 
  "Invar(" <> toc[2(b-S)+ss1] <> ", " <> toc[c /. T -> t1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{s1_, t1_}, diff1:{sz1_, tz1_}, 
            qn2:{s2_, t2_}, diff2:{sz2_, tz2_},
            {m_, n_}] :=
   CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}] *
   CG[{t2, tz2}, {n, tz1-tz2}, {t1, tz1}];
 
(* Relation between QNs in bra and QNs in ket. *)
(* X are quantum numbers in bra, Xp are quantum numbers in ket. *)   
QNrecalcop[delta:{deltas_, deltat_}] :=
  (Sp = S+deltas; Tp = T+deltat);

(* Transformation in bra *)
RULE1recalcop[{ms_, mt_}] := {Sz -> S, Tz -> T};

(* Transformation in ket *)
(* ms, mt are the quantum numbers of the operator *)
(* Since we are constructing the irreducible matrix elements, we take the maximum
   weight operator, i.e., S=ms, T=mt, Sz=S, Tz=T. *)
RULE2recalcop[{ms_, mt_}] := {S -> Sp, Sz -> S-ms, T -> Tp, Tz -> T-mt};

(* Bra: (S,T,Sz=S,Tz=T). *)
(* Op: (S=ms,Sz=ms,T=mt,Tz=mt) *)
(* Ket: (Sp=S+delta S, Tp=T+delta T, Sz=S-ms, Tz=T-mt *)

CG1recalcop[{ms_, mt_}] := CG[{Sp, S-ms}, {ms, ms}, {S, S}] * 
                           CG[{Tp, T-mt}, {mt, mt}, {T, T}];

(* Fed to QNrecalcop[] *)
diffsD = {
  {{1/2, -1}, FNDOUBLET<>"p-1.dat"},
  {{-1/2, -1}, FNDOUBLET<>"m-1.dat"},
  {{1/2, 0}, FNDOUBLET<>"p0.dat"},
  {{-1/2, 0}, FNDOUBLET<>"m0.dat"},
  {{1/2, +1}, FNDOUBLET<>"p+1.dat"},
  {{-1/2, +1}, FNDOUBLET<>"m+1.dat"}
};

diffsS = { {{0, 0}, FNSINGLET<>".dat"} };

diffsT = { 
  {{0, 0}, FNTRIPLET<>"s.dat"}, 
  {{1, 0}, FNTRIPLET<>"p.dat"},
  {{-1, 0}, FNTRIPLET<>"m.dat"}
  };

(* No sum over channels. Performed in function[]. *)
dogeneralloopSUMCH[function_, prefix_, filename_, upper_:False, INCLUDEZEROONDIAG_:False] := Module[
  {code, in1, startin2, in2, mx, koef, string},
  MyPrint["dogeneralloop[", function, " ", prefix, " ", filename, "]"];
  code = {};
  For[in1 = 1, in1 <= nrstates, in1++, (* in1 <= nrstates *)
    startin2 = If[upper, in1, 1];
    For[in2 = startin2, in2 <= nrstates, in2++, (* in2 <= nrstates *)
      MyPrint["##### dogeneralloopSUMCH ##### ", in1, " ", in2];
      mx = function[in1, in2];
      If[mx =!= 0 || (INCLUDEZEROONDIAG && in1 == in2),
        koef = koefstr1[mx];
        string = prefix <> "(" <>
          tos[in1] <> ", " <> tos[in2] <> ", " <> koef <> ");";
        Print[string];
        AppendTo[code,string];
      ];
    ];
  ];
  store[code, filename];
];

(* Sum over channels performed in the inner-most loop, i.e., 
in matrixniSUMCH[], rather than in dogeneralloop[]. *)
dooffdiagSUMCH[] := Module[{},
 MyPrint["dooffdiagSUMCH[]"];
 dogeneralloopSUMCH[matrixniSUMCH, (* different low-level routine! *)
 "OFFDIAG",
 PREFIX <> "-offdiag.dat", False]; (* False here!!! *)
];

matrixopdiagSUMCH[i_, i_] := matrixop1[i, i, Sum[ number[d[tz]], {tz, -1, 1}] ];

dodiagSUMCH[] := Module[{},
 code = {};
 For[in1 = 1, in1 <= nrstates, in1++,
   mx = matrixopdiagSUMCH[in1, in1];
   string = "DIAG(" <> tos[in1] <> ", " <> toc[mx] <> ");";
   Print[string];
   AppendTo[code,string];
 ];
 store[code,PREFIX<>"-diag.dat"];
];

matrixanomalous[i1_, i2_, tz_] := Module[{expr=0},
  expr += matrixel[i1, i2, nc[f[CR,tz,UP], d[CR,tz,DO]]];
  expr -= matrixel[i1, i2, nc[f[CR,tz,DO], d[CR,tz,UP]]];
  If[expr =!= 0,
  MyPrint["a #### ", expr];
  ];
  simpl4 @ expr
];

(* sum over tz and sz: that's the final coefficient *)
matrixanomalousSUMCH[i1_, i2_] := Module[{expr},
  expr = Sum[matrixanomalous[i1,i2,tz], {tz, -1, 1}];
  expr = simpl4[expr];
  MyPrint["matrixanomalousSUMCH === ", expr];
  expr
];

doanomalousSUMCH[] := Module[{},
  MyPrint["doanomalousSUMCH[]"];
  dogeneralloopSUMCH[matrixanomalousSUMCH, 
  "ANOMALOUS",
  PREFIX <> "-anomalous.dat", True];
];

(* fixed *)
matrixopisospinx[i1_, i2_, tz_] := Module[{Ix},
  Ix = simpl5 @ isospinx[d[tz]]; (* !! *)
  matrixop1[i1, i2, Ix] (* note, matrixop1, just like for dodiag[]!! *)
];

(* sum over tz and sz: that's the final coefficient *)
matrixisospinxSUMCH[i1_, i2_] := Module[{expr},
  expr = Sum[matrixopisospinx[i1,i2,tz], {tz, -1, 1}];
  expr = simpl5[expr];
  MyPrint["matrixisospinxSUMCH === ", expr];
  expr
];

doisospinxSUMCH[] := Module[{},
  MyPrint["doisospinxSUMCH[]"];
  dogeneralloopSUMCH[matrixisospinxSUMCH,
  "ISOSPINX",
  PREFIX <> "-isospinx.dat", True];
];


ClearAll[checksignfn];
checksignfn[{m_,___}] := !(m \[Element] Integers); (* True for fermionic operators *)

dorecalcfLOOP["SPSU2T"] := Module[{ch, op},
  MyPrint["dorecalcfLOOP[SPSU2T]"];

  (* d is the new site. We are thus computing the irreducible matrix elements
     <i||d||k> of the newly added site. We don't need the matrix elements of
     the previous site for this purpose!! *)
  makerecalcf[d[CR, +1, UP], {1/2, 1}, (* Sz, Tz of d[CR,+1,UP] *)
              "spsu2t-spinup+1"];

  makerecalcf[d[CR, 0, UP], {1/2, 0},
              "spsu2t-spinup0"];

  makerecalcf[d[CR, -1, UP], {1/2, -1},
              "spsu2t-spinup-1"];

  makerecalcf[d[CR, +1, DO], {-1/2, 1},
              "spsu2t-spindo+1"];

  makerecalcf[d[CR, 0, DO], {-1/2, 0},
              "spsu2t-spindo0"];

  makerecalcf[d[CR, -1, DO], {-1/2, -1},
              "spsu2t-spindo-1"];
];

doall[] := Module[{},
  donew[];

  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;

  dooffdiagSUMCH[];
  dodiagSUMCH[];
  doisospinxSUMCH[];
  dorecalcf["SPSU2T"];
  recalcop[{1/2, 1}, diffsD];
  recalcop[{0, 0}, diffsS];
  recalcop[{1, 0}, diffsT];

  doanomalousSUMCH[];
];
