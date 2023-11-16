(* 
SYMTYPE=QST
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2015-2023
*)

(* CHANGE LOG *)
(* 21.4.2016 - another attempt, using Mathematica 10.4.0 *)

$RecursionLimit = 409600;

symtype = "QST";
PREFIX = "qst";
fnprep[];

(* Make basis *)
basis = Get["one_site_basis.data"];

basis = basis /. {{n_, s_, t_}, l_List} :> {{n-3, s, t}, l};

mult[{q_, s_, t_}] := (2s+1)(2t+1);
basisprep[];

(* Problem dependent utility functions *)
<<"SSzTTz-mod.m";

fixstate[s_, sz_, t_, tz_, state_] := fixspinspin[state, s, sz, t, tz];

(* For -In.dat file: differences *)
Invar[{diffq_, diffs_, difft_}] := 
  "Invar(" <> ToString[diffq] <> ", " <> ToString[2diffs] <> ", " <> ToString[difft] <> ");";

(* For QN.dat file *)
InvarQN[{q_, s_, t_}] := 
  "Invar(" <> ToString[q] <> ", " <> ToString[2s+1] <> ", " <> ToString[t] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

(* For expressions with S, Sz, T, Tz. *)
(* Used in getexpr1[] when constructing the expressions for the NEW states *)
simpl2[expr_] := Simplify[expr /. rule, S+Sz \[Element] Integers &&
                                2S  \[Element] Integers && 
                                S >= 0 &&                 (* !!!*)
                                2Sz \[Element] Integers && 
                                T \[Element] Integers &&
                                Tz \[Element] Integers && 
                                T >= 0                     (* !!! *)
];

(* For expressions with S and T only. *)
(* Called from matrixniNEWCR, matrixni[]. (Applied on the final result in all routines, thus must
simplify expressions as far as possible.) *)
simpl4[expr_] := Simplify[expr, 2S \[Element] Integers && S >= 0 &&
                                T \[Element] Integers && T >= 0]; (* !!! *)

(* Called from matrixel[], calculation of offdiagonal matrix elements. *)
(* Note that simpl4 is also called! *)
(* simpl4 is used in more general routines once simpl5[] had already done
   the simplifications arising from Sz->Z, Tz->T. *)
simpl5[expr_] := simpl4[expr /. {Sz->S, Tz->T}]; (* Take maximal Sz and Tz!*)

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


(*  If[SIMPLIFYNEW, expr = simplifynew @ expr]; *)
(*
SIMPLIFYNEW = True;
simplifynew[x_] := Simplify[x /. Piecewise[{{val_, cond_}}, 0] -> val, S > 0 && T > 0]; (* !!! *)
*)

newstates1[qn:{q_, s_, t_}, state_] := 
  Do[newstates2[qn, {-q, diffs, difft}, state], {diffs, -s, s}, {difft, -t, t}];

(* Produce the expression for the new state *)

(* fixstate -> xyz *)

getexpr1[qn:{q_, s_, t_}, d:{diffq_, diffs_, difft_}, state_] := 
 Sum[ nc[coef @ simpl2[ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}]
                        CG[{T+difft, Tz-tz}, {t, tz}, {T, Tz}] ],
    dirpdt[{{Q+diffq, S+diffs, T+difft}, {Sz-sz,Tz-tz}, r}, xyz[s, sz, t, tz, state]]],
    {sz, -s, s}, {tz, -t, t}];

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely 
   symmetry reasons. Called from matrixniNEWCR[]. *)
(* TODO: third argument??? *)
MATRIXNICHECKSYM[i1_, i2_, _] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1 ||
                                 Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] != 1/2 ||
                                 Abs[ newqn[i1] [[3]] - newqn[i2] [[3]] ] > 1;
                                 
(* d - old site *)
(* f - new site *)

matrixniNEWCR[i1_, i2_, tz_, spin_] := Module[{expr=0},
  If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
  expr += matrixel[i1, i2, nc[f[CR, tz, spin], d[AN, tz, spin]] ];
  expr = simpl4 @ expr;
  MyPrint["matrixNEWCR === ", expr];
  expr
];

(* h.c. of the above *)
matrixniNEWAN[i1_, i2_, tz_, spin_] := Module[{expr=0},
  If[MATRIXNICHECKSYM[i1, i2], Print["skip"]; Return[0]];
  expr += matrixel[i1, i2, nc[d[CR, tz, spin], f[AN, tz, spin]] ];
  expr = simpl4 @ expr;
  MyPrint["matrixNEWAN === ", expr];
  expr
];

(* sum: hop one way and back *)
matrixniSUMCH[i1_, i2_, tz_, spin_] := Module[{},
  matrixniNEWCR[i1,i2,tz,spin] + matrixniNEWAN[i1,i2,tz,spin]
];

(* sum over tz and sz: that's the final coefficient *)
matrixniSUMCH[i1_, i2_] := Module[{expr},
  expr = Sum[matrixniSUMCH[i1,i2,tz,UP] + matrixniSUMCH[i1,i2,tz,DO], {tz, -1, 1}];
  expr = simpl4[expr];
  MyPrint["matrixniSUMCH === ", expr];
  expr
];

(* Matrix element for the HAMILTONIAN. *)
(* Called from matrixel[] *)
CGmatrixel[qn1:{q1_, s1_, t1_}, diff1:{sz1_, tz1_},
           qn2:{q2_, s2_, t2_}, diff2:{sz2_, tz2_},
           f[CR, tzn_, szn_]] := CG[{s2, sz2}, {1/2, If[szn == UP, 1/2, -1/2]}, {s1, sz1}] *
                                 CG[{t2, tz2}, {1, tzn}, {t1, tz1}];
(* TO DO: is this ever called? *)
CGmatrixel[qn1:{q1_, s1_, t1_}, diff1:{sz1_, tz1_},
           qn2:{q2_, s2_, t2_}, diff2:{sz2_, tz2_},
           f[AN, tzn_, szn_]] := CG[{s2, sz2}, {1/2, -If[szn == UP, 1/2, -1/2]}, {s1, sz1}] *
                                 CG[{t2, tz2}, {1, -tzn}, {t1, tz1}];

(* recalcf[] stuff *)
(* applied to bra in ireducel[] *)
irrule[opsz_, optz_] = 
 {S -> S + opsz, Sz -> Sz + opsz, 
  T -> T + optz, Tz -> Tz + optz,
  Q -> Q+1};
ircg[opsz_, optz_] := CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}] *
                      CG[{T, Tz}, {1, optz}, {T+optz, Tz+optz}];

outmake[{a_, b_, c_}] := 
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-S)+ss1] <> ", " <> toc[c /. T -> t1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* CG to obtain the matrix element from the irreducible matrix element *)
CG2recalcop[qn1:{q1_, s1_, t1_}, diff1:{sz1_, tz1_}, 
            qn2:{q2_, s2_, t2_}, diff2:{sz2_, tz2_},
            {m_, n_}] :=
   CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}] * CG[{t2, tz2}, {n, tz1-tz2}, {t1, tz1}];

(* Relation between QNs in bra and QNs in ket. *)
(* X are quantum numbers in bra, Xp are quantum numbers in ket. *)   
QNrecalcop[delta:{deltaq_, deltas_, deltat_}] :=
  (Qp = Q+deltaq; Sp = S+deltas; Tp = T+deltat);

(* <S,Sz=S,T,Tz=T| Sop=ms,Szop=ms,Top=mt,Tzop=mt | Sp,Spz=S-ms,Tp,Tpz=T-mt > *)

(* Transformation in bra *)
RULE1recalcop[{ms_, mt_}] := {Sz -> S, Tz -> T};

(* Transformation in ket *)
(* ms, mt are the quantum numbers of the operator *)
(* Since we are constructing the irreducible matrix elements, we take the maximum
   weight operator, i.e., S=Sz=ms, T=Tz=mt. *)
RULE2recalcop[{ms_, mt_}] := {Q -> Qp, S -> Sp, Sz -> S-ms, T -> Tp, Tz -> T-mt};

(* Bra: (Q,S,T,Sz=S,Tz=T). *)
(* Op: (S=ms,Sz=ms,T=mt,Tz=mt) *)
(* Ket: (Qp=Q+delta Q, Sp=S+delta S, Tp=T+delta T, Sz=S-ms, Tz=T-mt *)

(* CG to obtain the irreducible matrix element from the matrix element. *)
CG1recalcop[{ms_, mt_}] := CG[{Sp, S-ms}, {ms, ms}, {S, S}] * 
                           CG[{Tp, T-mt}, {mt, mt}, {T, T}];

(* Fed to QNrecalcop[] *)
diffsD = {
  {{-1, 1/2, -1}, FNDOUBLET<>"p-1.dat"},
  {{-1, -1/2, -1}, FNDOUBLET<>"m-1.dat"},
  {{-1, 1/2, 0}, FNDOUBLET<>"p0.dat"},
  {{-1, -1/2, 0}, FNDOUBLET<>"m0.dat"},
  {{-1, 1/2, +1}, FNDOUBLET<>"p+1.dat"},
  {{-1, -1/2, +1}, FNDOUBLET<>"m+1.dat"}
};

diffsS = { {{0, 0, 0}, FNSINGLET<>".dat"} };

diffsT = { 
  {{0, 0, 0}, FNTRIPLET<>"s.dat"}, 
  {{0, 1, 0}, FNTRIPLET<>"p.dat"},
  {{0, -1, 0}, FNTRIPLET<>"m.dat"}
  };

FNORBTRIPLET = "qst-orb-triplet";
  
diffsOT = { 
  {{0, 0, 0}, FNORBTRIPLET<>"s.dat"}, 
  {{0, 0, 1}, FNORBTRIPLET<>"p.dat"},
  {{0, 0, -1}, FNORBTRIPLET<>"m.dat"}
  };
 
(* dooffdiagSUMCH[] in nrgcoef.m *)

matrixopdiagSUMCH[i_, i_] := matrixop1[i, i, Sum[ number[d[tz]], {tz, -1, 1}] ];
(* dodiagSUMCH[] in nrgcoef.m *)

ClearAll[checksignfn];
checksignfn[{m_,___}] := !(m \[Element] Integers); (* True for fermionic operators *)

dorecalcfLOOP["QST"] := Module[{ch, op},
  MyPrint["dorecalcfLOOP[QST]"];

  (* d is the new site. We are thus computing the irreducible matrix elements
     <i||d||k> of the newly added site. We don't need the matrix elements of
     the previous site for this purpose!! *)
  makerecalcf[d[CR, +1, UP], {1/2, 1}, (* Sz, Tz of d[CR,+1,UP] *)
              "qst-spinup+1"];
  makerecalcf[d[CR, 0, UP], {1/2, 0},
              "qst-spinup0"];
  makerecalcf[d[CR, -1, UP], {1/2, -1},
              "qst-spinup-1"];
  makerecalcf[d[CR, +1, DO], {-1/2, 1},
              "qst-spindo+1"];
  makerecalcf[d[CR, 0, DO], {-1/2, 0},
              "qst-spindo0"];
  makerecalcf[d[CR, -1, DO], {-1/2, -1},
              "qst-spindo-1"];
];

recalcopX[m_, list_] := Module[{diffs},
  MyPrint["recalcopX[]"];
  (* Differences of conserved quantum numbers (i.e. those that characterize 
     the invariant subspaces. *)
  diffs = list[[All, 1]];
  MyPrint["diffs=", diffs];
  Scan[recalcopONEX[m, diffs, #]&, list];
];

check[{a_, b_, {d_, e_, f_}, {g_, h_, i_}, j_}] := True;
check[arg___] := (Print["Check failed:", args[arg]]; Exit[1]);

doall[] := Module[{},
(*  donew[]; Exit[]; *)
(*  loadnew[]; *)
  donew[];

  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;

  recalcopX[{0, 0}, diffsS];
  recalcopX[{0, 1}, diffsOT];
  recalcopX[{1/2, 1}, diffsD];
  recalcopX[{1, 0}, diffsT];

  dooffdiagSUMCH[];
  dodiagSUMCH[];
  dorecalcf["QST"];
];
