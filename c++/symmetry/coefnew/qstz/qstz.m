(* 
SYMTYPE=QSTZ
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2016-2023
*)

symtype = "QSTZ";
PREFIX = "qstz";
fnprep[];

(* Make basis *)
basis = Get["one_site_qstz_basis.data"];

basis = basis /. {{n_, s_, tz_}, l_List} :> {{n-3, s, tz}, l};

mult[{q_, s_, tz_}] := (2s+1);
basisprep[];

(* Problem dependent utility functions *)
<<"SSzTTz-mod-v2.m";

SIMPLIFYNEW = True;
simplifynew[x_] := Simplify[x /. Piecewise[{{val_, cond_}}, 0] -> val, S > 0]; (* !!! *)

fixstate[s_, sz_, state_] := fixspins[state, s, sz];

(* For -In.dat file: differences *)
Invar[{diffq_, diffs_, difftz_}] :=
  "Invar(" <> ToString[diffq] <> ", " <> ToString[2diffs] <> ", " <> ToString[difftz] <> ");";

(* For QN.dat file *)
InvarQN[{q_, s_, tz_}] :=
  "Invar(" <> ToString[q] <> ", " <> ToString[2s+1] <> ", " <> ToString[tz] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 5]; (* ! *)
SIMPLIFYNEW = True;

(* Used in getexpr1[] when constructing the expressions for the NEW states *)
simpl2[expr_] := Simplify[expr /.rule, S+Sz \[Element] Integers &&
                                2S  \[Element] Integers && 
                                S >= 5 &&                 (* !!!*)
                                2Sz \[Element] Integers && 
                                Tz \[Element] Integers
];

(* Called from matrixniNEWCR, matrixni[]. (Applied on the final result in all routines, thus must
simplify expressions as far as possible.) *)
simpl4[expr_] := Simplify[expr, 2S \[Element] Integers && S >= 5]; (* !!! *)

(* Called from matrixel[], calculation of offdiagonal matrix elements. *)
(* Note that simpl4 is also called! *)
(* simpl4 is used in more general routines once simpl5[] had already done
   the simplifications arising from Sz->Z. *)
simpl5[expr_] := simpl4[expr /. {Sz->S}]; (* Take maximal Sz!*)

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


newstates1[qn:{q_, s_, tz_}, state_] := 
  Do[newstates2[qn, {-q, diffs, -tz}, state], {diffs, -s, s}];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, s_, tz_}, d:{diffq_, diffs_, difftz_}, state_] := 
 Sum[ nc[coef @ simpl2[ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}] ],
    dirpdt[{{Q+diffq, S+diffs, Tz+difftz}, {Sz-sz}, r}, fixstate[s, sz, state]]], 
    {sz, -s, s}];

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
CGmatrixel[qn1:{q1_, s1_, tz1_}, diff1:{sz1_},
           qn2:{q2_, s2_, tz2_}, diff2:{sz2_},
           f[CR, tzn_, szn_]] := CG[{s2, sz2}, {1/2, If[szn == UP, 1/2, -1/2]}, {s1, sz1}];
           
(* TO DO: is this ever called? *)
CGmatrixel[qn1:{q1_, s1_, tz1_}, diff1:{sz1_},
           qn2:{q2_, s2_, tz2_}, diff2:{sz2_},
           f[AN, tzn_, szn_]] := (Print["ERROR"]; Exit[1];);

(* recalcf[] stuff *)
(* applied to bra in ireducel[] *)
irrule[opsz_, optz_] = 
 {S -> S + opsz, Sz -> Sz + opsz, Tz -> Tz + optz,
  Q -> Q+1};

ircg[opsz_, optz_] := CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

outmake[{a_, b_, c_}] := 
  "Invar(" <> toc[a /. Q -> q1] <> ", "<> toc[2(b-S)+ss1] <> ", " <> toc[c /. Tz -> tz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

CG2recalcop[qn1:{q1_, s1_, tz1_}, diff1:{sz1_}, 
            qn2:{q2_, s2_, tz2_}, diff2:{sz2_},
            {m_}] :=
   CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];

(* Relation between QNs in bra and QNs in ket. *)
(* X are quantum numbers in bra, Xp are quantum numbers in ket. *)   
QNrecalcop[delta:{deltaq_, deltas_, deltatz_}] :=
  (Qp = Q+deltaq; Sp = S+deltas; Tzp = Tz+deltatz);

(* Transformation in bra *)
RULE1recalcop[{ms_}] := {Sz -> S};

(* Transformation in ket *)
(* ms, mt are the quantum numbers of the operator *)
(* Since we are constructing the irreducible matrix elements, we take the maximum
   weight operator, i.e., S=ms, T=mt, Sz=S, Tz=T. *)
RULE2recalcop[{ms_}] := {Q -> Qp, S -> Sp, Sz -> S-ms, Tz -> Tzp};

(* Bra: (Q,S,T,Sz=S,Tz=T). *)
(* Op: (S=ms,Sz=ms,T=mt,Tz=mt) *)
(* Ket: (Qp=Q+delta Q, Sp=S+delta S, Tp=T+delta T, Sz=S-ms, Tz=T-mt *)

CG1recalcop[{ms_}] := CG[{Sp, S-ms}, {ms, ms}, {S, S}];

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
 
(* dooffdiagSUMCH[] in nrgcoef.m *)

matrixopdiagSUMCH[i_, i_] := matrixop1[i, i, Sum[ number[d[tz]], {tz, -1, 1}] ];
(* dodiagSUMCH[] in nrgcoef.m *)

ClearAll[checksignfn];
checksignfn[{m_,___}] := !(m \[Element] Integers); (* True for fermionic operators *)

dorecalcfLOOP["QSTZ"] := Module[{ch, op},
  MyPrint["dorecalcfLOOP[QSTZ]"];

  (* d is the new site. We are thus computing the irreducible matrix elements
     <i||d||k> of the newly added site. We don't need the matrix elements of
     the previous site for this purpose!! *)
  makerecalcf[d[CR, +1, UP], {1/2, 1}, (* Sz, Tz of d[CR,+1,UP] *)
              "qstz-spinup+1"];

  makerecalcf[d[CR, 0, UP], {1/2, 0},
              "qstz-spinup0"];

  makerecalcf[d[CR, -1, UP], {1/2, -1},
              "qstz-spinup-1"];

  makerecalcf[d[CR, +1, DO], {-1/2, 1},
              "qstz-spindo+1"];

  makerecalcf[d[CR, 0, DO], {-1/2, 0},
              "qstz-spindo0"];

  makerecalcf[d[CR, -1, DO], {-1/2, -1},
              "qstz-spindo-1"];
];

recalcopX[m_, list_] := Module[{diffs},
  MyPrint["recalcopX[]"];

  (* Differences of conserved quantum numbers (i.e. those that characterize 
     the invariant subspaces. *)
  diffs = list[[All, 1]];

  Scan[recalcopONEX[m, diffs, #]&, list];
];

recalcopONEX[m_, diffs_, {diffone_, filename_}] := Module[
  {seznam, checksign, i, j, expr, string},
  MyPrint["recalcopONEX[] ", m, diffs, {diffone, filename}, {fn2, fn3}];

  seznam = {}; (* In this array we accumulate the results *)

  checksign = checksignfn[m]; (* True for fermionic operators *)

  QNrecalcop[ diffone ]; (* Relation between QNs in bra and QNs in ket. *)

  For[i = 1, i <= nrstates, i++,  (* All combinations are considered! *)
    st1 = new[i] /. RULE1recalcop[m];
    If[checksign,
      st1 = prepbrasign @ st1,
      st1 = prepbra @ st1
    ];
    st1 = st1 /. nc[coef[0], ___] -> 0;

    For[j = 1, j <= nrstates, j++,
      st2 = prepket[ new[j] /. RULE2recalcop[m]];
      st2 = st2 /. nc[coef[0], ___] -> 0;

      expr = factorcoef @ nc[st1, st2];
      expr = expr //. HoldPattern[ nc[b___, old[x__], c___, old[y__], d___] ] :>
        pdt[{x}, {y}] mynewbk[nc[b], nc[c], nc[d]];

      MyPrint[{i, j}];

      cg1 = CG1recalcop[m];

      expr = expr /. pdt[{qn1_, diff1_, r1_}, {qn2_, diff2_, r2_}] :>
        checkdiff[diffqn[qn1, qn2], diffs] ireduc[i, j, qn1, qn2] *
        CG2recalcop[qn1, diff1, qn2, diff2, m];

      expr = expr/cg1;
      expr = simpl4[expr];

      If[expr =!= 0,
        string = expr /. faktor_. ireduc[i1_, ip_, IN1_, INp_] :>
          "{ " <> tos[i1] <> ", " <> tos[ip] <> ", " <>
          outmake[IN1] <> ", " <> outmake[INp] <>
          ", " <> koefstr3[ faktor ] <> "},";
        Print[string];

        res = expr /. {
          faktor_. ireduc[i1_, ip_, IN1_, INp_] :> {i1, ip, IN1, INp, faktor}
        };

        AppendTo[seznam, string];
      ];
    ];  (* j *)
  ]; (* i *)
  store[seznam, filename]; (* SIDEEFFECT: store[] also computes LEN! *)
];

doall[] := Module[{},
  donew[];

  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;

  dorecalcf["QSTZ"];

  dooffdiagSUMCH[];
  dodiagSUMCH[];

  recalcopX[{0}, diffsS];
  recalcopX[{1/2}, diffsD];
  recalcopX[{1}, diffsT];
];
