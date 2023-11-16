(*
SYMTYPE=QSLR
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "qslr-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
myops = Take[allops, channels]
basis = qsbasis[myops];
basis = transformtoLRvc[basis, myops, vacuum[]];
basis = bzvc2bzop @ basis;
mult[{q_, s_, p_}] := 2s+1;
basisprep[];

(* Problem dependent utility functions *)
fixstate[s_, sz_, state_] := fixspin[s, sz, state];

Invar[{diffq_, diffs_, diffp_}] := "Invar" <> tos3[diffq, 2diffs, diffp] <> ";"; 
InvarQN[{q_, s_, p_}] := "Invar" <> tos3[q, 2s+1, p] <> ";";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule, S+Sz \[Element] Integers];

simpl4[expr_] := FullSimplify[expr, 2S \[Element] Integers && S>2];
simpl5[expr_] := simpl4[expr /. Sz->S]; (* Take maximal Sz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str, "S"->"S(ss)"] 
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

newstates1[qn:{q_, s_, p_}, state_] := 
  Do[newstates2[qn, {-q, diffs, p}, state], {diffs, -s, s}];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, s_, p_}, d:{_, diffs_, _}, state_] := 
  Sum[ nc[coef @ simpl2 @ CG[{S+diffs, Sz-sz}, {s, sz}, {S, Sz}],
    dirpdt[{{Q-q, S+diffs, p P}, {Sz-sz}, r}, fixstate[s, sz, state]]], 
    {sz, -s, s}];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{q1_, s1_, p1_}, diff1:{sz1_}, qn2:{q2_, s2_, p2_}, diff2:{sz2_}, 
  f[izn_, _, szn_]] := CG[{s2, sz2}, {1/2, spinofop[izn,szn]}, {s1, sz1}];

(* recalcf[] stuff *)
irrule[opsz_, opp_] = {S -> S + opsz, Sz -> Sz + opsz, Q -> Q+1, P -> P opp};
ircg[opsz_, opp_] = CG[{S, Sz}, {1/2, opsz}, {S+opsz, Sz+opsz}];

outmake[{a_, b_, c_}] := 
  "Invar" <> toc3[a /. Q -> q1, 2(b-S)+ss1, c /. P->p1];

diffqn[qn1:{q1_, s1_, p1_}, qn2:{q2_, s2_, p2_}] := Module[{diff},
  diff = {q2-q1, s2-s1, p1 p2};
  Simplify[diff, P^2 == 1]
];

CG2recalcop[qn1:{q1_, s1_, p1_}, diff1:{sz1_}, 
            qn2:{q2_, s2_, p2_}, diff2:{sz2_}, 
            {m_}] :=
  CG[{s2, sz2}, {m, sz1-sz2}, {s1, sz1}];

QNrecalcop[delta:{deltaq_, deltas_, deltap_}] := 
  (Qp = Q+deltaq; Sp = S+deltas; Pp = P deltap);

RULE1recalcop[{m_}] := {Sz -> S};
RULE2recalcop[{m_}] := {Q -> Qp, S -> Sp, Sz -> S-m, P -> Pp};

CG1recalcop[{m_}] := CG[{Sp, S-m}, {m, m}, {S, S}];

diffsD = {
  {{-1, 1/2, 1}, FNDOUBLET<>"p.dat"},
  {{-1, -1/2, 1}, FNDOUBLET<>"m.dat"}
};

diffsDodd = {
  {{-1, 1/2, -1}, FNDOUBLET<>"diffp.dat"},
  {{-1, -1/2, -1}, FNDOUBLET<>"diffm.dat"}
};

diffsS = {
  {{0, 0, 1}, FNSINGLET<>".dat"}
};

diffsSodd = {
  {{0, 0, -1}, FNSINGLET<>"-odd.dat"}
};

diffsT = {
  {{0,  0, 1}, FNTRIPLET<>"s.dat"},
  {{0,  1, 1}, FNTRIPLET<>"p.dat"},
  {{0, -1, 1}, FNTRIPLET<>"m.dat"}
};

diffsTodd = {
  {{0,  0, -1}, FNTRIPLET<>"diffs.dat"},
  {{0,  1, -1}, FNTRIPLET<>"diffp.dat"},
  {{0, -1, -1}, FNTRIPLET<>"diffm.dat"}
};

doall[] := Module[{},
  donew[];
  dooffdiagNEW[];
  dodiag[];
  dorecalcf["QSLR"];

  recalcop[{1/2}, diffsD];
  recalcop[{1/2}, diffsDodd];
  recalcop[{0}, diffsS];
  recalcop[{0}, diffsSodd]; 
  recalcop[{1}, diffsT];
  recalcop[{1}, diffsTodd];
];
