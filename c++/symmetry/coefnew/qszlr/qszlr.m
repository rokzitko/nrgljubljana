(*
SYMTYPE=QSZLR
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "qszlr-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
myops = Take[allops, channels];
basis = qszbasis[myops];
basis = transformtoLRvc[basis, myops, vacuum[]];
basis = bzvc2bzop @ basis;
mult[{__}] := 1;
basisprep[];

(* Problem dependent utility functions *)

Invar[{diffq_, diffsz_, diffp_}] := 
  "Invar" <> tos3[diffq, 2diffsz, diffp] <> ";";

InvarQN[{q_, sz_, p_}] := 
  "Invar" <> tos3[q, 2sz, p] <> ";";

simpl2[expr_] := Simplify[expr, 2Sz \[Element] Integers];

simpl4[expr_] := Simplify[expr, 2Sz \[Element] Integers];
simpl5[expr_] := simpl4[expr];

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

newstates1[qn:{q_, sz_, p_}, state_] := 
  newstates2[qn, {-q, -sz, p}, state];

(* Produce the expression for the new state *)

getexpr1[qn:{q_, sz_, p_}, d:{diffq_, diffsz_, diffp_}, state_] := 
  nc[coef[1], dirpdt[{{Q+diffq, Sz+diffsz, P p}, {}, r}, state]];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{q1_, sz1_, p1_}, diff1:{}, qn2:{q2_, sz2_, p2_}, diff2:{}, 
  f[izn_, _, szn_]] := 1;

(* recalcf[] stuff *)
irrule[opsz_, opp_] = {Sz -> Sz + opsz, Q -> Q+1, P -> P opp};
ircg[opsz_, opp_] = 1;

outmake[{a_, b_, c_}] := 
  "Invar" <> toc3[a /. Q -> q1, 2(b-Sz)+ssz1, c/. P -> p1];

diffqn[qn1:{q1_, sz1_, p1_}, qn2:{q2_, sz2_, p2_}] := Module[{diff},
  diff = {q2-q1, sz2-sz1, p1 p2};
  Simplify[diff, P^2 == 1]
];

CG2recalcop[qn1:{q1_, sz1_, p1_}, diff1:{},
            qn2:{q2_, sz2_, p2_}, diff2:{}, 
            {m_}] := 1;

(* Left is (Q,Sz), right is (Q', Sz')=(Qp, Szp). *)
QNrecalcop[delta:{deltaq_, deltasz_, deltap_}] := 
  (Qp = Q+deltaq; Szp = Sz+deltasz; Pp = P deltap);

RULE1recalcop[{m_}] := {};

(* Just renamed according to QNrecalcop[] definitions! *)
RULE2recalcop[{m_}] := {Q -> Qp, Sz -> Szp, P -> Pp};

CG1recalcop[{m_}] := 1;

diffsS = {
  {{0, 0, 1}, FNSINGLET<>".dat"} 
};

diffsSodd = {
  {{0, 0, -1}, FNSINGLET<>"-odd.dat"} 
};

diffsD = {
  {{-1, 1/2, 1}, FNDOUBLET<>"p.dat"},
  {{-1, -1/2, 1}, FNDOUBLET<>"m.dat"}
};

diffsDodd = {
  {{-1, 1/2, -1}, FNDOUBLET<>"diffp.dat"},
  {{-1, -1/2, -1}, FNDOUBLET<>"diffm.dat"}
};

diffsT = {
  {{0, 0, 1}, FNTRIPLET<>"s.dat"}, 
  {{0, 1, 1}, FNTRIPLET<>"p.dat"},
  {{0, -1, 1}, FNTRIPLET<>"m.dat"}
};

diffsTodd = {
  {{0, 0, -1}, FNTRIPLET<>"diffs.dat"},
  {{0, 1, -1}, FNTRIPLET<>"diffp.dat"},
  {{0, -1, -1}, FNTRIPLET<>"diffm.dat"}
};

doall[] := Module[{},
  donew[];
  dooffdiagNEW[];
  dodiag[];
  dorecalcf["QSLR"];
  recalcop[{0}, diffsS];
  recalcop[{0}, diffsSodd];
  recalcop[{1/2}, diffsD];
  recalcop[{1/2}, diffsDodd];
  recalcop[{1}, diffsT];
  recalcop[{1}, diffsTodd];
];
