(*
SYMTYPE=ISOSZLR
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2009-2023
*)

PREFIX = "isoszlr-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
myops = Take[allops, channels];
basis = Simplify @ qszbasisvc @ myops;

snegrealconstants[n, NN];
snegrealconstants[nnn];

(* IMPORTANT: the following choice of nnop corresponds
to a particular partition of the lattice sites into two
sets! *)

nnop[a[]] = n;
nnop[b[]] = n;
nnop[c[]] = n;
nnop[d[]] = n;
nnop[e[]] = n;

Tminus = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, myops];
Print[Tminus];

basis = bzvc2bzop @ transformQStoIS[basis];

(* Added for LR *)
basis = transformtoLRvc[basis, myops, vacuum[]];
basis = bzvc2bzop @ basis;

mult[{i_, __}] := 2i+1;
basisprep[];

(* Simplify n: simply putting n \[Element] Integers interferes with 
snegrealconstants[n] *)

simpln[expr_] :=
  Simplify[(expr /. n->nnn), nnn \[Element] Integers] /. nnn->n;

simplNN[expr_] :=
  Simplify[(expr /. NN->nnn), nnn \[Element] Integers] /. nnn->NN;

(* Reduce isospin by one unit *)

isospindown[op_] := Module[{tmp},
  tmp = nc[Tminus, op];
  tmp = zeroonvac[tmp];
  If[tmp =!= 0,
    tmp = tmp / Sqrt[simpln[scalarproductop[tmp,tmp]]],
    0
  ]
];

(* Problem dependent utility functions *)

(* This will apply to f on site N+1 (that's being added), 
therefore n -> (NN+1). *)

fixstate[i_, iz_, state_] := 
  (simpln @ fixisospin[i, iz, state]) /. n -> NN+1;

Invar[{diffi_, diffsz_, diffp_}] := 
  "Invar(" <> ToString[2 diffi] <> ", " <> ToString[2 diffsz] <> ", " <> ToString[diffp] <> ");";

InvarQN[{i_, sz_, p_}] := 
  "Invar(" <> ToString[2i+1] <> ", " <> ToString[2sz] <> ", " <> ToString[p] <> ");";

rule = { Piecewise[{{value_, cond_}}, _] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr /.rule];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule,
  2Sz \[Element] Integers && Iq+Iz \[Element] Integers];

rules4 = 2Sz \[Element] Integers && 
  2Iq \[Element] Integers && Iq > 4 && Iq+Iz \[Element] Integers &&
  2Sz+Iq+Iz \[Element] Integers &&
  2Iq+NN \[Element] Integers;

(* simpl4[] is used in matrixniNEW[], thus here is where we tell
it to do (-1)^2N = 1, etc. *)
(* ratsimpl[] simplifies rational expressions. Defined in nrgcoef.m *)

simpl4[expr_] := Simplify[
ratsimpl[
  simplNN @ Simplify[expr /. rule, rules4],
  rules4]
] /. rule;

simpl5[expr_] := simpl4[expr /. rule /. {Iz->Iq}];
(* Take maximal Iz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"Iq"->"ISO(ii)"] ;
  str
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"Iq"->"ISO(iip)"] ;
  str
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"Iq"->"ISO(ii1)"];
  str 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

(* Recall: newstates2[qn,d,state] where qn are the preserved
quantum numbers, d the additional quantum numbers. see ../nrgcoef.m *)

newstates1[qn:{i_, sz_, p_}, state_] := 
  Do[newstates2[qn, {diffi, -sz, p}, state], {diffi, -i, i}];

(* Produce the expression for the new state *)

getexpr1[qn:{i_, sz_, p_}, d:{diffi_, diffsz_, diffp_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{Iq+diffi, Iz-iz}, {i, iz}, {Iq, Iz}] ],
    dirpdt[{{Iq+diffi, Sz+diffsz, P p}, {Iz-iz}, r}, 
      fixstate[i, iz, state]] ], 
    {iz, -i, i}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{i1_, sz1_, p1_}, diff1:{iz1_}, 
           qn2:{i2_, sz2_, p2_}, diff2:{iz2_},
           f[izn_, _, szn_]] := Module[{opiz, cg1, opsz, faktor},
   opiz = isoof[izn];
   cg1 = CG[{i2, iz2}, {1/2, opiz}, {i1, iz1}];
   opsz = spinofop[izn, szn];
   faktor = If[opiz == 1/2, 1, (-1)^NN (2opsz)]; (* opsz matters here! *)
   faktor cg1 ];

(* recalcf[] stuff *)

irrule[opiz_, opsz_, opp_] = {Sz -> Sz + opsz,
                              Iq -> Iq + opiz, Iz -> Iz + opiz,
                              P -> P opp};

ircg[opiz_, opsz_, opp_] :=
  CG[{Iq, Iz}, {1/2, opiz}, {Iq+opiz, Iz+opiz}] *
  If[opiz == 1/2, 1, (-1)^(NN+1) (2opsz)];

outmake[{a_, b_, c_}] :=
  "Invar(" <> toc[2(a-Iq)+ii1] <> ", "<> toc[2(b-Sz)+ssz1] <> ", " <> toc[c /. P->p1] <> ")";

diffqn[qn1:{i1_, sz1_, p1_}, qn2:{i2_, sz2_, p2_}] := Module[{diff},
  diff = {i2-i1, sz2-sz1, p1 p2};
  Simplify[diff, P^2 == 1]
];

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i1_, sz1_, p1_}, diff1:{iz1_},
            qn2:{i2_, sz2_, p2_}, diff2:{iz2_}, 
            {mi_}] := Module[{opiz, opsz},
  opiz = iz1-iz2;
  opsz = sz1-sz2;
  CG[{i2, iz2}, {mi, opiz}, {i1, iz1}] *
  If[mi == 1/2 && opiz == -1/2, 1 (2opsz), 1]
];

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai_, deltasz_, deltap_}] :=
  (Ip = Iq+deltai; Szp = Sz+deltasz; Pp = P deltap);

RULE1recalcop[{mi_}] := {Iz -> Iq};
RULE2recalcop[{mi_}] := {Iq -> Ip, Iz -> Iq-mi, Sz -> Szp, P -> Pp};

CG1recalcop[{mi_}] := 
  CG[{Ip, Iq-mi}, {mi, mi}, {Iq, Iq}] /. rule;
  (* Here no faktor is required for doublet operators,
     since mi will always be 1/2. *)

diffsD = {
  {{1/2, 1/2, 1},   FNDOUBLET <> "pp.dat"},
  {{1/2, -1/2, 1},  FNDOUBLET <> "pm.dat"},
  {{-1/2, 1/2, 1},  FNDOUBLET <> "mp.dat"},
  {{-1/2, -1/2, 1}, FNDOUBLET <> "mm.dat"}
};

diffsDodd = {
  {{1/2, 1/2, -1},   FNDOUBLET <> "diffpp.dat"},
  {{1/2, -1/2, -1},  FNDOUBLET <> "diffpm.dat"},
  {{-1/2, 1/2, -1},  FNDOUBLET <> "diffmp.dat"},
  {{-1/2, -1/2, -1}, FNDOUBLET <> "diffmm.dat"}
};

diffsS = {
  {{0, 0, 1}, FNSINGLET <> ".dat"} 
};

diffsSodd = {
  {{0, 0, -1}, FNSINGLET <> "-odd.dat"}
};

diffsT = {
  {{0, 0, 1}, FNTRIPLET <> "s.dat"}, 
  {{0, 1, 1}, FNTRIPLET <> "p.dat"},
  {{0, -1, 1}, FNTRIPLET <> "m.dat"}
};

diffsTodd = {
  {{0, 0, -1}, FNTRIPLET <> "diffs.dat"}, 
  {{0, 1, -1}, FNTRIPLET <> "diffp.dat"},
  {{0, -1, -1}, FNTRIPLET <> "diffm.dat"}
};

checksignfn[{mi_}] := ! (mi \[Element] Integers);

doall[] := Module[{},
  donew[];

  recalcop[{1/2}, diffsD];
  recalcop[{1/2}, diffsDodd];
  recalcop[{0}, diffsS];
  recalcop[{0}, diffsSodd];
  recalcop[{0}, diffsT];
  recalcop[{0}, diffsTodd];

  dooffdiagNEW[];
  (* dodiag[]; *) (* PARTICLE-HOLE SYMMETRY!! *)
  dorecalcf["ISOLR"];
];
