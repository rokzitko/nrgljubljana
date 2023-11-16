(*
SYMTYPE=ISO
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "iso-" <> ToString[channels] <> "ch";
fnprep[];

debug = 3;

(* Make basis *)
basisops = Take[allops, channels];
basis = Simplify @ qsbasisvc @ basisops;

snegrealconstants[n, NN];

(* IMPORTANT: the following choice of nnop corresponds
to a particular partition of the lattice sites into two
sets! *)

nnop[a[]] = n;
nnop[b[]] = n;
nnop[c[]] = n;
nnop[d[]] = n;
nnop[e[]] = n;

Tminus = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops];
Print[Tminus];

basis = bzvc2bzop @ transformQStoIS[basis];
mult[{i_, s_}] := (2i+1)(2s+1);
basisprep[];

(* Simplify n: simply putting n \[Element] Integers interferes with 
snegrealconstants[n] *)

simpln[expr_] :=
  Simplify[(expr /. n->nnn), nnn \[Element] Integers] /. nnn->n;

simplNN[expr_] :=
  (Simplify[(expr /. NN->nnn), nnn \[Element] Integers] /. nnn->NN) /. {Exp[I Pi NN]->(-1)^NN, Exp[-I Pi NN]->(-1)^NN};

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

fixstate[i_, iz_, s_, sz_, state_] :=
  (simpln @ fixisospin[i, iz, fixspin[s, sz, state]]) /. n -> NN+1;

Invar[{diffi_, diffs_}] :=
  "Invar(" <> ToString[2 diffi] <> ", " <> ToString[2 diffs] <> ");";

InvarQN[{i_, s_}] :=
  "Invar(" <> ToString[2i+1] <> ", " <> ToString[2s+1] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ Simplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr //. rule,
  S+Sz \[Element] Integers && Iq+Iz \[Element] Integers] //. rule;

rules4 = 
  2S \[Element] Integers && S > 4 && 
  2Iq \[Element] Integers && Iq > 4 &&
  S+Sz \[Element] Integers && Iq+Iz \[Element] Integers &&
  S+Sz+Iq+Iz \[Element] Integers;

simpl4[expr_] := simplNN[ ratsimpl[ simplNN @ Simplify[expr /. rule, rules4], rules4] /. rule ];

simpl5[expr_] := simpl4[expr /. {Sz->S, Iz->Iq}] /. rule;
(* Take maximal Sz, Iz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss)"] ;
  str = StringReplace[str,"Iq"->"ISO(ii)"] ;
  str
];

koefstr2[koef_] := Module[{str=koefstr[koef]},
  (* recalcf[] *)
  str = StringReplace[str,"S"->"S(ssp)"] ;
  str = StringReplace[str,"Iq"->"ISO(iip)"] ;
  str
];

koefstr3[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss1)"];
  str = StringReplace[str,"Iq"->"ISO(ii1)"];
  str 
];

(* Loop over all possible ways of combining states (cf. angular 
momentum addition rules, etc. *)

newstates1[qn:{i_, s_}, state_] := 
  Do[newstates2[qn, {diffi, diffs}, state], {diffs, -s, s}, {diffi, -i, i}];

(* Produce the expression for the new state *)

getexpr1[qn:{i_, s_}, d:{diffi_, diffs_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{S+diffs, Sz-sz},  {s, sz}, {S, Sz}] *
    CG[{Iq+diffi, Iz-iz}, {i, iz}, {Iq, Iz}] ],
    dirpdt[{{Iq+diffi, S+diffs}, {Iz-iz, Sz-sz}, r}, 
      fixstate[i, iz, s, sz, state]] ], 
    {sz, -s, s}, {iz, -i, i}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{i1_, s1_}, diff1:{iz1_, sz1_}, 
           qn2:{i2_, s2_}, diff2:{iz2_, sz2_},
           f[izn_, _, szn_]] := Module[{opiz, cg1, opsz, cg2, faktor},
   opiz = isoof[izn];
   cg1 = CG[{i2, iz2}, {1/2, opiz}, {i1, iz1}];
   opsz = spinofop[izn, szn];
   cg2 = CG[{s2, sz2}, {1/2, opsz}, {s1, sz1}];
   faktor = If[opiz == 1/2, 1, (-1)^NN (2opsz)];
   faktor cg1 cg2 ];

(* recalcf[] stuff *)

irrule[opiz_, opsz_] = {S -> S + opsz, Sz -> Sz + opsz,
                        Iq -> Iq + opiz, Iz -> Iz + opiz};

ircg[opiz_, opsz_] :=
  CG[{S, Sz},  {1/2, opsz}, {S+opsz, Sz+opsz}] *
  CG[{Iq, Iz}, {1/2, opiz}, {Iq+opiz, Iz+opiz}] *
  If[opiz == 1/2, 1, (-1)^(NN+1) (2opsz)];

outmake[{a_, b_}] :=
  "Invar(" <> toc[2(a-Iq)+ii1] <> ", "<> toc[2(b-S)+ss1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i1_, s1_}, diff1:{iz1_, sz1_},
            qn2:{i2_, s2_}, diff2:{iz2_, sz2_},
            {mi_, ms_}] := Module[{opiz, opsz},
  opiz = iz1-iz2;
  opsz = sz1-sz2;
  CG[{i2, iz2}, {mi, opiz}, {i1, iz1}] *
  CG[{s2, sz2}, {ms, opsz}, {s1, sz1}] *
  If[mi == 1/2 && opiz == -1/2, 1 (2opsz), 1]
] //. rule;

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai_, deltas_}] :=
  (Ip = Iq+deltai; Sp = S+deltas);

RULE1recalcop[{mi_, ms_}] := {Iz -> Iq, Sz -> S};
RULE2recalcop[{mi_, ms_}] := {Iq -> Ip, Iz -> Iq-mi, S -> Sp, Sz -> S-ms};

CG1recalcop[{mi_, ms_}] :=
  CG[{Ip, Iq-mi}, {mi, mi}, {Iq, Iq}] *
  CG[{Sp, S-ms},  {ms, ms}, {S, S}] //. rule;
  (* Here no faktor is required for doublet operators,
     since mi will always be 1/2. *)

diffsS = {
  {{0, 0}, FNSINGLET <> ".dat"}
};

diffsD = {
  {{1/2, 1/2},   FNDOUBLET <> "pp.dat"},
  {{1/2, -1/2},  FNDOUBLET <> "pm.dat"},
  {{-1/2, 1/2},  FNDOUBLET <> "mp.dat"},
  {{-1/2, -1/2}, FNDOUBLET <> "mm.dat"}
};

diffsT = {
  {{0, 0}, FNTRIPLET <> "s.dat"},
  {{0, 1}, FNTRIPLET <> "p.dat"},
  {{0, -1}, FNTRIPLET <> "m.dat"}
};

checksignfn[{mi_, ms_}] := ! (ms \[Element] Integers);

doall[] := Module[{},
  donew[];
  dooffdiagNEW[];
  dorecalcf["ISO"];

  recalcop[{0, 0}, diffsS];
  recalcop[{1/2, 1/2}, diffsD];
  recalcop[{0, 1}, diffsT];
];
