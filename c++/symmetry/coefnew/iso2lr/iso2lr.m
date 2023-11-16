(* 
SYMTYPE=ISO2LR  (symmetry type for an even number of impurities between conduction leads)
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2008-2023
*)

PREFIX = "iso2lr-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
basis = Simplify @ qsbasisvc @ basisops;

snegrealconstants[n, NN];

(* IMPORTANT: the following choice of nnop corresponds
to a particular partition of the lattice sites into two
sets! *)

(* ### *)  (* ### denotes places where the phase factor (-1)^n appears *)
nnop[a[]] = n;
nnop[b[]] = n+1; (* ************ ISO2LR difference.. *)
nnop[c[]] = n;
nnop[d[]] = n;
nnop[e[]] = n;

(* ### *)
Tminus = Simplify @ Total @ Map[isospinminus[#, nnop[#]]&, basisops];
Print[Tminus];

basis = bzvc2bzop @ transformQStoIS[basis];

Print[basis];

Print["lrmap=", lrmap[basisops]];

LRMAP = {a[CR, UP] -> b[CR, UP], a[CR, DO] -> b[CR,DO],
         b[CR, UP] -> a[CR, UP], b[CR, DO] -> a[CR,DO]};

Print["LRMAP=", LRMAP];

(* WORKAROUND for MMA6 and sneg < 1.192 *)
snegorthog[m_] := {} /; Tr[Abs[m]] == 0;

transformtoLRvcMY[bz_, l_, vak_, snextrarule_:{}] := Module[
  {bvc, bz2, bzsim, bzasim, bvcsim, bvcasim},
  bvc = bzop2bzvc[bz, vak];
  bz2 = bz /. Join[LRMAP, snextrarule];
  bzsim  = mapbasis[bz, bz2, 1/2 Plus[##]& ];
  bzasim = mapbasis[bz, bz2, 1/2 Subtract[##]& ];
  bvcsim  = bzop2bzvc[bzsim, vak];
  bvcasim = bzop2bzvc[bzasim, vak];
  {bvcsim, bvcasim} =
  Map[dropemptysubspaces@orthogbasisvc[#, bvc]&, {bvcsim, bvcasim}];
  bvcsim  = appendquantumnumber[bvcsim, 1];
  bvcasim = appendquantumnumber[bvcasim, -1];
  Sort @ Join[bvcsim, bvcasim]
];


(* **** CHECK THIS! How do we handle the index n when doing the LEFT <-> RIGHT transformation? *)
lrchain = basisops;
Print["lrchain=", basisops];
basis = transformtoLRvcMY[basis, lrchain, vacuum[]];

Print[basis];

basis = bzvc2bzop @ basis;
mult[{i_, s_, p_}] := (2i+1)(2s+1);
basisprep[];

(* Simplify n: simply putting n \[Element] Integers interferes with 
snegrealconstants[n] *)

simpln[expr_] := 
  Simplify[(expr /. n->nnn), nnn \[Element] Integers] /. nnn->n;

simplNN[expr_] := 
  Simplify[(expr /. NN->nnn), nnn \[Element] Integers] /. nnn->NN;

simplNN[expr_, rules_] := 
  FullSimplify[(expr /. NN->nnn), nnn \[Element] Integers && rules] /. nnn->NN;

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
  (simpln @ fixisospin[i, iz, fixspin[s, sz, state]]) /. n -> NN+1; (* ### *)

Invar[{diffi_, diffs_, diffp_}] := 
  "Invar" <> tos3[2diffi, 2diffs, diffp] <> ";";

InvarQN[{i_, s_, p_}] := 
  "Invar" <> tos3[2i+1, 2s+1, p] <> ";";

rule = { Piecewise[{{value_, cond_}}, _] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule,
  S+Sz \[Element] Integers && Iq+Iz \[Element] Integers];

(* 
NOTE: using the "trick" that consists of using a simplification rule
S > [some large value] can lead to results which involve imaginary 
coefficients, because Mathematica 6 is quite happy to simplify expression 
so that they involve square roots of negative numbers. 
Way to go, Mathematica! :-(
So we play by the book and only use S >= 0 rule.
*)

rules4 = 2S \[Element] Integers && S >= 0 &&
  2Iq \[Element] Integers && Iq >= 0 &&
  S+Sz \[Element] Integers && Iq+Iz \[Element] Integers &&
  S+Sz+Iq+Iz \[Element] Integers;
  
simpl4[expr_] := simplNN[expr /. rule, rules4] /. rule;
  
simpl5[expr_] := simpl4[expr /. {Sz->S, Iz->Iq}]; 
(* Take maximal Sz, Iz!*)

koefstr1[koef_] := Module[{str=koefstr[koef]},
  str = StringReplace[str,"S"->"S(ss)"] ; (* Careful about Sign! *)
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

newstates1[qn:{i_, s_, p_}, state_] := 
  Do[newstates2[qn, {diffi, diffs, p}, state], {diffs, -s, s}, {diffi, -i, i}];

(* Produce the expression for the new state *)

getexpr1[qn:{i_, s_, p_}, d:{diffi_, diffs_, _}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{S+diffs, Sz-sz},  {s, sz}, {S, Sz}] *
    CG[{Iq+diffi, Iz-iz}, {i, iz}, {Iq, Iz}] ],
    dirpdt[{{Iq+diffi, S+diffs, p P}, {Iz-iz, Sz-sz}, r}, 
      fixstate[i, iz, s, sz, state]] ], 
    {sz, -s, s}, {iz, -i, i}];

pbconj[expr_] := simplNN @ conj[expr];

(* NEW! *)
MATRIXNICHECKSYM[i1_, i2_, _] := 
  Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1/2 &&
  Abs[ newqn[i1] [[2]] - newqn[i2] [[2]] ] != 1/2;

(* Matrix element for the HAMILTONIAN *)
(* called from matrixel[], which computes -offdiag.dat files! *)
  
CGmatrixel[qn1:{i1_, s1_, p1_}, diff1:{iz1_, sz1_}, 
           qn2:{i2_, s2_, p2_}, diff2:{iz2_, sz2_},
           f[izn_, nr_, szn_]] := Module[{opiz, cg1, opsz, cg2, faktor, NNfact},
   opiz = isoof[izn];
   cg1 = CG[{i2, iz2}, {1/2, opiz}, {i1, iz1}];
   opsz = spinofop[izn, szn];
   cg2 = CG[{s2, sz2}, {1/2, opsz}, {s1, sz1}];
   (* ************ ISO2LR changes ************ *)
   NNfact = Switch[nr, 0, NN, 1, NN+1, _, Print["CGmatrixel error"]; Exit[]]; (* #### *)
   faktor = If[opiz == 1/2, 1, (-1)^NNfact (2opsz)];
   faktor cg1 cg2 ];

(* recalcf[] stuff *)

(* ### ??? *)
irrule[opiz_, opsz_, opp_, nn_] = {S -> S + opsz, Sz -> Sz + opsz,
                                   Iq -> Iq + opiz, Iz -> Iz + opiz,
                                   P -> P opp
};

ircg[opiz_, opsz_, opp_, nn_] := 
  CG[{S, Sz},  {1/2, opsz}, {S+opsz, Sz+opsz}] *
  CG[{Iq, Iz}, {1/2, opiz}, {Iq+opiz, Iz+opiz}] *
  If[opiz == 1/2, 1, (-1)^(NN+1+nn) (2opsz)]; (* ******** HACK!! *) (* ### *)

outmake[{a_, b_, c_}] := 
  "Invar" <> toc3[2(a-Iq)+ii1, 2(b-S)+ss1, c /. P->p1];

diffqn[qn1:{i1_, s1_, p1_}, qn2:{i2_, s2_, p2_}] := Module[{diff},
  diff = {i2-i1, s2-s1, p1 p2};
  Simplify[diff, P^2 == 1]
];

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i1_, s1_, p1_}, diff1:{iz1_, sz1_}, 
            qn2:{i2_, s2_, p2_}, diff2:{iz2_, sz2_}, 
            {mi_, ms_}] := Module[{opiz, opsz},
  opiz = iz1-iz2;
  opsz = sz1-sz2;
  CG[{i2, iz2}, {mi, opiz}, {i1, iz1}] *
  CG[{s2, sz2}, {ms, opsz}, {s1, sz1}] *
  If[mi == 1/2 && opiz == -1/2, 1 (2opsz), 1] /. rule
];

(* HMM.. In case of doublet operators, the result of
CG2recalcop[] formally depends on the index of the
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai_, deltas_, deltap_}] := 
  (Ip = Iq+deltai; Sp = S+deltas; Pp = P deltap);

RULE1recalcop[{mi_, ms_}] := {Iz -> Iq, Sz -> S};
RULE2recalcop[{mi_, ms_}] := {Iq -> Ip, Iz -> Iq-mi, S -> Sp, Sz -> S-ms, P -> Pp};

CG1recalcop[{mi_, ms_}] := 
  CG[{Ip, Iq-mi}, {mi, mi}, {Iq, Iq}] *
  CG[{Sp, S-ms},  {ms, ms}, {S, S}] /. rule;
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

checksignfn[{mi_, ms_}] := ! (ms \[Element] Integers);

dorecalcfLOOP["ISO2LR", p_] := Module[{str, ch, op, nn},
  str = If[p == 1, "", "diff"];
  For[ch = 0, ch <= 1, ch++,
    op = ch2op[ch];
    nn = nnop[op[]]-n; (* ********** JUST THE SHIFT! *)
    Print["op=", op, " nn=", nn];
    makerecalcf[op[CR, UP], {1/2, 1/2, p, nn},
    FNSPINUPISOUP <> str <> tos[op]];
    makerecalcf[op[CR, DO], {1/2, -1/2, p, nn},
    FNSPINDOWNISOUP <> str <> tos[op] ];
    makerecalcf[op[AN, DO], {-1/2, 1/2, p, nn},
    FNSPINUPISODOWN <> str <> tos[op]];
    makerecalcf[op[AN, UP], {-1/2, -1/2, p, nn},
    FNSPINDOWNISODOWN <> str <> tos[op] ];
  ];
];

dorecalcfLOOP["ISO2LR"] := Module[{},
  dorecalcfLOOP["ISO2LR", 1];
  dorecalcfLOOP["ISO2LR", -1];
];

doall[] := Module[{},
  donew[];
  dooffdiagNEW[];
  dorecalcf["ISO2LR"];

  recalcop[{1/2, 1/2}, diffsD];
  recalcop[{1/2, 1/2}, diffsDodd];
  recalcop[{0, 0}, diffsS];
  recalcop[{0, 0}, diffsSodd];
  recalcop[{0, 1}, diffsT];
  recalcop[{0, 1}, diffsTodd];
];
