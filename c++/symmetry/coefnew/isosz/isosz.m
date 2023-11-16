(*
SYMTYPE=ISOSZ
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2009-2023
*)

PREFIX = "isosz-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
basis = Simplify @ qszbasisvc @ basisops;

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
mult[{i_, _}] := 2i+1;
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

Invar[{diffi_, diffsz_}] := 
  "Invar(" <> ToString[2 diffi] <> ", " <> ToString[2 diffsz] <> ");";

InvarQN[{i_, sz_}] := 
  "Invar(" <> ToString[2i+1] <> ", " <> ToString[2sz] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ Simplify[expr, 2S \[Element] Integers && S >= 3];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule,
  2Sz \[Element] Integers && Iq+Iz \[Element] Integers && NN \[Element] Integers];

rules4 = 2Sz \[Element] Integers && 
  2Iq \[Element] Integers && Iq > 4 && Iq+Iz \[Element] Integers &&
  2Sz+Iq+Iz \[Element] Integers;

simpl4[expr_] := Simplify[ (ratsimpl[ simplNN @ Simplify[expr, rules4], rules4]) /. rule ];

simpl5[expr_] := simpl4[expr /. {Iz->Iq}]; 
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

newstates1[qn:{i_, sz_}, state_] := 
  Do[newstates2[qn, {diffi, -sz}, state], {diffi, -i, i}];

(* Produce the expression for the new state *)

getexpr1[qn:{i_, sz_}, d:{diffi_, diffsz_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{Iq+diffi, Iz-iz}, {i, iz}, {Iq, Iz}] ],
    dirpdt[{{Iq+diffi, Sz+diffsz}, {Iz-iz}, r}, 
      fixstate[i, iz, state]] ], 
    {iz, -i, i}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{i1_, sz1_}, diff1:{iz1_}, 
           qn2:{i2_, sz2_}, diff2:{iz2_},
           f[izn_, _, szn_]] := Module[{opiz, cg1, opsz, faktor},
   opiz = isoof[izn];
   cg1 = CG[{i2, iz2}, {1/2, opiz}, {i1, iz1}];
   opsz = spinofop[izn, szn];
   faktor = If[opiz == 1/2, 1, (-1)^NN (2opsz)]; (* opsz matters here! *)
   faktor cg1 ];

(* recalcf[] stuff *)

irrule[opiz_, opsz_] = {Sz -> Sz + opsz,
                        Iq -> Iq + opiz, Iz -> Iz + opiz};

ircg[opiz_, opsz_] := 
  CG[{Iq, Iz}, {1/2, opiz}, {Iq+opiz, Iz+opiz}] *
  If[opiz == 1/2, 1, (-1)^(NN+1) (2opsz)];

outmake[{a_, b_}] := 
  "Invar(" <> toc[2(a-Iq)+ii1] <> ", "<> toc[2(b-Sz)+ssz1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i1_, sz1_}, diff1:{iz1_}, 
            qn2:{i2_, sz2_}, diff2:{iz2_}, 
            {mi_, ms_}] := Module[{opiz, opsz},
  opiz = iz1-iz2;
  opsz = sz1-sz2;
  CG[{i2, iz2}, {mi, opiz}, {i1, iz1}] *
  If[mi == 1/2 && opiz == -1/2, 1 (2opsz), 1]
];

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai_, deltasz_}] :=
  (Ip = Iq+deltai; Szp = Sz+deltasz);

RULE1recalcop[{mi_, ms_}] := {Iz -> Iq};
RULE2recalcop[{mi_, ms_}] := {Iq -> Ip, Iz -> Iq-mi, Sz -> Szp};

CG1recalcop[{mi_, ms_}] := 
  CG[{Ip, Iq-mi}, {mi, mi}, {Iq, Iq}];
  (* Here no faktor is required for doublet operators,
     since mi will always be 1/2. *)

diffsD = {
  {{1/2, 1/2},   FNDOUBLET <> "pp.dat"},
  {{1/2, -1/2},  FNDOUBLET <> "pm.dat"},
  {{-1/2, 1/2},  FNDOUBLET <> "mp.dat"},
  {{-1/2, -1/2}, FNDOUBLET <> "mm.dat"}
};

diffsS = { 
  {{0, 0}, FNSINGLET <> ".dat"} 
};

diffsT = { 
  {{0, 0}, FNTRIPLET <> "s.dat"}, 
  {{0, 1}, FNTRIPLET <> "p.dat"},
  {{0, -1}, FNTRIPLET <> "m.dat"}
};

checksignfn[{mi_, ms_}] := ! (ms \[Element] Integers);

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

doall[] := Module[{},
  donew[];

  (* see ../qsz/qsz.m *)
  (* ATTENTION: upper=False here! These are not used to construct the Hamiltonian... *)
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];

  dooffdiagNEW[];
  (* dodiag[]; *) (* PARTICLE-HOLE SYMMETRY!! *)

  dorecalcf["ISO"]; (* Note: for QSZ we use QS (i.e. the default). *)

  recalcop[{1/2, 1/2}, diffsD];
  recalcop[{0, 0}, diffsS];
  recalcop[{0, 1}, diffsT];
];
