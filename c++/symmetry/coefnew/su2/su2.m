(*
SYMTYPE=SU2
Part of "NRG Ljubljana"
Rok Zitko, rok.zitko@ijs.si, 2009-2023
$Id$ 
*)

PREFIX = "su2-" <> ToString[channels] <> "ch";
fnprep[];

(* Make basis *)
basisops = Take[allops, channels];
basis = Simplify @ qbasisvc @ basisops;

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
mult[{i_}] := 2i+1;
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

Invar[{diffi_}] := 
  "Invar(" <> ToString[2 diffi] <> ");";

InvarQN[{i_}] := 
  "Invar(" <> ToString[2i+1] <> ");";

rule = { Piecewise[{{value_, cond_}}, 0] :> value }; (* new 2023 *)

simplifynew[expr_] := simpl2 @ FullSimplify[expr];
SIMPLIFYNEW = True;

simpl2[expr_] := Simplify[expr /. rule, Iq+Iz \[Element] Integers && NN \[Element] Integers];

rules4 = 2Iq \[Element] Integers && Iq > 4 && Iq+Iz \[Element] Integers;

simpl4[expr_] := ratsimpl[ simplNN @ FullSimplify[expr, rules4], rules4];

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

newstates1[qn:{i_}, state_] := 
  Do[newstates2[qn, {diffi}, state], {diffi, -i, i}];

(* Produce the expression for the new state *)

getexpr1[qn:{i_}, d:{diffi_}, state_] := 
  Sum[ nc[coef @ simpl4[ 
    CG[{Iq+diffi, Iz-iz}, {i, iz}, {Iq, Iz}] ],
    dirpdt[{{Iq+diffi}, {Iz-iz}, r}, fixstate[i, iz, state]] ], 
  {iz, -i, i}];

pbconj[expr_] := simplNN @ conj[expr];

(* Matrix element for the HAMILTONIAN *)
CGmatrixel[qn1:{i1_}, diff1:{iz1_}, 
           qn2:{i2_}, diff2:{iz2_},
           f[izn_, _, szn_]] := Module[{opiz, cg1, opsz, faktor},
   opiz = isoof[izn];
   cg1 = CG[{i2, iz2}, {1/2, opiz}, {i1, iz1}];
   opsz = spinofop[izn, szn];
   faktor = If[opiz == 1/2, 1, (-1)^NN (2opsz)]; (* !!! Restored !!! *)
   faktor cg1 ];

(* irrule[] and ircg[] are called from recalcf[]. The results are then used
as parameters to ireducel[] in nrgcoef.m. irrule is used to prepare the 
bra when computing brakets. ircg is used as the denominator in the resulting
expression. *)

irrule[opiz_] = {Iq -> Iq + opiz, Iz -> Iz + opiz};

ircg[opiz_] := 
  CG[{Iq, Iz}, {1/2, opiz}, {Iq+opiz, Iz+opiz}] *
  If[opiz == 1/2, 1, (-1)^(NN+1)];

outmake[{a_}] := 
  "Invar(" <> toc[2(a-Iq)+ii1] <> ")";

diffqn[qn1:{__}, qn2:{__}] := qn2-qn1;

(* TRICK: just take the difference, iz1-iz2, etc. !!! *)
CG2recalcop[qn1:{i1_}, diff1:{iz1_}, 
            qn2:{i2_}, diff2:{iz2_}, 
            {mi_}] := Module[{opiz},
  opiz = iz1-iz2;
  CG[{i2, iz2}, {mi, opiz}, {i1, iz1}] 
];

(* HMM.. In case of doublet operators, the result of 
CG2recalcop[] formally depends on the index of the 
lattice site. Is this a problem in practice? *)

QNrecalcop[delta:{deltai__}] :=
  (Ip = Iq+deltai);

RULE1recalcop[{mi_}] := {Iz -> Iq};
RULE2recalcop[{mi_}] := {Iq -> Ip, Iz -> Iq-mi};

CG1recalcop[{mi_}] := 
  CG[{Ip, Iq-mi}, {mi, mi}, {Iq, Iq}];
  (* Here no faktor is required for doublet operators,
     since mi will always be 1/2. *)

(* NEW: Return True if the (i1, i2) OFFDIAG line is zero for purely
   symmetry reasons. *)
MATRIXNICHECKSYM[i1_, i2_, _] := Abs[ newqn[i1] [[1]] - newqn[i2] [[1]] ] != 1/2;

diffsD = {
  {{1/2},   FNDOUBLET <> "p.dat"},
  {{-1/2},  FNDOUBLET <> "m.dat"}
};

diffsS = {
  {{0}, FNSINGLET <> ".dat"} 
};

checksignfn[{mi_}] := ! (mi \[Element] Integers);

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},
  op = ch2op[ch];
  opspinz = spinz[ op[] ];
  matrixop1[i1, i2, opspinz]
];

matrixnihop1[i1_, i2_, ch_] := Module[{op, expr = 0},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR, ch, UP], op[AN, UP]]];
  expr += matrixel[i1, i2, nc[op[CR, DO], f[AN, ch, DO]]];
  simpl4 @ expr  
];

matrixnihop2[i1_, i2_, ch_] := Module[{op, expr = 0},
  If[MATRIXNICHECKSYM[i1, i2, ch], Print["skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR, ch, DO], op[AN, DO]]];
  expr += matrixel[i1, i2, nc[op[CR, UP], f[AN, ch, UP]]];
  simpl4 @ expr  
];

dorecalcfLOOP["SU2"] := Module[{ch, op, nn},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];

    (* makerecalcf[] -- First argument: operator. Second argument: izop. *)

    (* First doublet: [f^dag_UP, f_DO] szop=1/2 pair *)
    makerecalcf[op[CR, UP], {1/2},  PREFIX <> "-type1-isoup-" <> tos[op]];
    makerecalcf[op[AN, DO], {-1/2}, PREFIX <> "-type1-isodown-" <> tos[op]];

    (* Second doublet: [f^dag_DO, -f_UP] szop=-1/2 pair *)
    makerecalcf[op[CR, DO], {1/2},  PREFIX <> "-type2-isoup-" <> tos[op]];
    makerecalcf[-op[AN, UP], {-1/2}, PREFIX <> "-type2-isodown-" <> tos[op]];
  ];
];

matrixopspinz[i1_, i2_, ch_] := Module[{op, opspinz},  
 op = ch2op[ch];
 opspinz = spinz[ op[] ];
 matrixop1[i1, i2, opspinz]
];

matrixopspiny[i1_, i2_, ch_] := Module[{op, opspiny},
 op = ch2op[ch];
 opspiny = spiny[ op[] ];
 matrixop1[i1, i2, opspiny]
];

matrixopspinx[i1_, i2_, ch_] := Module[{op, opspinx},
 op = ch2op[ch];
 opspinx = spinx[ op[] ];
 matrixop1[i1, i2, opspinx]
];

doall[] := Module[{},
  donew[];
  (* ATTENTION: upper=False here! These are not used to construct the Hamiltonian... *)
  dogeneralloop[matrixopspinx, "SPINX", PREFIX <> "-spinx.dat", False];
  dogeneralloop[matrixopspiny, "SPINY", PREFIX <> "-spiny.dat", False];
  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-spinz.dat", False];

  dogeneralloop[matrixopspinz, "SPINZ", PREFIX <> "-H-spinz.dat", True];

  dogeneralloop[matrixnihop1, "OFFDIAG_1", PREFIX <> "-offdiag-1.dat", True];
  dogeneralloop[matrixnihop2, "OFFDIAG_2", PREFIX <> "-offdiag-2.dat", True];

  dorecalcf["SU2"];
  recalcop[{1/2}, diffsD];
  recalcop[{0}, diffsS];
];
