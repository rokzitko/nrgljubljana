(* Helper functions for SU(2)xSO(3) basis generation for three orbital
   problems. Rok Zitko, Aug 2015 *)

(* Total spin, total angular momentum *)
oS = 1/2;
oT = 1;

(* Operators for oT=1, oS=1/2 will be of the form d[CR,t,sz], where
   t={-1,0,1}, sz={UP,DO}. *)

opSM[op_[q___]] := Sum[spinminus[op[q, tz]], {tz, -oT, oT}];
opSP[op_[q___]] := Sum[spinplus[op[q, tz]], {tz, -oT, oT}];
opSX[op_[q___]] := Sum[spinx[op[q, tz]], {tz, -oT, oT}];
opSY[op_[q___]] := Sum[spiny[op[q, tz]], {tz, -oT, oT}]; 
opSZ[op_[q___]] := Sum[spinz[op[q, tz]], {tz, -oT, oT}];

(* Angular momentum operators for sites listes in list ops. *)

msz[1/2] = UP;
msz[-1/2] = DO;

tminus[op_[q___, sz_]] := 
 VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixM[oT], 
 Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];
tplus[op_[q___, sz_]] := 
 VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixP[oT], 
 Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];
txx[op_[q___, sz_]] := 
 VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixX[oT], 
 Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];
tyy[op_[q___, sz_]] := 
 VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixY[oT], 
 Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];
tzz[op_[q___, sz_]] := 
 VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixZ[oT], 
 Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];

opTM[op_[q___]] := Sum[tminus[op[q, msz[sz]]], {sz, -oS, oS}];
opTP[op_[q___]] := Sum[tplus[op[q, msz[sz]]], {sz, -oS, oS}];
opTX[op_[q___]] := Sum[txx[op[q, msz[sz]]], {sz, -oS, oS}]; 
opTY[op_[q___]] := Sum[tyy[op[q, msz[sz]]], {sz, -oS, oS}]; 
opTZ[op_[q___]] := Sum[tzz[op[q, msz[sz]]], {sz, -oS, oS}];
 
preponesite[] := Module[{},
 ops = {d[]};

 SM = Total @ Map[opSM, ops];
 SP = Total @ Map[opSP, ops];
 SX = Total @ Map[opSX, ops];
 SY = Total @ Map[opSY, ops];
 SZ = Total @ Map[opSZ, ops];
 SS = Expand[nc[SZ, SZ] + nc[SY, SY] + nc[SX, SX]];

 TM = Total@Map[opTM, ops];
 TP = Total@Map[opTP, ops];
 TX = Total@Map[opTX, ops];
 TY = Total@Map[opTY, ops];
 TZ = Total@Map[opTZ, ops];
 TT = Expand[nc[TZ, TZ] + nc[TY, TY] + nc[TX, TX]];
];
 
calcS2[op_] := vev[conj[op]~nc~zeroonvac[SS~nc~op]];

calcT2[op_] := vev[conj[op]~nc~zeroonvac[TT~nc~op]];

(* Resum lists of the form {{a,x},{b,x}...} into (a+b+..)x. *)
resum[l : {{_, a_} ..}] := Total[l[[All, 1]]] a;

(* Separate coefficient and operator parts. *)
(* ASSUMPTION: z nc form *or* op form. *)
ClearAll[splitcoef];
splitcoef[HoldPattern[z_. y : nc[_[___] ..]]] := {z, y};
splitcoef[z_. y : _?operatorQ[___]] := {z, y};
splitcoef::form = "Error in splitcoef: 1";
splitcoef[x_] := Message[splitcoef::form, x];

(* Factorize all coefficients in an expression involving operators. *)
(* Assumption: a linear combination of (products of) particle 
creation operators *)
ClearAll[factorizecoefs];
factorizecoefs[HoldPattern[Plus[l__]]] := Module[{p},
 p = Map[splitcoef, {l}];
 p = GatherBy[p, #[[2]] &];
 p = Map[resum, p];
 Total[p]
];
factorizecoefs[x_] := x;

(* SLOW *)
normalizeopOLD[op_] := Simplify[op/Simplify[normop[op]]];

(* POZOR: "ohranjati" je potrebno predznak!! *)

ClearAll[normalizeopNEW, normalizeop1];
normalizeop1[z_?NumericQ] := z/Abs[z];
normalizeop1[z_?NumericQ x_] := (z/Abs[z]) normalizeop1[x];
normalizeop1[x : _?operatorQ[___]] := x;
normalizeop1[HoldPattern[x : nc[_[___] ..]]] := x;
normalizeop1[HoldPattern[x : Plus[l__]]] := Simplify[x/Simplify[Norm[
 {l} /. {HoldPattern[nc[_[___] ..]] -> 1, 
 o_?operatorQ[___] -> 1}
]]];
normalizeopNEW[x_] := normalizeop1[factorizecoefs[x]];

ClearAll[normalizeop];
(* Caching! *)
(* Expand[] is necessary in general! *)
normalizeop[x_] := normalizeop[x] = normalizeopNEW[Expand[x]];

ClearAll[normopNEW, normop1];
normop1[z_?NumericQ] := Abs[z];
normop1[z_?NumericQ x_] := Abs[z] normop1[x];
normop1[x : _?operatorQ[___]] := 1;
normop1[HoldPattern[x : nc[_[___] ..]]] := 1;
normop1[HoldPattern[x : Plus[l__]]] := Simplify[Norm[{l} /. 
 {HoldPattern[nc[_[___] ..]] -> 1, o_?operatorQ[___] -> 1}
]];
normopNEW[x_] := normop1[factorizecoefs[Expand[x]]];

(* SLOW *)
sdownOLD[op_] := Module[{b},
 b = Expand@zeroonvac[nc[SM, op]];
 b = normalizeop[b];
 b = Expand[b]
];
tdownOLD[op_] := Module[{b},
 b = Expand @ zeroonvac[nc[TM, op]];
 b = normalizeop[b];
 b = Expand[b]
];

ClearAll[sdown, tdown];
(* Caching *)
sdown[op_] := sdown[op] = Module[{b, sm},
 sm = Plus @@ Map[spinminus, Union[Cases[op, 
  x_?operatorQ[_, j___, lz_, sz_] :> x[j, lz], {0, Infinity}]]];
 b = Expand@zeroonvac[nc[sm, op]];
 b = normalizeop[b];
 b = Expand[b]
];
tdown[op_] := tdown[op] = Module[{b, tm},
 tm = Plus @@ Map[tminus, Union[Cases[op, 
  x_?operatorQ[_, j___, lz_, sz_] :> x[j, sz], {0, Infinity}]]];
 b = Expand @ zeroonvac[nc[tm, op]];
 b = normalizeop[b];
 b = Expand[b]
];

(* Given maximal-weight state op (with sz=s and tz=t), generate a 
corresponding state with chosen sz<s and tz<t. *)
(* IMPORTANT: the output is normalized to 1, but we follow a determinate 
phase convention, thus normalizations need to conserve the overall sign!!! *)

ClearAll[fixspinspin];
(* Caching *)
fixspinspin[op_, s_, sz_, t_, tz_] := fixspinspin[op, s, sz, t, tz] = 
  Nest[sdown, Nest[tdown, op, t - tz], s - sz];

Off[ClebschGordan::phy];
(* Combine two states using angular momentum algebra for spin SU(2) 
   and angular momentum SO(3). It's a simple double sum over Sz and Tz 
   with corresponding Clebsch-Gordan coefficients. *)
  
SU2SU2combine[op1_, op2_, {q1 : {s1_, t1_}, q2 : {s2_, t2_}}, 
 fixqn1fnc_, fixqn2fnc_] := Module[{l},
  VPrintTemporary[op1, " # ", op2];
  l = Table[{{s, t}, Expand @ Sum[
   Module[{cg1, cg2},
    Off[ClebschGordan::phy];
    cg1 = ClebschGordan[{s1, sz1}, {s2, sz2}, {s, s}];
    cg2 = ClebschGordan[{t1, tz1}, {t2, tz2}, {t, t}];
    If[cg1 cg2 =!= 0,
     Simplify[ cg1 cg2 nc[fixqn1fnc[op1, s1, sz1, t1, tz1], 
                          fixqn2fnc[op2, s2, sz2, t2, tz2]]  ],
     0]
    ],
  {sz1, -s1, s1}, {sz2, -s2, s2}, {tz1, -t1, t1}, {tz2, -t2, t2}]},
  {s, Abs[s1 - s2], s1 + s2}, {t, Abs[t1 - t2], t1 + t2}];
  l = Flatten[l, 1];
  (* Drop non-existing states, arising from "destructive interference". 
     That's why the signs matter!! *)
  l = DeleteCases[l, {_, 0}];
  (* Normalize the resulting states. *)
  l = Map[{First[#], normalizeop@Last[#]} &, l]
];
  
sameuptosign[a_, b_] := (Simplify[a - b] === 0) || (Simplify[a + b] === 0);

(* Combine two sets of states using SU2SU2combine. *)

SU2SU2merge[{Qn1_, l1_}, {Qn2_, l2_}, qn1fnc_, qn2fnc_, combineqnfnc_,
 fixqn1fnc_, fixqn2fnc_] := Module[{qn1, qn2, states},
 VPrintTemporary["***", Qn1, Qn2];
 (* {s,t} for first set *)
 qn1 = qn1fnc[Qn1];
 (* {s,t} for second set *)
 qn2 = qn2fnc[Qn2];
 states = Flatten[Outer[
  SU2SU2combine[#1, #2, {qn1, qn2}, fixqn1fnc, fixqn2fnc] &,
  l1, l2, 1], 2];
 states = mergebasis@
  Map[{combineqnfnc[First[#], Qn1, Qn2], {Last[#]}} &, states];
 (* Drop non-existing states. *)
 states = Map[{First[#], DeleteCases[Last[#], 0]} &, states];
 (* Drop redundant states. 
    The equivalence relation is defined up to overall sign!! *)
 states = 
  Map[{First[#], Union[Last[#], SameTest -> sameuptosign]} &, states];
 (* Drop empty subspaces. *)
 states = Select[states, Last[#] =!= {} &]
];

(* Quantum number addition: we need to add particle number q. *)
combineaddq[{qns__}, {q1_, ___}, {q2_, ___}] := {q1 + q2, qns};

getst[{q_, s_, t_}] := {s, t};

(* Tensor product between two basis sets. *)
SU2SU2basistensorproduct[l1_, l2_, other___] := Module[{l},
 l = Outer[SU2SU2merge[#1, #2, other] &, l1, l2, 1];
 l = Flatten[l, 2];
 l = mergebasis[l];
 (* Drop duplicate states. *)
 l = Map[{First[#], Union[Last[#], SameTest -> sameuptosign]} &, l]; (* UP TO SIGN! *)
 l
];
 
preptwosites[] := Module[{},
 ops = {d[], f[]};

 SM = Total @ Map[opSM, ops];
 SP = Total @ Map[opSP, ops];
 SX = Total @ Map[opSX, ops];
 SY = Total @ Map[opSY, ops];
 SZ = Total @ Map[opSZ, ops];
 SS = Expand[nc[SZ, SZ] + nc[SY, SY] + nc[SX, SX]];

 TM = Total@Map[opTM, ops];
 TP = Total@Map[opTP, ops];
 TX = Total@Map[opTX, ops];
 TY = Total@Map[opTY, ops];
 TZ = Total@Map[opTZ, ops];
 TT = Expand[nc[TZ, TZ] + nc[TY, TY] + nc[TX, TX]];
];
 
