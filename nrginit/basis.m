(****************************** BASIS STATES ************************************)

(* We define a class Invar[] which wraps the quantum numbers, i.e.
   Invar[q,ss] for (Q,S) NRG or Invar[ii,ss] for (I,S) NRG, etc. For
   backward compatibility, we also allow simple List[] wrapping. *)

makebasis[basisops]; (* Occupation number representation *)
vak = vacuum[]; (* vc[0, 0, ..., 0] *)

SS2S[ss_] = (ss-1)/2;
II2I[ii_] = (ii-1)/2;
JJ2J[jj_] = (jj-1)/2;

basisfilename = "basis";
hamfilename = "ham";
opfilename =  "op";

writedir = paramdefault["writedir", ""];
MyPut[x_, fn_, forcewrite_:False] := If[option["WRITE"] || forcewrite, Put[x, writedir <> fn]];

GENERATEBASIS = If[option["READBASIS"], False, True];

If[GENERATEBASIS == False,

  (* NEW: If requested (i.e. if there is a line "options=READBASIS" in the
  parameters file), we read the base states from a file. Thus repeated
  generation of the Hamiltonian matrices is avoided which can save a
  considerable amount of time if parameter sweeps are performed.  Warning:
  no checks are performed at this time, so you better know what you are
  doing. *)

  MyPrint["Reading basis from " <> basisfilename];
  bvc = silentGet[basisfilename, Path -> dumppath];
  MyVPrint[2, "bvc=", bvc];

  (* If reading fails, fall back to generating the states. *)
  If[bvc === $Failed,
    GENERATEBASIS = True;
  ];
];


If[GENERATEBASIS == True,
  MyPrint["Generating basis"];
  timestart["basis"];

  (*** STEP 1: create basis operators ***)

  needqslist = {"QS", "QSLR", "ISO", "ISO2", "ISOLR", "ISO2LR"};
  needqszlist = {"QSZ", "QSZLR", "ISOSZ", "ISOSZLR"};
  needqlist = {"SU2", "U1", "SL"};

  If[MemberQ[needqslist, SYMTYPE],
    (* Returns creation operators for basis states in (Q,S) basis. In
    the case of ISO, this is transformed in step 2 to (I,S) basis. *)

    bz = bzvc2bzop[ qsbasisvc[basisops] ];
  ];

  If[MemberQ[needqszlist, SYMTYPE],
    bz = qszbasis[basisops];
  ];

  If[isP[],
    bz = qbasis[basisops];
    bz = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, bz]; (* Charge parity *)
    bz = mergebasis[bz];
  ];

  If[isPP[],
    If[CHANNELS != 2, MyError["SYMTYPE=PP requires two channels!"]];
    If[Length[basisops] != 2, MyError["General case not implemented yet."]];
    basis1 = qbasis[ {basisops[[1]]} ];
    basis2 = qbasis[ {basisops[[2]]} ];

    basis1 = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, basis1]; (* Charge parity *)
    basis1 = mergebasis[basis1];

    basis2 = Map[ {2(Mod[First[#], 2]-1/2), Last[#]}&, basis2]; (* Charge parity *)
    basis2 = mergebasis[basis2];

    makebasis[ basisops ]; (* !! *)

    bz = basistensorproduct[basis1, basis2, Function[{qn1,qn2}, {qn1[[1]],qn2[[1]]}]];
  ];

  If[isNONE[],
    bz = nonebasis[basisops];
    bz = Map[{{0}, #[[2]]}&, bz];
  ];

  If[MemberQ[needqlist, SYMTYPE],
    bz = qbasis[basisops];
  ];

  If[isSPSU2[] || isSPSU2LR[],
    bz = sbasis[basisops];
  ];

  If[isSPU1[] || isSPU1LR[],
    bz = szbasis[basisops];
  ];

  If[isDBLQSZ[],
    MyAssert[CHANNELS == 2];
    If[Length[basisops] == 2,
      basisops1 = {First[basisops]};
      basisops2 = {Last[basisops]},
    (* else *)
      MyAssert[Mod[Length[lrchain],2] == 0];
      basisops1 = Take[lrchain, Length[lrchain]/2];
      basisops2 = Take[lrchain, -Length[lrchain]/2];
    ];
    MyPrint["basisops1=", basisops1];
    MyPrint["basisops2=", basisops2];
    bz = quickDBLSZ[basisops1, basisops2, qszbasis];
  ];

  If[isDBLSU2[],
    MyAssert[CHANNELS == 2];
    MyAssert[Length[basisops] == 2];
    basisops1 = {First[basisops]};
    basisops2 = {Last[basisops]};
    bz = quickDBL[basisops1, basisops2, quickSU2basis];
  ];

  (* Important: we assume that lrchain is properly defined! *)
  If[isDBLISOSZ[],
    MyAssert[CHANNELS == 2];
    If[Length[basisops] == 2,
      basisops1 = {First[basisops]};
      basisops2 = {Last[basisops]},
    (* else *)
      MyAssert[Mod[Length[lrchain],2] == 0];
      basisops1 = Take[lrchain, Length[lrchain]/2];
      basisops2 = Take[lrchain, -Length[lrchain]/2];
    ];
    bz = quickDBLSZ[basisops1, basisops2, quickISOSZbasis];
  ];

  If[isSL3[],
    Module[{ops1,ops2,ops3,bz1,bz2,bz3,bz12},
      MyAssert[CHANNELS == 3];

      (* Operator indexes which correspond to the different channels *)
      If[NROPS == 6,
        ch1ops = {3,4};
        ch2ops = {1,5};
        ch3ops = {2,6};
      ];
      If[NROPS == 3,
        ch1ops = {1};
        ch2ops = {2};
        ch3ops = {3};
      ];

      (* Operators in each channel *)
      ops1 = basisops[[ch1ops]];
      ops2 = basisops[[ch2ops]];
      ops3 = basisops[[ch3ops]];
      MyPrint["ops1=", ops1];
      MyPrint["ops2=", ops2];
      MyPrint["ops3=", ops3];

      (* Basis sets for each sub-system *)
      bz1 = spinlessbasis @ qbasis[ops1];
      bz2 = spinlessbasis @ qbasis[ops2];
      bz3 = spinlessbasis @ qbasis[ops3];
      MyPrint["bz1=", bz1];
      MyPrint["bz2=", bz2];
      MyPrint["bz3=", bz3];

      makebasis[basisops]; (* Reset the definition!! *)

      bz12 = basistensorproduct[bz1, bz2, {#1[[1]], #2[[1]]} &];
      bz = basistensorproduct[bz12, bz3, {#1[[1]], #1[[2]], #2[[1]]} &];

      bz = Map[{First[#] + NROPS/3, Last[#]}&, bz]; (* Shift *)

      (* Projection to single spin is performed later; see below. *)
    ];
  ];

  basisexceptions = {"QST", "QSTZ", "QSZTZ", "SPSU2T", "QSC3", "SPSU2C3", "QJ"};
  If[MemberQ[basisexceptions, SYMTYPE],
    loadmodule[SYMTYPE <> "basis.m"];
  ];

  MyVPrint[2, "Baza:", bz];

  (*** STEP 2: Transform to correct naming convention, apply necessary
  transformations, etc. ***)

  If[isQS[],
    (* Transform to (Q,2S+1) naming convention as used in the NRG program. *)
    bz = Map[ {{#[[1,1]], 2*#[[1,2]]+1}, #[[2]]}&, bz];
    bvc = bzop2bzvc[bz, vak];
    degnr[{q_Integer, ss_Integer, i___}] := ss;
  ];

  If[isISO[] || isISOSZ[] || isSU2[],  (* Common to ISO, ISOSZ, and SU2 *)
    (* Generate isospin operators *)
    Tops = Map[isospin[#, nnop[#]]&, basisops];
    Tsum = Plus @@ Tops;
    Tminus = Simplify[Tsum[[1]] - I Tsum[[2]]];

    (* IMPORTANT: check if nnop[] is defined correctly. Wrong signs
    lead to results that are wrong in very subtle manner, so one
    should be very careful here! *)

    MyPrint["nnops=", Map[{#, nnop[#]}&, basisops]];
    MyPrint["Tminus=", Tminus];

    (* Perform the transformation *)
    bvc = bzop2bzvc[bz, vak]; (* ! *)
    bvc = transformQStoIS[bvc]; (* The same for QS, QSZ & Q input! *)
    bz = bzvc2bzop[bvc];

    If[isISO[], (* ISO specific *)
      (* Transform to (2I+1,2S+1) naming convention. *)
      bz = Map[ {{2*#[[1,1]]+1, 2*#[[1,2]]+1}, #[[2]]}&, bz];
      degnr[{ii_Integer, ss_Integer, j___}] := ii ss;
    ];

    If[isISOSZ[], (* ISOSZ specific *)
      (* Transform to (2I+1,2Sz+1) naming convention. *)
      bz = Map[ {{2*#[[1,1]]+1, 2*#[[1,2]]+1}, #[[2]]}&, bz];
      degnr[{ii_Integer, ssz_Integer, j___}] := ii;
    ];

    If[isSU2[], (* SU2 specific *)
      (* Transform to (2I+1) naming convention. *)
      bz = Map[ {{2*#[[1,1]]+1}, #[[2]]}&, bz];
      degnr[{ii_Integer, j___}] := ii;
    ];

    bvc = bzop2bzvc[bz, vak];
    MyVPrint[2, "ISO/ISOSZ/SU2 baza=", bz];
  ];

  If[isDBLQSZ[],
    bz = Map[ {{#[[1,1]], #[[1,2]], 2*#[[1,3]]+1}, #[[2]]}&, bz ];
    bvc = bzop2bzvc[bz, vak];
    MyVPrint[2, "DBLQSZ bz=", bz];
    MyVPrint[2, "DBLQSZ bvc=", bvc];
    degnr[{_, _, _}] = 1;
  ];

  If[isDBLSU2[],
    bz = Map[ {{2*#[[1,1]]+1, 2*#[[1,2]]+1}, #[[2]]}&, bz ];
    bvc = bzop2bzvc[bz, vak];
    MyVPrint[2, "DBLSU2 bz=", bz];
    MyVPrint[2, "DBLSU2 bvc=", bvc];
    degnr[{ii1_, ii2_}] = ii1 ii2;
  ];

  If[isDBLISOSZ[],
    bz = Map[ {{2*#[[1,1]]+1, 2*#[[1,2]]+1, 2*#[[1,3]]+1}, #[[2]]}&, bz ];
    bvc = bzop2bzvc[bz, vak];
    MyVPrint[2, "DBLSU2 bz=", bz];
    MyVPrint[2, "DBLSU2 bvc=", bvc];
    degnr[{ii1_, ii2_, _}] = ii1 ii2;
  ];

  If[isQSZ[],
    (* Transform to (Q,2Sz+1) naming convention as used in the QSZ program. *)
    bz = Map[ {{#[[1,1]], 2*#[[1,2]]+1}, #[[2]]}&, bz];
    bvc = bzop2bzvc[bz, vak];
    degnr[{q_Integer, sz_Integer, i___}] := 1;
  ];

  If[isP[] || isPP[] || isNONE[] || isU1[] || isSL[] || isSL3[],
    bvc = bzop2bzvc[bz, vak];
    degnr[{___}] := 1;
  ];

  If[isSPSU2[] || isSPSU2LR[],
    bz = Map[ {{2*#[[1,1]]+1}, #[[2]]}&, bz];
    bvc = bzop2bzvc[bz, vak];
    degnr[{ss_Integer, i___}] := ss;
  ];

  If[isSPU1[] || isSPU1LR[],
    bz = Map[ {{2*#[[1,1]]+1}, #[[2]]}&, bz];
    bvc = bzop2bzvc[bz, vak];
    degnr[{___}] := 1;
  ];

  (* Nothing to do for QST, QSTZ and SPSU2T *)

  nrstates1 = Plus @@ Map[(degnr @ #[[1]] * Length @ #[[2]])&, bvc];
  nrstates2 = Plus @@ Map[(degnr @ #[[1]] * Length @ #[[2]])&, bz];
  If[nrstates1 != nrstates2,
    MyError["BUG: bz !~ bvc"];
  ];

  MyVPrint[2, "Baza (step 2):", bz];

  (*** Step 3: Now it's time to apply eventual transformations of the
  basis states. ***)

  If[ BASISRULE =!= "",
    If[StringQ[BASISRULE],
      rule = ToExpression[BASISRULE];
      MyPrint["Transformation rule: ", BASISRULE, " -> ", rule],
    (* else *)
      rule = BASISRULE;
      MyPrint["Transformation rule: ", rule];
    ];

    (* Perform the transformation *)
    TT3 = Timing[ bvc = transformbasis[bvc, rule]; ];
    bz = bzvc2bzop[bvc];

    MyVPrint[2, "PROJECTED baza=", bz];

    (* Count all states *)
    nrstates = Plus @@ Map[(degnr[ #[[1]] ] Length[ #[[2]] ])&, bvc];
    MyPrint["PROJECTED NR=", nrstates, " TT3=", Chop[TT3, 4] ];
  ];


  (*** Step 4: Add phonons, if requested. ***)
  (* Note: phonons need to be added prior to generating parity-adapted
  basis, since in some models phonon kets transform under the mirror
  symmetry (for example in ONE/COM model). *)
  If[ MAKEPHONON =!= Null,
    MyVPrint[1, "Adding phonons"];

    bz = transformtoPH[bz, nph];
    MyVPrint[2, "PHONON baza (op)=", bz];

    bvc = bzop2bzvc[bz, vak];
    MyVPrint[2, "PHONON baza (vc)=", bvc];
  ];


  (*** Step 5: Generate parity-adapted basis ***)
  dolr[] := Module[{},
    If[ (isLR[] || option["LRTRICK"]) && CHANNELS == 2,
      If[lrchain === {},
        MyError["Error: LR symmetry but lrchain not defined in the model."];
      ];
      MyPrint["Generating parity-adapted basis. lrchain=", lrchain];

      bvc = transformtoLRvc[bz, lrchain, vak, lrextrarule];
      bz = bzvc2bzop[bvc];
      MyVPrint[2, "LR (parity) baza=", bz];

      If[option["LRTRICK"],
        droplrindex[{qn_, states_}] := {Drop[qn, -1], states};
        bvc = mergebasis[Map[droplrindex, bvc]];
        bz = mergebasis[Map[droplrindex, bz]];
        MyVPrint[2, "LRTRICK baza=", bvc];
        MyVPrint[2, "LRTRICK baza=", bz];
      ];
    ];
    If[ isLR[] && CHANNELS == 1, (* for testing purposes mostly! *)
      bvc = Map[{Append[#[[1]], 1], #[[2]]}&, bvc];
      bz = Map[{Append[#[[1]], 1], #[[2]]}&, bz];
    ];
  ];

  If[ !option["LRSPIN"], dolr[] ]; (* LR transform before adding spin kets *)

  (*** Step 6: Add spin kets if requested. ***)

  If[ MAKESPINKET =!= Null,
   MyPrint["MAKESPINKET=", MAKESPINKET];

   If[isQS[] || isSPSU2[] || isSPSU2LR[] || isQST[] || isQSTZ[] ||
     isSPSU2T[] || isQSC3[] || isSPSU2C3[],
      If[ BZSPIN === Null,
        bzspin = spinbasis[MAKESPINKET];
        MyPrint["spin basis=", bzspin];
        (* Maximum spin state *)
        bzspin = {First @ spinbasis[MAKESPINKET]};
        MyPrint["spin basis=", bzspin],
      (* else *)
        bzspin = BZSPIN;
        MyPrint["using manually set spin basis=", bzspin]
      ];

      Off[ClebschGordan::phy];

      (* Recall: isQS covers QS and QSLR *)
      If[SYMTYPE == "QS",
        spin1fnc = ((#[[2]]-1)/2 &);
        spincombinefncNRG[s_, {q1_, _}, {___}] := {q1, 2s+1};
      ];
      If[SYMTYPE == "QSLR",
        spin1fnc = ((#[[2]]-1)/2 &);
        spincombinefncNRG[s_, {q1_, _, p1_}, {___}] := {q1, 2s+1, p1};
      ];
      If[isSPSU2[],
        spin1fnc = ((#[[1]]-1)/2 &);
        spincombinefncNRG[s_, {___}, {___}] := {2s+1};
      ];
      If[isSPSU2T[] || isSPSU2LR[] || isSPSU2C3[],
        spin1fnc = ((#[[1]]-1)/2 &);
        spincombinefncNRG[s_, {s1_, p1_}, {___}] := {2s+1, p1};
      ];
      If[isQST[] || isQSTZ[] || isQSC3[],
        spin1fnc = ((#[[2]]-1)/2 &);
        spincombinefncNRG[s_, {q1_, s1_, t1_}, {___}] := {q1, 2s+1, t1};
      ];

      bz = SU2basistensorproduct[bz, bzspin,
        spin1fnc, First,
        spincombinefncNRG,
        fixspin, fixspinket];
      bvc = bzop2bzvc[bz, vak];
    ];

    (* Modified ???addspinket[] functions to conform to the (2Sz+1) convention. *)
    QSZaddspinketNRG[bvc_List, SP_] := Module[{spinbz},
      spinbz = spinbasis[SP];
      MyVPrint[2, "spinbz=", spinbz];
      basistensorproduct[spinbz, bvc, Function[{qn1, qn2},
        Module[{x=qn2}, x[[2]] += 2*First[qn1]; x] ]  (* Note the factor 2! *)
      ]
    ];

    If[isQSZ[] || isISOSZ[] || isQSZTZ[],
      bvc = QSZaddspinketNRG[bvc, MAKESPINKET];
      bz = bzvc2bzop[bvc];
    ];

    SPU1addspinketNRG[bvc_List, SP_] := Module[{spinbz},
      spinbz = spinbasis[SP];
      basistensorproduct[spinbz, bvc, Function[{qn1, qn2},
        Module[{x=qn2}, x[[1]] += 2*First[qn1]; x] ]  (* Note the factor 2! *)
      ]
    ];

    If[isSPU1[] || isSPU1LR[],
      bvc = SPU1addspinketNRG[bvc, MAKESPINKET];
      bz = bzvc2bzop[bvc];
    ];

    SZaddspinketNRG[bvc_List, SP_] := Module[{spinbz},
      spinbz = spinbasis[SP];
      basistensorproduct[spinbz, bvc, Function[{qn1, qn2},
        Module[{x=qn2}, x[[3]] += 2*First[qn1]; x] ]  (* Note the factor 2! *)
      ]
    ];

    If[isDBLQSZ[] || isDBLISOSZ[],
      bvc = SZaddspinketNRG[bvc, MAKESPINKET];
      bz = bzvc2bzop[bvc];
    ];

    If[isDBLSU2[] || isSU2[] || isU1[] || isP[] || isPP[] || isNONE[],
      QSZaddspinketU1[bvc_List, SP_] := Module[{spinbz},
        spinbz = spinbasis[SP];
        basistensorproduct[spinbz, bvc, Function[{qn1, qn2}, qn2]]
      ];
      bvc = QSZaddspinketU1[bvc, MAKESPINKET];
      bz = bzvc2bzop[bvc];
    ];

    MyVPrint[2, "SPINKET baza=", bz];
    MyVPrint[2, "SPINKET baza=", bvc];
  ];

  (*** Step 6b: Add orbital kets, if requested ***)
  If[ MAKEORBKET =!= Null,
    MyPrint["MAKEORBKET=", MAKEORBKET];

    If[isQST[],
     Module[{orbbz, tzfnc, oT, oS, opTM, TM, tzdown},
      orbbz = { First @ spinbasis[MAKEORBKET] };
      MyPrint["orb basis=", orbbz];

      (* Spin kets must already be defined. *)
      orbbz = orbbz /. ket[i_] :> ket[Null, i];
      bz = bz /. ket[i_] :> ket[i, Null];
      MyPrint["orbbz=", orbbz];
      MyPrint["bz=", bz];
      tz1fnc = (#[[3]] &);
      tzcombinefncNRG[t_, {q1_, s1_, t1_}, {___}] := {q1, s1, t};

      fixtzket[k_, t_, t_] := k;
      fixtzket[k_, t_, tz_] /; tz < t := k /. ket[s_, i_] :> ket[s, i - (t - tz)];

      fixtz[op_, t_, t_] := op;
      fixtz[op_, t_, tz_] /; tz < t := Nest[tzdown, op, t - tz];

      oT = 1;
      oS = 1/2;
      tminus[op_[q___, sz_]] :=
        VMV[Table[op[CR, q, tz, sz], {tz, oT, -oT, -1}], spinmatrixM[oT],
          Table[op[AN, q, tz, sz], {tz, oT, -oT, -1}]];
      msz[1/2] = UP;
      msz[-1/2] = DO;
      opTM[op_[q___]] := Sum[tminus[op[q, msz[sz]]], {sz, -oS, oS}];
      TM = Total@Map[opTM, {f[]}];
      MyPrint["TM=", TM];

      normalizeop[op_] := Simplify[op/Simplify[normop[op]]];
      tzdown[op_] := normalizeop @ zeroonvac[Expand[nc[TM, op]]]; (* normalize!! *)

      zeroonvac[x:nc[a___, y:op_?fermionQ[i_, j___], ket[___]]] :=
        If[isannihilation[y], 0, x];

      bz = SU2basistensorproduct[bz, orbbz,
        tz1fnc, First,
        tzcombinefncNRG,
        fixtz, fixtzket];
      bz = Expand[bz];
      bvc = bzop2bzvc[bz, vak];
     ];
    ];

    If[isQSTZ[] || isQSZTZ[],
     Module[{orbbz},
      orbbz = spinbasis[MAKEORBKET];
      MyPrint["orb basis=", orbbz];

      (* Spin kets must already be defined. *)
      orbbz = orbbz /. ket[i_] :> ket[Null, i];
      bvc = bvc /. ket[i_] :> ket[i, Null];

      bvc = basistensorproduct[orbbz, bvc, Function[{qn1, qn2}, {qn2[[1]], qn2[[2]], qn2[[3]]+qn1[[1]]}] ];
      bz = bzvc2bzop[bvc];
     ];
    ];

    If[isSPSU2T[],
      MyError["Not implemented."];
    ];

    If[!isORB[] && !isQJ[],
      orbbz = spinbasis[MAKEORBKET];
      MyPrint["orb basis=", orbbz];

      (* Spin kets must already be defined. *)
      orbbz = orbbz /. ket[i_] :> ket[Null, i];
      bvc = bvc /. ket[i_] :> ket[i, Null];

      bvc = basistensorproduct[orbbz, bvc, Function[{qn1, qn2}, qn2]];
      bz = bzvc2bzop[bvc];
    ];

    MyVPrint[2, "ORBKET baza bz=", bz];
    MyVPrint[2, "ORBKET baza bvc=", bvc];
  ];

  (*** Step 6c ***)
  If[ option["LRSPIN"], dolr[] ]; (* LR transform *after* adding spin kets *)

  (*** Step 7: spinless-fermions hack ***)

  (* Support for spinless fermions is implemented as a simple hack:
  we keep only spin-up fermions in the basis! *)

  If[parambool["spinless"] || isSL[],
    MyPrint["spinless=true"];
    bz = spinlessbasis[bz];
    If[ SYMTYPE == "QS" || SYMTYPE == "QSZ",
      bz = Map[{First[#] + {NROPS, 0}, Last[#]}&, bz]; (* Shift *)
    ];
    If[ SYMTYPE == "U1" || SYMTYPE == "SL",
      bz = Map[{First[#] + {NROPS}, Last[#]}&, bz]; (* Shift *)
    ];
    bvc = bzop2bzvc[bz, vak];
  ];

  (* Hook for possible changes of the basis states *)
  hookfile["hook_basis"];

  MyPut[bvc, basisfilename];
  timeadd["basis"];
];

(* At this point bvc contains the basis in the occupation number
representation. bz is not used directly past this point. *)

bz =.; (* Safety measure! *)
