(*
  NRG Ljubljana - Additional operator definitions 
  (c) Rok Zitko, rok.zitko@ijs.si, 2005-2016
*)

Print["operators.m started"];

setkm[t_, sp_, op_, opimp_] := 
 Module[{opud, w1, w2, w3, s1, s2, s3, komp, temp},
  opud = {op[AN, UP], op[AN, DO]};
  w1 = 1/2 PauliX . opud;
  w2 = 1/2 PauliY . opud;
  w3 = 1/2 PauliZ . opud;
  {s1, s2, s3} = spinxyz[opimp[]];
  komp = If[sp == UP, 1, 2];
  temp = w1[[komp]]~nc~s1 + w2[[komp]]~nc~s2 + w3[[komp]]~nc~s3 //Expand;
  If[t == CR, conj[temp], temp (-1)^sp]
];

Module[{t}, 
(* We start a new table to avoid interference. *)
t = {}; 
  (* === Misc. operators for testing === *)
  
  (* Identity operator. <I>=1. *)
  t = Join[t, mtSingletOp["I", 1] ];

  (* Parts of the total Hamiltonian *)
  
  If[calcopq["Himp"],
    Module[{op},
      If[!ValueQ[Himp],
        Himp = H1; (* Default impurity Hamiltonian *)
      ];
      op = Expand[ Himp /. params ];
      MPVCFAST = False; (* Workaround *)
      t = Join[t, mtSingletOp["Himp", op]];
      MPVCFAST = True;
    ];
  ];

  If[calcopq["Hhyb"],
    Module[{op},
      If[!ValueQ[Hhyb],
        Hhyb = Hc; (* Default hybridization Hamiltonian *)
      ];
      op = Expand[ Hhyb /. params ];
      MPVCFAST = False;
      t = Join[t, mtSingletOp["Hhyb", op]];
      MPVCFAST = True;
    ];
  ];

  (* Potential energy *)
  If[calcopq["Hpot"],
    Module[{op},
      If[!ValueQ[Hpot],
        Hpot = 0; (* Default is zero *)
      ];
      op = Expand[ Hpot /. params ];
      MPVCFAST = False; (* Workaround *)
      t = Join[t, mtSingletOp["Hpot", op]];
      MPVCFAST = True;
    ];
  ];
  (* === Electron gas  === *)

  (* 'Occupancy' of electron gas on the dot *)
  t = Join[t, mtSingletOp["n_f", number[f[0]] ]];
  t = Join[t, mtSingletOp["q_f", number[f[0]]-1 ]];

  (* Square of occupancy of electron gas *)
  t = Join[t, mtSingletOp["n_f^2", pow[number[f[0]], 2] ] ];
  t = Join[t, mtSingletOp["q_f^2", pow[number[f[0]]-1, 2] ] ];

  (* Pairing at the position of the impurity *)
  t = Join[t, mtSingletOp["pair_f", nc[f[CR, 0, UP], f[CR, 0, DO]] ] ];

  (* Spectral density of electron gas *)
  t = Join[t, mtDoubletOp["A_f", f[#1, 0, #2]& ]];

  t = Join[t, mtDoubletOp["A_F", f ]];
      
  (* Conjugated A_f operator for computing anomalous Green's
     function. See ireducTable[] and ireducMatrixSpeedy[] in initial.m *)
  If[calcopq["Ac_f"],
    AppendTo[t, {"dAc_f"}];
    t = Join[t, ireducTable[ ((-1)^#2 f[1-#1, 0, 1-#2])& ] ];
  ];
    
  If[calcopq["A_f_u"], (* Spectral density of electron gas (spin UP) *)
    AppendTo[t, {"dA_f_u"}];
    t = Join[t, ireducTable[ f[0], UP ] ];
  ];
    
  If[calcopq["A_f_d"], (* Spectral density of electron gas (spin DOWN) *)
    AppendTo[t, {"dA_f_d"}];
    t = Join[t, ireducTable[ f[0], DO ] ];
  ];

  t = Join[t, mtSingletOp["SXf0", spinx[f[0]] ] ];
  t = Join[t, mtSingletOp["SX2f0", pow[spinx[f[0]],2] ] ];
  t = Join[t, mtSingletOp["SYf0", spiny[f[0]] ] ];
  t = Join[t, mtSingletOp["SY2f0", pow[spiny[f[0]],2] ] ];
  t = Join[t, mtSingletOp["SZf0", spinz[f[0]] ] ];
  t = Join[t, mtSingletOp["SZ2f0", pow[spinz[f[0]],2] ] ];

  t = Join[t, mtTripletOp["sigma_f", f[] ]]; (* Spin operator *)

  t = Join[t, mtTripletOp["sigma_F", F ]]; (* Spin operator *)
  t = Join[t, mtOrbTripletOp["orbsigma_F", F]]; (* Orbital moment operator *)
  (* For these two, one needs to define spinz[F], spinplus[F], spinminus[F],
  and orbmomentz[F], orbmomentplus[F], orbmomentminus[F]. *)

  (* Additional sites along the Wilson chain *)

  If[calcopq["A_f_i1"],
    AppendTo[t, {"dA_f_i1"}];
    t = Join[t, ireducTable[ f[#1, 0, 1, #2]& ] ];
  ];
  If[calcopq["A_f_i2"],
    AppendTo[t, {"dA_f_i2"}];
    t = Join[t, ireducTable[ f[#1, 0, 2, #2]& ] ];
  ];
  If[calcopq["A_f_i3"],
    AppendTo[t, {"dA_f_i3"}];
    t = Join[t, ireducTable[ f[#1, 0, 3, #2]& ] ];
  ];
  If[calcopq["A_f_i4"],
    AppendTo[t, {"dA_f_i4"}];
    t = Join[t, ireducTable[ f[#1, 0, 4, #2]& ] ];
  ];

  (* === Add site 'd' === *)

  (* Occupancy of embedded dot *)
  t = Join[t, mtSingletOp["n_d", number[d[]] ]];
  t = Join[t, mtSingletOp["q_d", number[d[]]-1 ]];

  (* Occupancy of embedded dot, spin up *)
  t = Join[t, mtSingletOp["n_d_up", number[d[], UP] ] ];
  t = Join[t, mtSingletOp["n_d_u", number[d[], UP] ] ];

  (* Occupancy of embedded dot, spin down *)
  t = Join[t, mtSingletOp["n_d_down", number[d[], DO] ] ];
  t = Join[t, mtSingletOp["n_d_d", number[d[], DO] ] ];

  (* Occupancy of embedded dot, double occupancy *)
  t = Join[t, mtSingletOp["n_d_ud", number[d[], DO] ~ nc ~ number[d[], UP] ] ];

  (* "Local-moment fraction" *)
  (* flm=n-2n_u n_d=2n-n^2=Projector[1] *)
  t = Join[t, mtSingletOp["flm", number[d[]] - 2 number[d[], DO] ~ nc ~ number[d[], UP] ] ];

  (* Square of occupancy of embedded dot *)
  t = Join[t, mtSingletOp["n_d^2", pow[number[d[]],2] ] ];

  (* Square of occupancy of embedded dot, spin up *)
  t = Join[t, mtSingletOp["n_d_up^2", pow[number[d[], UP],2] ] ];

  (* Square of occupancy of embedded dot, spin down *)
  t = Join[t, mtSingletOp["n_d_down^2", pow[number[d[], DO],2] ] ];

  (* Charge fluctuations of embedded dot, for isospin symmetry! *)
  t = Join[t, mtSingletOp["q_d^2", pow[number[d[]]-1, 2] ] ];

  (* Hopping amplitude from impurity to the lead *)
  t = Join[t, mtSingletOp["hop0", hop[d[], f[0]] ] ];

  (* Hopping amplitude from impurity to the lead, squared *)
  t = Join[t, mtSingletOp["hop0^2", pow[hop[d[], f[0]],2] ] ];

  (* Hopping-n_d combination; used for testing the sum rules. *)
  t = Join[t, mtSingletOp["hop0-n_d", nc[hop[d[], f[0]], number[d[]]] ] ];

  (* Current operator *) 
  t = Join[t, mtSingletOp["cur0", current[d[], f[0]] ] ];
 
  (* Isospin operator, x component *)
  t = Join[t, mtSingletOp["pair_d", nc[d[CR,UP], d[CR, DO]] ] ];
  t = Join[t, mtSingletOp["Ix_d", isospinx[d[]] ] ];
  t = Join[t, mtSingletOp["Ix_d^2", pow[isospinx[d[]],2] ] ];

  (* Isospin operator, y component *)
  t = Join[t, mtSingletOp["Iy_d", isospiny[d[]] ] ];
  t = Join[t, mtSingletOp["Iy_d^2", pow[isospiny[d[]],2] ] ];

  (* Isospin operator, z component (Iz_d is equal to 1/2 q_d) *)
  t = Join[t, mtSingletOp["Iz_d", isospinz[d[]] ] ];
  t = Join[t, mtSingletOp["Iz_d^2", pow[isospinz[d[]],2] ] ];
  
  (* Majorana operators *)
  (* If fermion parity is conserved, they all have zero expectation value. *)
  t = Join[t, mtSingletOp["maj1up", (d[CR,UP]+d[AN,UP])/2 ] ];
  t = Join[t, mtSingletOp["maj1do", (d[CR,DO]+d[AN,DO])/2 ] ];
  t = Join[t, mtSingletOp["maj2up", I (d[CR,UP]-d[AN,UP])/2 ] ];
  t = Join[t, mtSingletOp["maj2do", I (d[CR,DO]-d[AN,DO])/2 ] ];
  
  (* Transverse isospin squared *)
  t = Join[t, mtSingletOp["Ixy_d^2", pow[isospinx[d[]],2] + pow[isospiny[d[]],2] ] ];
  
  t = Join[t, mtSingletOp["pnup_d", 
    nc[number[d[], UP], d[CR,UP], d[CR, DO]] ] ];
  t = Join[t, mtSingletOp["pndo_d", 
    nc[number[d[], DO], d[CR,UP], d[CR, DO]] ] ];

  t = Join[t, mtSingletOp["pn_d",
    nc[number[d[]], d[CR,UP], d[CR, DO]] ] ];
  
    (* Recall that <AB> = \int_-\infty^0 (-1/Pi)Im <<A;B>>(omega) at
       zero temperature (GS). *)
    
  (* Spin-resolved current/hopping operators *)
  (* R/I = Re vs. Im, U/D = UP vs. DO *)
  t = Join[t, mtSingletOp["jRU",     hop[d[], f[0], UP] ] ];
  t = Join[t, mtSingletOp["jRD",     hop[d[], f[0], DO] ] ];
  t = Join[t, mtSingletOp["jIU", current[d[], f[0], UP] ] ];
  t = Join[t, mtSingletOp["jID", current[d[], f[0], DO] ] ];

  t = Join[t, mtSingletOp["jIU^2", pow[current[d[], f[0], UP],2] ] ];
  
  (* Note that jRU+jRD <=> hop0 *)

  (* Equivalent operators for the 2nd channel *)
  t = Join[t, mtSingletOp["jRUr",     hop[d[], f[1], UP] ] ];
  t = Join[t, mtSingletOp["jRDr",     hop[d[], f[1], DO] ] ];
  t = Join[t, mtSingletOp["jIUr", current[d[], f[1], UP] ] ];
  t = Join[t, mtSingletOp["jIDr", current[d[], f[1], DO] ] ];

  (* Current from one lead to another, spin UP only. Use this 
     for spin-less models! *)
  t = Join[t, mtSingletOp["cur", 
    current[d[], f[0], UP] - current[d[], f[1], UP]
  ] ];

  (* Both spin orientations. Use this for the spin-full models! *)
  t = Join[t, mtSingletOp["curUD",
    current[d[], f[0]] - current[d[], f[1]]
  ] ];
  

 If[CHANNELS >= 2,
    (* 'Occupancy' of electron gas on the dot *)
    t = Join[t, mtSingletOp["q_f1", number[f[1]]-1 ] ];

    (* Square of occupancy of electron gas *)
    t = Join[t, mtSingletOp["q_f1^2", pow[number[f[1]]-1, 2] ] ];

    t = Join[t, mtSingletOp["hop1", hop[d[], f[1]] ] ];

    (* Symmetric and antisymmetric combinations *)
    t = Join[t, mtSingletOp["hops",
      1/Sqrt[2] (hop[d[], f[0]] + hop[d[], f[1]]) ] ];

    If[isLR[],
      t = Join[t, mtGeneralOp["hopa",
        1/Sqrt[2] (hop[d[], f[0]] - hop[d[], f[1]]) ] ],
    (* else *)
      t = Join[t, mtSingletOp["hopa",
        1/Sqrt[2] (hop[d[], f[0]] - hop[d[], f[1]]) ] ]
    ];

    (* Squares of the above. Both are even parity operators! *)
    t = Join[t, mtSingletOp["hops^2",
      pow[1/Sqrt[2] (hop[d[], f[0]] + hop[d[], f[1]]), 2] ] ];
    t = Join[t, mtSingletOp["hopa^2",
      pow[1/Sqrt[2] (hop[d[], f[0]] - hop[d[], f[1]]), 2] ] ];
  ];
  
  If[CHANNELS >= 3,
    (* 'Occupancy' of electron gas on the dot *)                                             
    t = Join[t, mtSingletOp["q_f2", number[f[2]]-1 ] ];
  
    (* Square of occupancy of electron gas *)                                                
    t = Join[t, mtSingletOp["q_f2^2", pow[number[f[2]]-1, 2] ] ];                            
  
    t = Join[t, mtSingletOp["hop2", hop[d[], f[2]] ] ]; 

    t = Join[t, mtSingletOp["q_f012",
      number[f[0]] + number[f[1]] + number[f[2]] - 3 ] ];
  ];
  
  (* Spectral density of embedded dot *)
  t = Join[t, mtDoubletOp["A_d", d ]];
    
  (* NOTE: A_d_u, A_d_d etc. are to be used with U1 symmetry type. There
     are also U_d, D_d etc. operators which correspont to spin-projected
     operators. rz, 14 Aug 2014 *)
  If[calcopq["A_d_u"], (* Spectral density (spin UP) *)
    AppendTo[t, {"dA_d_u"}];
    t = Join[t, ireducTable[ d[#1, #2]&, UP ] ];
  ];
  
  If[calcopq["A_d_d"], (* Spectral density (spin DOWN) *)
    AppendTo[t, {"dA_d_d"}];
    t = Join[t, ireducTable[ d[#1, #2]&, DO ] ];
  ];
  
  (* Spin-projected doublet operators. Note that these are not SU(2)_spin doublets,
     thus they should never be used in calculations with such symmetry types. They
     are fine for QSZ and SPU1 symtypes. *)
  If[calcopq["U_d"],
    AppendTo[t, {"dU_d"}];
    t = Join[t, ireducTable[ (If[#2 == UP, 1, 0] d[#1, #2])& ] ];
  ];

  If[calcopq["D_d"],
    AppendTo[t, {"dD_d"}];
    t = Join[t, ireducTable[ (If[#2 == DO, 1, 0] d[#1, #2])& ] ];
  ];

  (* Version of the above, explicitly defined to be the creation operators.
     Should be exactly equivalent to U_d & D_d for QSZ and SPU1 symmetries, 
     since #1 is always CR for those cases. *)
  If[calcopq["CU_d"],
    AppendTo[t, {"dCU_d"}];
    t = Join[t, ireducTable[ (If[#1 == CR, 1, 0] If[#2 == UP, 1, 0] d[#1, #2])& ] ];
  ];

  If[calcopq["CD_d"],
    AppendTo[t, {"dCD_d"}];
    t = Join[t, ireducTable[ (If[#1 == CR, 1, 0] If[#2 == DO, 1, 0] d[#1, #2])& ] ];
  ];

  (* Important implementation note: argument #2 is the spin of the excitation described by 
     the operator, i.e., which invariant subspaces it links. For SPU1 #1 is always CR.
     Creation operator will increase spin, annihilation operator will decrease it. 
     We therefore need to invert (1-#2) the argument to obtain the correct operator for 
     particle annihilation. *)
  If[calcopq["AD_d"],
    AppendTo[t, {"dAD_d"}];
    t = Join[t, ireducTable[ (If[#1 == CR, 1, 0] If[#2 == UP, 1, 0] d[1-#1, 1-#2])& ] ];
  ];

  If[calcopq["AU_d"],
    AppendTo[t, {"dAU_d"}];
    t = Join[t, ireducTable[ (If[#1 == CR, 1, 0] If[#2 == DO, 1, 0] d[1-#1, 1-#2])& ] ];
  ];

  (* For SYMTYPE == "NONE" *)
  If[calcopq["A_d_0"], (* Spectral density (type 0) *)
    AppendTo[t, {"dA_d_0"}];
    t = Join[t, ireducTable[ d[#1, #2]&, 0 ] ];
  ];
  If[calcopq["A_d_1"], (* Spectral density (type 1) *)
    AppendTo[t, {"dA_d_1"}];
    t = Join[t, ireducTable[ d[#1, #2]&, 1 ] ];
  ];
  If[calcopq["A_d_2"], (* Spectral density (type 2) *)
    AppendTo[t, {"dA_d_2"}];
    t = Join[t, ireducTable[ d[#1, #2]&, 2 ] ];
  ];
  If[calcopq["A_d_3"], (* Spectral density (type 3) *)
    AppendTo[t, {"dA_d_3"}];
    t = Join[t, ireducTable[ d[#1, #2]&, 3 ] ];
  ];

  (* For SYMTYPE == "NONE" *)
  (* Note that this is not consistent with A_d & Ac_d operators for other
     symmetry types. In particular, A_d is the creation operator, while
     sA_d are annihilation operators. TO DO: fix this?! *)
  t = Join[t, mtSingletOp["sA_d_u",  d[AN, UP] ] ];
  t = Join[t, mtSingletOp["sAc_d_u", d[CR, UP] ] ];
  t = Join[t, mtSingletOp["sA_d_d",  d[AN, DO] ] ];
  t = Join[t, mtSingletOp["sAc_d_d", d[CR, DO] ] ];
  
  (* See symmetry_test__self_operators_for_ISO.nb. B_d is the operator which 
  results from [Iz^2,d], i.e. from the Coulomb repulsion term expressed in an
  isospin-symmetric way. *)
  If[calcopq["B_d"],
    AppendTo[t, {"dB_d"}];
    t = Join[t, ireducTable[
      d[#1, #2] ~ nc ~ (d[CR,1-#2] ~ nc ~ d[AN,1-#2] - 1/2)  If[#1==AN,-1,1]&,
    nnop[d[]] ] ];
  ];

  (* Conjugated A_d operator for computing anomalous Green's
     function. See ireducTable[] and ireducMatrixSpeedy[] in initial.m *)
  (* Recall: UP=1, DO=0 *)
  (* For spin-UP --> -d[1-#2, DO] *)
  (* For spin-DO --> +d[1-#2, UP] *)
  (* WARNING: the minus sign is conventional, some people use different definitions! 
     Be careful! *)

  If[calcopq["Ac_d"],
    AppendTo[t, {"dAc_d"}];
    t = Join[t, ireducTable[ ((-1)^#2 d[1-#1, 1-#2])& ] ];
  ];
  
  (* QSZ only! *)
  If[calcopq["nAc_d"],
    AppendTo[t, {"dnAc_d"}];
    t = Join[t, ireducTable[ (nc[number[d[], #2], d[1-#1, 1-#2]])& ] ];
  ];

  (* As Ac_d, but with no sign flipping!! Don't use this for SU(2)
     symmetric cases, since it doesn't have the correct symmetry properties!! *)
  If[calcopq["AC_d"],
    AppendTo[t, {"dAC_d"}];
    t = Join[t, ireducTable[ (d[1-#1, 1-#2])& ] ];
  ];

  t = Join[t, mtTripletOp["sigma_d", d[] ]]; (* Spin operator *)

  (* Support for automatically generated operator for calculating
     the generalized Green's function for the self-energy trick.
     Added 6. 6. 2012 *)
  If[calcopq["self_d"],
    AppendTo[t, {"dself_d"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopd ]];
    MPVCFAST = True;
  ];

  If[calcopq["selfc_d"],
    AppendTo[t, {"dselfc_d"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopcd ]];
    MPVCFAST = True;
  ];

  (* To be used with SYMTYPE=U1, etc. *)
  (* Added 6. 9. 2012 *)
  If[calcopq["self_d_u"],
    AppendTo[t, {"dself_d_u"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopd, UP ]];
    MPVCFAST = True;
  ];

  If[calcopq["self_d_d"],
    AppendTo[t, {"dself_d_d"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopd, DO ]];
    MPVCFAST = True;
  ];
    
  (* By analogy with self_d *)  

  If[calcopq["self_a"],
    AppendTo[t, {"dself_a"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopa ]];
    MPVCFAST = True;
  ];

  If[calcopq["selfc_a"],
    AppendTo[t, {"dselfc_a"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopca ]];
    MPVCFAST = True;
  ];

  If[calcopq["self_a_u"],
    AppendTo[t, {"dself_a_u"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopa, UP ]];
    MPVCFAST = True;
  ];

  If[calcopq["self_a_d"],
    AppendTo[t, {"dself_a_d"}];
    MPVCFAST = False; (* selfopd can be a complicated expression *)
    t = Join[t, ireducTable[ selfopa, DO ]];
    MPVCFAST = True;
  ];

  (* Generalized operator for computing the T-matrix (i.e.,
  conductance and other transport properties) for arbitrary
  hybridization operators. 7.3.2013 *)

  (*
  If[calcopq["hyb_f"],
    AppendTo[t, {"dhyb_f"}];
    MPVCFAST = False;
    t = Join[t, ireducTable[ hybopf ]];
    MPVCFAST = True;
    ];
  *)

  (** Projection operators **)

  t = Join[t, mtSingletOp["P0", projector0  @ d[] ] ];
  t = Join[t, mtSingletOp["P1", projector1  @ d[] ] ];
  t = Join[t, mtSingletOp["P2", projector2  @ d[] ] ];
  t = Join[t, mtSingletOp["P0a", projector0  @ a[] ] ];
  t = Join[t, mtSingletOp["P1a", projector1  @ a[] ] ];
  t = Join[t, mtSingletOp["P2a", projector2  @ a[] ] ];

  (* This only works for (Q, S_z) basis, since these projectors do NOT
  transform as SU(2) multiplets. The same holds, obviously, for the
  components of the spin operator. *)
  
  If[isQSZ[] || isNONE[] || isSU2[] || isU1[] || isANYJ[] || isISOSZ[] || isSPU1[] 
    || isDBLSU2[] || isDBLISOSZ[] || isP[] || isPP[],
    (* Diagonal projectors. Recall that P^2=P *)
    t = Join[t, mtSingletOp["Pu", projectorUP @ d[] ] ];
    t = Join[t, mtSingletOp["Pd", projectorDO @ d[] ] ];

    t = Join[t, mtSingletOp["Pua", projectorUP @ a[] ] ];
    t = Join[t, mtSingletOp["Pda", projectorDO @ a[] ] ];

    t = Join[t, mtSingletOp["SXd",  spinx[d[]] ] ];    
    t = Join[t, mtSingletOp["SX2d", pow[spinx[d[]],2] ] ];
    t = Join[t, mtSingletOp["SYd",  spiny[d[]] ] ];
    t = Join[t, mtSingletOp["SY2d", pow[spiny[d[]],2] ] ];
    t = Join[t, mtSingletOp["SZd",  spinz[d[]] ] ];
    t = Join[t, mtSingletOp["SZ2d", pow[spinz[d[]],2] ] ];
    t = Join[t, mtSingletOp["SPd",  spinplus[d[]] ] ];
    t = Join[t, mtSingletOp["SP2d", pow[spinplus[d[]],2] ] ];
    t = Join[t, mtSingletOp["SMd",  spinminus[d[]] ] ];
    t = Join[t, mtSingletOp["SM2d", pow[spinminus[d[]],2] ] ];

    (* For model=KONDO *)
    t = Join[t, mtSingletOp["SX",  spinketbraX[SPIN] ] ];
    t = Join[t, mtSingletOp["SX2", pow[spinketbraX[SPIN],2] ] ];    
    t = Join[t, mtSingletOp["SY",  spinketbraY[SPIN] ] ];
    t = Join[t, mtSingletOp["SY2", pow[spinketbraY[SPIN],2] ] ];
    t = Join[t, mtSingletOp["SZ",  spinketbraZ[SPIN] ] ];
    t = Join[t, mtSingletOp["SZ2", pow[spinketbraZ[SPIN],2] ] ];
    t = Join[t, mtSingletOp["SP",  spinketbraP[SPIN] ] ];
    t = Join[t, mtSingletOp["SM",  spinketbraM[SPIN] ] ];
    t = Join[t, mtSingletOp["Ruu",  nc[ket[+1/2], bra[+1/2]] ] ];
    t = Join[t, mtSingletOp["Rdd",  nc[ket[-1/2], bra[-1/2]] ] ];
    t = Join[t, mtSingletOp["Rud",  nc[ket[+1/2], bra[-1/2]] ] ];
    t = Join[t, mtSingletOp["Rdu",  nc[ket[-1/2], bra[+1/2]] ] ];
    
    (* For model=TWOKONDO *)
    (* Operators sx1, etc. must be defined in the model definition itself. *)
    t = Join[t, mtSingletOp["S1X",  sx1 ] ];
    t = Join[t, mtSingletOp["S1X2", pow[sx1,2] ] ];    
    t = Join[t, mtSingletOp["S1Y",  sy1 ] ];
    t = Join[t, mtSingletOp["S1Y2", pow[sy1,2] ] ];
    t = Join[t, mtSingletOp["S1Z",  sz1 ] ];
    t = Join[t, mtSingletOp["S1Z2", pow[sz1,2] ] ];

    t = Join[t, mtSingletOp["S2X",  sx2 ] ];
    t = Join[t, mtSingletOp["S2X2", pow[sx2,2] ] ];    
    t = Join[t, mtSingletOp["S2Y",  sy2 ] ];
    t = Join[t, mtSingletOp["S2Y2", pow[sy2,2] ] ];
    t = Join[t, mtSingletOp["S2Z",  sz2 ] ];
    t = Join[t, mtSingletOp["S2Z2", pow[sz2,2] ] ];

    t = Join[t, mtSingletOp["S1S2", S1S2 ] ];
    
    (* Projection operators for two sites *)
    t = Join[t, mtSingletOp["P00", nc[projector0[d],  projector0[a]  ] ]];
    t = Join[t, mtSingletOp["P0u", nc[projector0[d],  projectorUP[a] ] ]];
    t = Join[t, mtSingletOp["P0d", nc[projector0[d],  projectorDO[a] ] ]];
    t = Join[t, mtSingletOp["P02", nc[projector0[d],  projector2[a]  ] ]];
    t = Join[t, mtSingletOp["Pu0", nc[projectorUP[d], projector0[a]  ] ]];
    t = Join[t, mtSingletOp["Puu", nc[projectorUP[d], projectorUP[a] ] ]];
    t = Join[t, mtSingletOp["Pud", nc[projectorUP[d], projectorDO[a] ] ]];
    t = Join[t, mtSingletOp["Pu2", nc[projectorUP[d], projector2[a]  ] ]];
    t = Join[t, mtSingletOp["Pd0", nc[projectorDO[d], projector0[a]  ] ]];
    t = Join[t, mtSingletOp["Pdu", nc[projectorDO[d], projectorUP[a] ] ]];
    t = Join[t, mtSingletOp["Pdd", nc[projectorDO[d], projectorDO[a] ] ]];
    t = Join[t, mtSingletOp["Pd2", nc[projectorDO[d], projector2[a]  ] ]];
    t = Join[t, mtSingletOp["P20", nc[projector2[d],  projector0[a]  ] ]];
    t = Join[t, mtSingletOp["P2u", nc[projector2[d],  projectorUP[a] ] ]];
    t = Join[t, mtSingletOp["P2d", nc[projector2[d],  projectorDO[a] ] ]];
    t = Join[t, mtSingletOp["P22", nc[projector2[d],  projector2[a]  ] ]];
    
    (* Projection operators for three sites *)
    If[calcopq["Pthree"],
     Print["Pthree"];
     Module[{ops, proj, comb, Pstr, ll},
      ops = {d, a, b};
      proj = {projector0, projectorUP, projectorDO, projector2};
      comb = Outer[#2[#1] &, ops, proj];
      Pstr[i_] := Switch[i, 1, "0", 2, "u", 3, "d", 4, "2"];
      ll = Table[{"P" <> Pstr[i] <> Pstr[j] <> Pstr[k],
                  nc[comb[[1, i]], comb[[2, j]], comb[[3, k]]]},
                  {i, 4}, {j, 4}, {k, 4}];
      ll = Flatten[ll, 2];
      ll = Map[(Print[First[#]]; mtdoSingletOp[First[#],Last[#]])&, ll];
      ll = Join @@ ll;
      t = Join[t, ll];            
     ];
    ];
  ];

 (** Correlations *)

  (* Spin on the dot *)
  t = Join[t, mtSingletOp["SdSd", spinspin[d[], d[]]] ];

  (* Scalar product of spin operators on dot and in electron liquid *)
  t = Join[t, mtSingletOp["SdSf", spinspin[d[], f[0]]] ];

  t = Join[t, mtSingletOp["n_dn_f", nc[number[d[]], number[f[0]]] ] ];

  t = Join[t, mtSingletOp["q_dq_f", nc[number[d[]]-1, number[f[0]]-1] ] ];

  (* ntot = \sum_i n_i *)
  t = Join[t, mtSingletOp["ntot", ntot] ];

  (* ntot^2 = (\sum_i n_i)^2 *)
  t = Join[t, mtSingletOp["ntot^2", ntot2] ];

  Module[{ntotUP, ntotDO, ntotUP2, ntotDO2},
    ntotUP = Plus @@ Map[number[#,UP]&, basopsdot];
    ntotDO = Plus @@ Map[number[#,DO]&, basopsdot];
    ntotUP2 = pow[ntotUP, 2];
    ntotDO2 = pow[ntotDO, 2];
    t = Join[t, mtSingletOp["ntotUP", ntotUP] ];
    t = Join[t, mtSingletOp["ntotDO", ntotDO] ];
    t = Join[t, mtSingletOp["ntotUP^2", ntotUP2] ];
    t = Join[t, mtSingletOp["ntotDO^2", ntotDO2] ];
  ];

  (* s^2 = (\sum_i S_i)^2 *)
  t = Join[t, mtSingletOp["s^2", ops2] ];

  (* even/odd combinations *)
  t = Join[t, mtSingletOp["neven",  neven]];
  t = Join[t, mtSingletOp["nodd",   nodd]];
  t = Join[t, mtSingletOp["neven^2",  pow[neven,2] ]];
  t = Join[t, mtSingletOp["nodd^2",   pow[nodd,2] ]];
  t = Join[t, mtSingletOp["s2even", s2even]];
  t = Join[t, mtSingletOp["s2odd",  s2odd]];
  t = Join[t, mtSingletOp["seso",   seso]];

  (* Phonon operators *)
  t = Join[t, mtSingletOp["acr",    phononplus[nph]]];
  t = Join[t, mtSingletOp["aan",    phononminus[nph]]];

  t = Join[t, mtSingletOp["nph",    phononnumber[nph] ] ];
  t = Join[t, mtSingletOp["nph^2",  pow[phononnumber[nph],2] ] ];
  
  t = Join[t, mtSingletOp["ph0",    nc[ket[0], bra[0]] ] ];
  t = Join[t, mtSingletOp["ph1",    nc[ket[1], bra[1]] ] ];
  t = Join[t, mtSingletOp["ph2",    nc[ket[2], bra[2]] ] ];
  t = Join[t, mtSingletOp["ph3",    nc[ket[3], bra[3]] ] ];
  t = Join[t, mtSingletOp["ph4",    nc[ket[4], bra[4]] ] ];
  t = Join[t, mtSingletOp["ph5",    nc[ket[5], bra[5]] ] ];
  t = Join[t, mtSingletOp["ph6",    nc[ket[6], bra[6]] ] ];
  t = Join[t, mtSingletOp["ph7",    nc[ket[7], bra[7]] ] ];
  t = Join[t, mtSingletOp["ph8",    nc[ket[8], bra[8]] ] ];
  t = Join[t, mtSingletOp["ph9",    nc[ket[9], bra[9]] ] ];

  displop = phononplus[nph] + phononminus[nph];
  t = Join[t, mtSingletOp["displ",    displop ] ];
  t = Join[t, mtSingletOp["displ^2",  pow[displop, 2] ] ];

  t = Join[t, mtSingletOp["displq",     nc[displop, number[d[]]-1]  ] ];
  t = Join[t, mtSingletOp["displhop0",  nc[displop, hop[d[], f[0]]] ] ];


 (* NEW, 18.9.2007: for computing spectral functions in models
     with electron-phonon coupling using the self-energy trick. *)

  If[calcopq["M_d"], (* M_sigma = (a+a^\dag) d_\sigma *)
    AppendTo[t, {"dM_d"}];
    opfn = nc[displop, d[#1, #2]] &;
    t = Join[t, ireducTable[ opfn ]];
  ];

  If[calcopq["N_d"], (* N_sigma = (a+a^\dag) f_{0\sigma} *)
    AppendTo[t, {"dN_d"}];
    opfn = nc[displop, f[#1, 0, #2]] &;
    t = Join[t, ireducTable[ opfn ]];
  ];

  (* Spectral function of an arbitrary linear combination of orbitals. 
     Specified as, for example, A_x(d[]+a[]) or A_x(d[]-2a[]+b[]), etc.*)

  Module[{lst, expr, koefs},
    lst = calcoplist["A_x"];
    If[Length[lst] != 0, MyPrint["A_x list: ", lst] ];
    Scan[ (
      AppendTo[t, { "d" <> First[#] }];
      expr = ToExpression @ Last[#];
      koefs = Map[Function[op,{Coefficient[expr, op], op}], basops];
      koefs = Select[koefs, First[#] != 0 &];
      t = Join[t, ireducTable[koefs] ]
      )&, lst];
  ];

  (* Local energies: Hinit will correspond to the initial NRG cluster, 
     while Hloc is the truly local part (i.e., without the hybridisation). *)
  Module[{op},
    op = H /. params;
    t = Join[t, mtSingletOp["Hinit", op]];
    op = op /. f[CR|AN, _, _]->0;
    t = Join[t, mtSingletOp["Hloc", op]];
  ];
    
  (* === Add site 'a' === *)
  If[NRDOTS >= 2,
  
    (* Occupancy of side dot *)
    t = Join[t, mtSingletOp["n_a", number[a[]] ] ];

    (* Occupancy of side dot *)
    t = Join[t, mtSingletOp["q_a", number[a[]]-1 ] ];
    
    (* Square of occupancy of side dot *)
    t = Join[t, mtSingletOp["n_a^2", pow[number[a[]], 2] ] ];
    
    (* Charge fluctuations of side dot, for isospin symmetry! *)
    t = Join[t, mtSingletOp["q_a^2", pow[number[a[]]-1, 2] ] ];

    (* Occupancy of side dot, double occupancy *)
    t = Join[t, mtSingletOp["n_a_ud", number[a[], DO] ~ nc ~ number[a[], UP] ] ];

    (* Spin on the dot *)
    t = Join[t, mtSingletOp["SaSa", spinspin[a[], a[]]] ];

    (* Correlation side dot <-> electron liquid *)
    t = Join[t, mtSingletOp["SaSf", spinspin[a[], f[0]]] ];

    (* Correlation side dot <-> embedded dot *)
    t = Join[t, mtSingletOp["SaSd", spinspin[a[], d[]]] ];

    (* s1s2sym = 1/N \sum_<ij> S_i . S_j, also for N >= 2 *)
    t = Join[t, mtSingletOp["s1s2sym", spinspinsymmetric[basopsdot] ] ];

    (* Isospin correlation side dot <-> embedded dot *)
    t = Join[t, mtSingletOp["IaId", isoiso[a[], d[]]] ];

    (* Out-of-diagonal number correlation *)
    t = Join[t, mtSingletOp["n_an_d", nc[number[a[]], number[d[]]] ] ];

    (* Out-of-diagonal number correlation *)
    t = Join[t, mtSingletOp["q_aq_d", nc[number[a[]]-1, number[d[]]-1] ] ];

    (* Out-of-diagonal number correlation with the Fermi sea *)
    t = Join[t, mtSingletOp["n_an_f", nc[number[a[]], number[f[0]]] ] ];

    If[calcopq["A_a"], (* Spectral density of side dot *)
      AppendTo[t, {"dA_a"}];
      t = Join[t, ireducTable[ a[] ] ];
    ];
      
    If[calcopq["A_even"], (* even combination *)  
      AppendTo[t, {"dA_even"}];
      t = Join[t, ireducTable[ {{1, d[]}, {1, a[]}} ] ];
    ];
    If[calcopq["A_odd"], (* odd combination *)
      AppendTo[t, {"dA_odd"}];
      t = Join[t, ireducTable[ {{1, d[]}, {-1, a[]}} ] ];
    ];

    If[calcopq["sigma_a"], (* Dynamic spin-susceptibility on the other dot *)
      AppendTo[t, {"tsigma_a"}];
      t = Join[t, ireducsigmaTable[ a[] ] ];
    ];

    If[calcopq["self_a"],
      AppendTo[t, {"dself_a"}];
      MPVCFAST = False; (* selfopa can be a complicated expression *)
      t = Join[t, ireducTable[ selfopa ]];
      MPVCFAST = True;
    ];
      
    If[calcopq["selfc_a"],
      AppendTo[t, {"dselfc_a"}];
      MPVCFAST = False; (* selfopd can be a complicated expression *)
      t = Join[t, ireducTable[ selfopca ]];
      MPVCFAST = True;
    ];

    (* Isospin operator, x component *)
    t = Join[t, mtSingletOp["pair_a", nc[a[CR,UP], a[CR, DO]] ] ];
    t = Join[t, mtSingletOp["Ix_a", isospinx[a[]] ] ];
    t = Join[t, mtSingletOp["Ix_a^2", pow[isospinx[a[]],2] ] ];

    (* Isospin operator, y component *)
    t = Join[t, mtSingletOp["Iy_a", isospiny[a[]] ] ];
    t = Join[t, mtSingletOp["Iy_a^2", pow[isospiny[a[]],2] ] ];

    (* Isospin operator, z component (Iz_a is equal to 1/2 q_a) *)
    t = Join[t, mtSingletOp["Iz_a", isospinz[a[]] ] ];
    t = Join[t, mtSingletOp["Iz_a^2", pow[isospinz[a[]],2] ] ];

    (** For SIDE and RING models **)
    (* This is also a-d bond operator in general *)
    hp = hop[a[], d[]];

    (* Hopping operator: n_e = 1/2(n_a+n_d+hop) *)
    (* n_o = 1/2(n_a+n_d-hop) *)
    t = Join[t, mtSingletOp["hop", hp] ];

    (* Square of hopping operator *)
    t = Join[t, mtSingletOp["hop^2", pow[hp,2] ] ];

    (* nc[n,hop] -> n_e^2 = 1/4(ntot^2+hop^2+2nhop) *)
    (* n_o^2 = 1/4(ntot^2+hop^2-2nhop *)
    t = Join[t, mtSingletOp["nhop", nc[ntot, hp ] ] ];
    
    (** Projector operators for two-sites **)

    (* For (Q, S) basis, Ptwo = Ppar + Pperp can be defined. *)
    t = Join[t, mtSingletOp["Ptwo",
      number[d[]] ~ nc ~ number[a[]] -
      2 hubbard[a[]] ~ nc ~ number[d[]] -
      2 hubbard[d[]] ~ nc ~ number[a[]] +
      4 hubbard[d[]] ~ nc ~ hubbard[a[]]
    ] ];

    (* The following only works for (Q, S_z) basis, since individual
    projectors do NOT transform as SU(2) multiplets. *)
    (* Mind the order: the first one is 'd', the second is 'a' *)
    If[isQSZ[] || isISOSZ[],
      t = Join[t, mtSingletOp["Puu", 
        nc[projectorUP[d[]], projectorUP[a[]] ] ]];
      t = Join[t, mtSingletOp["Pud", 
        nc[projectorUP[d[]], projectorDO[a[]] ] ]];
      t = Join[t, mtSingletOp["Pdu", 
        nc[projectorDO[d[]], projectorUP[a[]] ] ]];
      t = Join[t, mtSingletOp["Pdd", 
        nc[projectorDO[d[]], projectorDO[a[]] ] ]];
      
      (* Probability that two electrons with parallel spin occupy the dots *)
      t = Join[t, mtSingletOp["Ppar", 
        nc[projectorDO[d[]], projectorDO[a[]] ] +
        nc[projectorUP[d[]], projectorUP[a[]] ] 
      ]];

      (* Probability that two electrons with opposite spin occupy the dots *)
      t = Join[t, mtSingletOp["Pperp",
        nc[projectorUP[d[]], projectorDO[a[]] ] +
        nc[projectorDO[d[]], projectorUP[a[]] ]
      ]];

      t = Join[t, mtSingletOp["Sz",  spinspinz [d[], a[]] ]];
      t = Join[t, mtSingletOp["Sxy", spinspinxy[d[], a[]] ]];
      t = Join[t, mtSingletOp["Spm", spinspinpm[d[], a[]] ]];
      t = Join[t, mtSingletOp["Smp", spinspinmp[d[], a[]] ]];

      t = Join[t, mtSingletOp["SZa", spinz[a[]] ] ];
      t = Join[t, mtSingletOp["SZ2a", pow[spinz[a[]],2] ] ];
    ];
  ];


  (* === Add site 'b' === *)
  If[ NRDOTS >= 3,
    (* Occupancy of 2nd side dot *)
    t = Join[t, mtSingletOp["n_b", number[b[]] ] ];

    (* Occupancy of 2nd side dot *)
    t = Join[t, mtSingletOp["q_b", number[b[]]-1 ] ];
    
    (* Square of occupancy of the2nd side dot *)
    t = Join[t, mtSingletOp["n_b^2", pow[number[b[]], 2] ] ];

    (* Charge fluctuations of side dot, for isospin symmetry! *)
    t = Join[t, mtSingletOp["q_b^2", pow[number[b[]]-1, 2] ] ];

    (* Symmetrized fluctuations of side dots, for QSLR code *)
    t = Join[t, mtSingletOp["n_ab^2", 
      1/2 (pow[number[b[]], 2] + pow[number[b[]], 2]) ] ];

    (* Symmetrized fluctuations of side dots, for ISOLR code *)
    t = Join[t, mtSingletOp["q_ab^2", 
      1/2 (pow[number[b[]]-1, 2] + pow[number[b[]]-1, 2]) ] ];

    (* Out-of-diagonal number correlations *)
    t = Join[t, mtSingletOp["n_bn_d", nc[number[b[]], number[d[]]] ] ];
    t = Join[t, mtSingletOp["n_an_b", nc[number[a[]], number[b[]]] ] ];
    
    (* Correlation 2nd side dot <-> embedded dot *)
    t = Join[t, mtSingletOp["SbSd", spinspin[b[], d[]]] ];

    (* Correlation 1st side dot <-> 2nd side dot *)
    t = Join[t, mtSingletOp["SaSb", spinspin[a[], b[]]] ];

    (* Symmetrized correlation for LR symmetric triangular TQD *)
    t = Join[t, mtSingletOp["SabSd", 
      1/2 (spinspin[a[], d[]] + spinspin[b[], d[]]) ] ];

    (* Bonding orbital *)
    mob = 1/2 (a[CR, SIGMA] + Sqrt[2] d[CR, SIGMA] + b[CR, SIGMA]);
    (* Non-bonding orbital *)
    mon = 1/Sqrt[2] (a[CR, SIGMA] - b[CR, SIGMA]); (* Sign ??? *)
    (* Anti-bonding orbital *)
    moa = 1/2 (a[CR, SIGMA] - Sqrt[2] d[CR, SIGMA] + b[CR, SIGMA]);

    (* Symmetric combination of a-b *)
    mox = 1/Sqrt[2] (a[CR, SIGMA] + b[CR, SIGMA]); (* Sign ??? *)
    
    (* Occupancy of bonding orbital *)
    nmob = Sum[nc[mob, conj[mob]], {SIGMA, 0, 1}];
    nmob2 = pow[nmob, 2]; (* occupancy squared *)

    (* Occupancy of non-bonding orbital *)
    nmon = Sum[nc[mon, conj[mon]], {SIGMA, 0, 1}];
    nmon2 = pow[nmon, 2];

    (* Occupancy of anti-bonding orbital *)
    nmoa = Sum[nc[moa, conj[moa]], {SIGMA, 0, 1}];
    nmoa2 = pow[nmoa, 2];

    (* Occupancy of symmetric combination of a-b *)
    nmox = Sum[nc[mox, conj[mox]], {SIGMA, 0, 1}];
    nmox2 = pow[nmox, 2];
    
    t = Join[t, mtSingletOp["nmob", nmob] ];
    t = Join[t, mtSingletOp["nmob^2", nmob2] ];
    t = Join[t, mtSingletOp["nmon", nmon] ];
    t = Join[t, mtSingletOp["nmon^2", nmon2] ];
    t = Join[t, mtSingletOp["nmoa", nmoa] ];
    t = Join[t, mtSingletOp["nmoa^2", nmoa2] ];
    t = Join[t, mtSingletOp["nmox", nmox] ];
    t = Join[t, mtSingletOp["nmox^2", nmox2] ];

    (* Hop a-b *)
    hpab = hop[a[], b[]];
    t = Join[t, mtSingletOp["hopab", hpab] ];

    (* Square of hopping operator *)
    t = Join[t, mtSingletOp["hopab^2", pow[hpab,2] ] ];

    (* Hop b-d *)
    hpbd = hop[b[], d[]];
    t = Join[t, mtSingletOp["hopbd", hpbd] ];

    (* Square of hopping operator *)
    t = Join[t, mtSingletOp["hopbd^2", pow[hpbd,2] ] ];

    t = Join[t, mtSingletOp["SZb", spinz[b[]] ] ];
    t = Join[t, mtSingletOp["SZ2b", pow[spinz[b[]],2] ] ];

    t = Join[t, mtSingletOp["pair_b", nc[b[CR,UP], b[CR, DO]] ] ];
    t = Join[t, mtSingletOp["Ix_b", isospinx[b[]] ] ];
    t = Join[t, mtSingletOp["Ix_b^2", pow[isospinx[b[]],2] ] ];

    If[calcopq["A_b"], (* Spectral density of 2nd side dot *)
      AppendTo[t, {"dA_b"}];
      t = Join[t, ireducTable[ b[] ] ];
    ];
  ];


  (* === Add site 'e' === *)
  If[ NRDOTS >= 4,
    t = Join[t, mtSingletOp["n_e", number[e[]] ] ];
    t = Join[t, mtSingletOp["q_e", number[e[]]-1 ] ];
    t = Join[t, mtSingletOp["n_e^2", pow[number[e[]], 2] ] ];
    t = Join[t, mtSingletOp["q_e^2", pow[number[e[]]-1, 2] ] ];
    t = Join[t, mtSingletOp["n_an_e", nc[number[a[]], number[e[]]] ] ];
    t = Join[t, mtSingletOp["n_bn_e", nc[number[b[]], number[e[]]] ] ];
    t = Join[t, mtSingletOp["n_dn_e", nc[number[d[]], number[e[]]] ] ];
    t = Join[t, mtSingletOp["SaSe", spinspin[a[], e[]]] ];
    t = Join[t, mtSingletOp["SbSe", spinspin[b[], e[]]] ];
    t = Join[t, mtSingletOp["SdSe", spinspin[d[], e[]]] ];
    If[calcopq["A_e"],
      AppendTo[t, {"dA_e"}];
      t = Join[t, ireducTable[ e[] ] ];
    ];
  ];

  (* === Add site 'g' === *)
  If[ NRDOTS >= 5,
    t = Join[t, mtSingletOp["n_g", number[g[]] ] ];
    t = Join[t, mtSingletOp["q_g", number[g[]]-1 ] ];
    t = Join[t, mtSingletOp["n_g^2", pow[number[g[]], 2] ] ];
    t = Join[t, mtSingletOp["q_g^2", pow[number[g[]]-1, 2] ] ];
    t = Join[t, mtSingletOp["SdSg", spinspin[d[], g[]]] ];
    If[calcopq["A_g"],
      AppendTo[t, {"dA_g"}];
      t = Join[t, ireducTable[ g[] ] ];
    ];
  ];

  texportable = t;
];

MyPrint["operators.m done"];

texportable
