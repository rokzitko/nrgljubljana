(*
  NRG Ljubljana - Hamiltonian definitions
*)

MyPrint["models started"];

(* NOTE: One can add new models to this file or to "custommodels.m". *)

(* Mapping of parameters as parsed from input files (prefixed by
   extra) to parameters as used in the Hamiltonian definition. *)

params = Join[params, {
    (* Hybridization to the 1st channel (!) *)
    gammaPol2 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma2) gammaA],

    (* Hybridization to the 2nd channel (!) *)
    gammaPol2to2 -> Sqrt[(1/Pi) thetaCh[2] (extraGamma2to2) gammaA],

    (* Symmetric definition for asymmetric couplings *)
    gammaPolch1 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma1) gammaA],
    gammaPolch2 -> Sqrt[(1/Pi) thetaCh[2] (extraGamma2) gammaA],
    gammaPolch3 -> Sqrt[(1/Pi) thetaCh[3] (extraGamma3) gammaA],

    (* Coupling terms: Alambda correction required! *)
    Jspin    ->  gammaA extraJspin,
    Jcharge  ->  gammaA extraJcharge,
    Jcharge1 ->  gammaA extraJcharge1,
    Jcharge2 ->  gammaA extraJcharge2,
    Jkondo   ->  gammaA extraJkondo,
    Jkondo1  ->  gammaA extraJkondo1,
    Jkondo2  ->  gammaA extraJkondo2,
    Jkondo3  ->  gammaA extraJkondo3,
    Jkondo1P ->  gammaA extraJkondo1P,
    Jkondo2P ->  gammaA extraJkondo2P,
    Jkondo1Z ->  gammaA extraJkondo1Z,
    Jkondo2Z ->  gammaA extraJkondo2Z,
    JkondoP  ->  gammaA extraJkondoP,
    JkondoZ  ->  gammaA extraJkondoZ,
    Jkondo1ch2  ->  gammaA extraJkondo1ch2,
    Jkondo2ch2  ->  gammaA extraJkondo2ch2,

    gep      -> extrag,
    dd       -> extrad
 }];
(* be careful not to introduce trailing comma in the list above! *)   

(* Scalar product of isospin operators for op1 and op2 *)

isoiso[op1_[j1___], op2_[j2___]] := Module[{isosp1, isosp2},
  isosp1 = isospin[op1[j1], nnop[op1[j1]]];
  isosp2 = isospin[op2[j2], nnop[op2[j2]]];
  inner[isosp1, isosp2]
];

(* d^\dag_\uparrow d^\dag_\downarrow + d_\downarrow d_\uparrow *)
(* Use the SPSU2 symmetry type for Hamiltonian with pairing terms! *)
pairing[op_] := 2 isospinx[op];


(* ===== Some examples ==== *)

(**** [ONE-CHANNEL] Two impurities ****)

(* Side-coupled quantum dots. Used in PRB 73, 035332 (2006) and PRB
   74, 241305(R) (2006). The additional site is denoted by 'a', while
   'd' is directly coupled to the conduction channel. *)

If[ MODEL == "SIDE",
  def1ch[2];

  nnop[a[]] = 0; (* Recall: d is -1 ! *)

  (* Inter-impurity Hamiltonian *)
  H1a = t hop[d[], a[]];

  (* Magnetic field. For QSZ only! *)
  If[VARIANT == "MAGFIELD",
    checkQSZ[];
    H1 = H1 + B (spinz[d[]] + spinz[a[]]);
  ];

  (* Various inter-impurity perturbation terms. *)
  If[VARIANT == "Uad",
    H1a += Uad nc[number[a[]]-1, number[d[]]-1]; (* Charge-coupling *)
  ];
  If[VARIANT == "Jad",
    H1a += Jad spinspin[a[], d[]]; (* Exchange interaction *)
  ];
  If[VARIANT == "Tad",
    H1a += Tad twohop[a[], d[]]; (* Two-electron hopping *)
  ];

  (* Unequal gate voltages *)
  If[VARIANT == "KAPPA",
    Ha = (delta-kappa) number[a[]] + U/2 pow[number[a[]]-1, 2];
    H1 = (delta+kappa) number[d[]] + U/2 pow[number[d[]]-1, 2];
  ];
  If[VARIANT == "KAPPA1",
    Ha = (delta-kappa) number[a[]] + U/2 pow[number[a[]]-1, 2];
  ];

  (* Non-interacting central impurity *)
  If[VARIANT == "DNI",
    H1 = delta number[d[]]; (* No U term! *)
  ];
    
  (* Non-interacting side-coupled impurity, separate epsilon *)
  If[VARIANT == "FANO",
    Ha = eps number[a[]];
  ];

  (* Interacting side-coupled impurity *)
  If[VARIANT == "FANOU1",
    Ha = eps number[a[]] + U1 hubbard[a[]];
  ];

  (* Interacting side-coupled impurity + dot-dot repulstion *)
  If[VARIANT == "FANOU1V",
    Ha = eps number[a[]] + U1 hubbard[a[]];
    Ha += V chargecharge[d[], a[]];
  ];

  If[VARIANT == "FANOU1FIELD",
    checkQSZ[];
    Ha = eps number[a[]] + U1 hubbard[a[]];
    Ha += B spinz[a[]];
    H1 += B spinz[d[]];
  ];
    
  H = H0 + H1 + Hc + Ha + H1a;
];

(* Two levels with inter-level repulsion *)
If[ MODEL == "TWOLEVELS",
  def1ch[2];

  H1 = delta1 number[d[]];
  Ha = delta2 number[a[]];
  H1a = Uad nc[number[d[]], number[a[]]];
  Hca = gammaPol hop[f[0], a[]];
  
  If[VARIANT == "ASYM",
    Hca = gammaPol2 hop[f[0], a[]];
  ];
  If[VARIANT == "ASYMHOP",
    Hca = gammaPol2 hop[f[0], a[]];
    H1a = H1a + tad hop[d[], a[]];
  ];
  (* NOTE: this parametrization corresponds to the Hamiltonian used
  in testing the NRG for master equations. *)
  (* NOTE: for spinless problems, spin UP components are kept. *)
  If[VARIANT == "TESTS",  
    Hca = gammaPol2 hop[f[0], a[]];
    H1a = Uad nc[number[d[],UP]-1/2, number[a[],UP]-1/2];
    H1a = H1a + tad hop[d[], a[]];
  ];  
    
  H = H0 + Hc + H1 + Ha + H1a + Hca;
];

(* Ring made out of dots 'a' and 'd', both coupled to the conduction
   channel. Used in PRB 74, 045312 (2006). *)

If[ MODEL == "RING", 
  def1ch[2];
  
  Had = 0;
  Hca = gammaPol hop[f[0], a[]];

  (* Energy shift *)
  Heps = -deltaeps number[a[]] + deltaeps number[d[]];

  (* Magnetic field. For QSZ only! *)
  If[VARIANT == "MAGFIELD",
    checkQSZ[];
    H1 = H1 + B (spinz[d[]] + spinz[a[]]);
  ];

  (* Magnetic field and inter-impurity hopping. *)
  (* For the entanglement paper PRB 74, 241305(R) (2006). *)
  If[VARIANT == "MAGFIELDtad",    
    checkQSZ[];
    H1 = H1 + B (spinz[d[]] + spinz[a[]]);
    Had = tad hop[a[], d[]];
  ];

  (* One dot couples via magnetic exchange, the other via isospin
     exchange. States 0-2 on dot 1 and 1 on dot 2 are projected
     out. *) 
  (* For the spin-charge separation paper PRB 74, 224411 (2006). *)
  If[VARIANT == "KONDOSPINISO",
    Hc  = Jspin spinspin[f[0], d[]];
    Hca = Jcharge isoiso[f[0], a[]];
    H1 = 0;
    Ha = 0;
    BASISRULE = "nc[ projector[d[], 5], projector[a[], 6] ]";
    ];
    
  (* One positive and one negative e-e interaction. This model is  
  studied in the spin-charge paper, PRB 74, 224411 (2006). *)
  If[VARIANT == "INVERSESIGN",
    Ha = delta number[a[]] - U/2 pow[number[a[]]-1, 2]; (* MINUS SIGN! *)
  ];
    
  (* Positive/negative U, hand-tuned hybridisation strengths. *)  
  If[VARIANT == "TWOGAMMAINVU", (* As above, with inverted U *)
    Hca = gammaPol2 hop[f[0], a[]];
    Ha = delta number[a[]] - U/2 pow[number[a[]]-1, 2]; 
  ];
  
  (* Perturbation: two-particle hopping between the impurities. Studied  
  in the paper on parallel quantum dots, PRB 74, 045312 (2006). *)
  If[VARIANT == "Tad",  
    Had = Tad twohop[a[], d[]];
  ];
  If[VARIANT == "tad",
    Had = tad hop[a[], d[]];
  ];
  If[VARIANT == "Uad",
    Had = Uad nc[number[a[]]-1, number[d[]]-1];
  ];
  If[VARIANT == "Jad",  
    Had = Jad spinspin[a[], d[]];
  ];
    
  (* Different hybridization strength Gamma2 for impurity 'a'. *)  
  (* Assimetric coupling of equal dots to the band... *)
  If[VARIANT == "ASIM", (* ... using a ratio parameter beta *)
    Hc  = gammaPol (1+beta) hop[f[0], d[]]; (* V_k = V_k^0 (1+beta) *)
    Hca = gammaPol (1-beta) hop[f[0], a[]];
  ];
  If[VARIANT == "TWOGAMMA", (* ... using an extra parameter Gamma2 *)
    Hca = gammaPol2 hop[f[0], a[]];
  ];

  H = H0 + H1 + Hc + Ha + Had + Hca + Heps;
];  


(**** [TWO-CHANNELS] No impurity ****)

If[ MODEL == "ZERO",
  def2ch[0];
  H = H0;

  (* Inter-band coupling term *)
  If[ VARIANT == "COUPLED",
    H = H0 + t hop[f[0], f[1]];
  ];
    
  lrchain = {f[0], f[1]};
];


(**** [TWO-CHANNELS] Single-impurity problems ****)

(* A test for LR symmetry. TO DO: When and where should this test be performed?? *)
checkgammasym[] := Module[{test},
  test = Abs[gammaPolCh[1] - gammaPolCh[2] /. params] < 10^-10;
  If[!test,
    MyError["Coupling to conduction channels is not symmetric."];
  ];
];
   

(* Single impurity Anderson model coupled to two conduction channels.
   Note that only symmetric channels couples to the impurity, while
   the antisymmetric channel is totally decoupled. This model can be
   used to study accidental (numerical) symmetry breaking by the NRG
   iteration. *)


If[ MODEL == "ONE" && VARIANT == "SIAM",
  def2ch[1];

  (* CHANGE 21.3.07 -> Sqrt[2] factor removed. *) 
  Hc = gammaPolCh[1] hop[f[0], d[]] + gammaPolCh[2] hop[f[1], d[]];
  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)

  lrchain = {f[0], d[], f[1]};
];

If[ MODEL == "ONE" && VARIANT == "SIAMFIELD",
  def2ch[1];

  Hc = gammaPolCh[1] hop[f[0], d[]] + gammaPolCh[2] hop[f[1], d[]];
  H = H0 + H1 + Hc + B spinz[d[]];

  nnop[d[]] = -1; (* Override back *)

  lrchain = {f[0], d[], f[1]};
];


(* As above, but with broken reflection symmetry. *)
If[ MODEL == "ONE" && VARIANT == "SIAMLR",
  def2ch[1];

  (* lr is the left-right asymmetry *)
  (* 
     Gamma_L = Gamma (1-lr)
     Gamma_R = Gamma (1+lr) 
  *)
     
  Hc = gammaPol Sqrt[1 - lr] hop[f[0], d[]] + 
       gammaPol Sqrt[1 + lr] hop[f[1], d[]];

  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];

(* As above, but using a multiplicative factor in hopping! *)
If[ MODEL == "ONE" && VARIANT == "SIAMk",
  def2ch[1];

  Hc = gammaPol k hop[f[0], d[]] + 
       gammaPol (1/k) hop[f[1], d[]];

  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];

(* As above, but using a multiplicative factor in hopping! *)
If[ MODEL == "ONE" && VARIANT == "SIAMkphi",
  def2ch[1];

  Hc = gammaPol k hopphi[f[0], d[], phi] + 
       gammaPol (1/k) hopphi[f[1], d[], -phi];

  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];

If[ MODEL == "ONE" && VARIANT == "SIAM12",
  def2ch[1];

  Hc = gammaPolch1 hop[f[0], d[]] + 
       gammaPolch2 hop[f[1], d[]];

  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];

If[ MODEL == "ONE" && VARIANT == "SIAM12FIELD",
  def2ch[1];

  Hc = gammaPolch1 hop[f[0], d[]] + 
       gammaPolch2 hop[f[1], d[]];

  H1 = H1 + B spinz[d[]];

  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];
   
(* Use this model to study the two-channel Kondo model and the
   associated 2CK non-Fermi liquid fixed point. *)
If[ MODEL == "ONE" && VARIANT == "KONDO",
  def2ch[1];

  Hc = Jkondo1 spinspin[f[0], d[]] + Jkondo2 spinspin[f[1], d[]];
  H1 = 0;

  H = H0 + H1 + Hc;

  BASISRULE = "projector1[ d[] ]";
  nnop[d[]] = -1; (* Override back *)

  lrchain = {f[0], d[], f[1]};
];

(**** [THREE-CHANNELS] One-impurity problems ****)

If[ MODEL == "3CH1" && VARIANT == "KONDO",
  def3ch[1];

  Hc = Jkondo1 spinspin[f[0], d[]] + 
       Jkondo2 spinspin[f[1], d[]] +
       Jkondo3 spinspin[f[2], d[]];
  H1 = 0;

  H = H0 + H1 + Hc;

  BASISRULE = "projector1[ d[] ]";
  nnop[d[]] = -1; (* Override back *)
];


(**** [TWO-CHANNELS] Two-impurity problems ****)

(* Two Kondo impurities in series: this is 2IK model with a 2IK critical point
   separating the impurity singlet and the double Kondo phase. *)

If[ MODEL == "TWO" && (VARIANT == "DQDKONDO" || VARIANT == "JDQDKONDO"),
  def2ch[2];
  
  Hc = Jkondo1 spinspin[f[0], a[]] + Jkondo2 spinspin[f[1], d[]];
  
  Had = Jad spinspin[a[], d[]];

  If[ VARIANT == "DQDKONDO", 
    H1 = 0;
    Ha = 0;
    BASISRULE = "nc[ projector[d[], 5], projector[a[], 5] ]";
  ];

  (* JDQDKONDO variant does not project out 0 and double occupancy. *)
  If[ VARIANT == "JDQDKONDO",
    (* We define these for definiteness. The values are not important,
       we have isospin symmetry no matter how the signs are defined. *)
    (* See nrg/test_initial_symmetries_JDQDKONDO.nb *)
    nnop[d[]] = -1;
    nnop[a[]] = -1;
  ];

  H = H0 + H1 + Hc + Ha + Had;
];

(* Three impurities embedded in series between two conduction leads.
   Studied in cond-mat/0606287. *)
If[ MODEL == "THREE" && (VARIANT == "NEWTQD" || VARIANT == "NEWTQD3"),
  def2ch[3];

  Hc = gammaPolCh[1] hop[f[0], a[]] + gammaPolCh[2] hop[f[1], b[]];

  Had = tpp hop[a[], d[]];
  Hbd = tpp hop[b[], d[]];

  nnop[a[]] = -1;
  nnop[b[]] = -1;

  If[VARIANT === "NEWTQD3",
    Ha = delta1 number[a[]] + U1/2 pow[number[a[]]-1, 2];
    H1 = delta2 number[d[]] + U2/2 pow[number[d[]]-1, 2];
    Hb = delta3 number[b[]] + U3/2 pow[number[b[]]-1, 2];
    Had = tpp12 hop[a[], d[]];
    Hbd = tpp23 hop[b[], d[]];
  ];

  lrchain = {f[0], a[], d[], b[], f[1]};

  (* Careful: must be the last thing done! *)
  H = H0 + H1 + Hc + Ha + Hb + Had + Hbd;
];
  

(* cf. Galpin et al. PRB 81, 075437 (2010) *)
If[ MODEL == "SU4",
  def2ch[2];

  Hc = gammaPolCh[1] hop[f[0], d[]] + gammaPolCh[2] hop[f[1], a[]];

  H1 = (delta-Borb) number[d[]] + B spinz[d[]];
  Ha = (delta+Borb) number[a[]] + B spinz[a[]];

  Had = U/2 pow[number[d[]] + number[a[]] - phi, 2];
  
  If[VARIANT == "MIXING",
    Hc = (gammaPolCh[1] hop[f[0], d[]] + gammaPolCh[2] hop[f[1], a[]]) Sqrt[Cos[beta]] +
         (gammaPolCh[1] hop[f[0], a[]] + gammaPolCh[2] hop[f[1], d[]]) Sqrt[Sin[beta]];
  ];
    
  H = H0 + H1 + Ha + Hc + Had;

  lrchain = {f[0], d[], a[], f[1]};
];
