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

(* Single Impurity Anderson Model *)
If[ MODEL == "SIAM" && (VARIANT != "" && VARIANT != "MAGFIELD"),
  def1ch[1];

  If[VARIANT == "SIAMKONDO",
    Hc = gammaPol hop[f[0], d[]] + Jkondo2 spinspin[f[0], d[]];
  ];

  If[VARIANT == "FM",
    Hc = gammaPol (hop[f[0], d[], UP] Sqrt[1+Pl] + hop[f[0], d[], DO] Sqrt[1-Pl]);
    H1 += Bz * spinz[d[]];
  ];

  If[VARIANT == "MAGFIELDXYZ",
    H1 += Bx * spinx[d[]] + By * spiny[d[]] + Bz * spinz[d[]];
  ];

  If[VARIANT == "MAGFIELDXYZnc",
    H1 += Bx * spinx[d[]] + By * spiny[d[]] + Bz * spinz[d[]];
    Hc = gammaPol hopphi[f[0], d[], UP, phi] k + 
         gammaPol hopphi[f[0], d[], DO, phi] (1/k);
  ];

  If[VARIANT == "INFU_EPS", 
    (* Note: we use eps rather than delta in this case. *)  
    H1 = eps*(number[d[]]-1);
    BASISRULE = "1-projector2[ d[] ]";
  ];

  If[VARIANT == "INFU_DELTA",
    BASISRULE = "1-projector2[ d[] ]";
  ];
  
  If[VARIANT == "INFU_GAP",
    H1 = eps*(number[d[]]-1); (* !! *)
    H1 = H1 + pairdelta pairing[d[]];
    BASISRULE = "1-projector2[ d[] ]";
  ];

  (* If it is preferrable to specify epsilon rather than delta *)
  If[VARIANT == "EPS",
    H1 = eps number[d[]] + U hubbard[d[]];
  ];

  If[VARIANT == "MAGFIELDEPS",
    H1 = eps number[d[]] + U hubbard[d[]] + B spinz[d[]];
  ];
    
  If[VARIANT == "DECOUPLED",
    Hc = 0;
  ];

  If[VARIANT == "HOPCHARGE",
    H00 = nc[hop[f[0], d[]], number[d[]]-1] + 
          nc[number[d[]]-1, hop[f[0], d[]]];
    Hc = gammaPol (hop[f[0], d[]] + beta H00);
  ];

  (* Jkondo is expressed in units of D. Therefore, rho*J=J/(2D) =Jkondo[here]/2. *)
  
  If[VARIANT == "KONDO", (* Pure Kondo model *)
    H1 = 0; 
    Hc = Jkondo spinspin[f[0], d[]];
    BASISRULE = "projector1[ d[] ]";
  ];

  If[VARIANT == "KONDOMAGFIELD", (* Pure Kondo model in magnetic field *)
    H1 = 0; 
    Hc = Jkondo spinspin[f[0], d[]] + B spinz[d[]];
    BASISRULE = "projector1[ d[] ]";
  ];
    
  If[VARIANT == "ISOKONDO", (* Pure isospin Kondo model *)  
    H1 = 0;
    Hc = Jcharge isoiso[f[0], d[]]; (* Analog of spinspin[] *)
    BASISRULE = "projector[d[], 6]";
  ];    

  If[VARIANT == "KONDO2", (* "Kondo model" with delta and U (no projection!) *)
    Hc = Jkondo spinspin[f[0], d[]];
  ];
    
  If[VARIANT == "ISOKONDO2",
    Hc = Jcharge isoiso[f[0], d[]]; (* Analog of spinspin[] *)
  ];    

  If[VARIANT == "TRUEKONDO", Module[{sz, sp, sm, oz, op, om, ss},
    sz = spinketbraZ[1/2];
    sp = spinketbraP[1/2];
    sm = spinketbraM[1/2];

    oz = nc[sz, 1/2(number[f[0], UP] - number[f[0], DO])];
    op = nc[sp, f[CR, 0, DO], f[AN, 0, UP]];
    om = nc[sm, f[CR, 0, UP], f[AN, 0, DO]];

    ss = oz + 1/2 (op + om) // Expand;
    MyVPrint[2, "TRUEKONDO ss=", ss];
    
    Hc = Jkondo ss;
    H1 = 0; (* Important, since there will be no d[] operator! *)
    BASISRULE = "projector[d[], 5]";
    MAKESPINKET = d[];
  ] ];
    
  If[VARIANT == "SW", (* Schrieffer-Wolff equivalent Kondo model *)
    H1 = 0;
    (* Note: rho=1/(2D) *)
    JJ = 4 realGamma/Pi (1/Abs[realU/2-realdelta] + 1/Abs[realU/2+realdelta]); 
    Hc = JJ spinspin[f[0], d[]];
    BASISRULE = "projector[d[], 5]";
  ];  

  If[VARIANT == "SWU", (* Schrieffer-Wolff equivalent Kondo model (+) on-site U *)
    JJ = 4 realGamma/Pi (1/Abs[realU/2-realdelta] + 1/Abs[realU/2+realdelta]); 
    Hc = JJ spinspin[f[0], d[]];
  ];  

  If[VARIANT == "HOLSTEIN",
    nph = ToExpression @ optionvalue["Nph"];
    MyPrint["nph=", nph];
    H1 = H1 + omega phononnumber[nph] 
            + gep nc[number[d[]]-1, phononplus[nph] + phononminus[nph]];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  If[VARIANT == "HOLSTEIN_MAGFIELD",
    nph = ToExpression @ optionvalue["Nph"];
    MyPrint["nph=", nph];
    H1 = H1 + omega phononnumber[nph] 
            + gep nc[number[d[]]-1, phononplus[nph] + phononminus[nph]];
    H1 = H1 + B spinz[d[]];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  If[VARIANT == "HOLSTEINSHIFT",
    nph = ToExpression @ optionvalue["Nph"];
    MyPrint["nph=", nph];
    H1 = H1 + omega phononnumber[nph] 
            + gep nc[number[d[]]-n0, phononplus[nph] + phononminus[nph]];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  If[VARIANT == "SZAA", (* Coupling of phonons to spin degree of freedom *)
    nph = ToExpression @ optionvalue["Nph"];
    MyPrint["nph=", nph];
    H1 = H1 + omega phononnumber[nph] 
            + gep nc[spinz[d[]], phononplus[nph] + phononminus[nph]];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  If[VARIANT == "COM", (* Coupling to centre-of-mass phonon mode *)
    nph = ToExpression @ optionvalue["Nph"];

    With[{displ = phononplus[nph] + phononminus[nph]},
      Hc = gammaPol nc[1 - ggg displ, hop[f[0], d[]] ];
    ];
    H1 = H1 + omega phononnumber[nph];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  If[VARIANT == "COM2", (* Centre-of-mass coupling to second-order *)
    MyPrint["SIAM/COM2"];
    nph = ToExpression @ optionvalue["Nph"];

    With[{displ = phononplus[nph] + phononminus[nph]},
      displ2 = pow[displ, 2];
      Print["displ2=", displ2];
      Hc = gammaPol nc[1 - ggg displ + 1/2 ggg^2 displ2, hop[f[0], d[]] ];
      Print["Hc=", Hc];
    ];
    H1 = H1 + omega phononnumber[nph];
    MAKEPHONON = 1; (* One phonon mode *)
  ];
    
  ClearAll[fastpow];  
  fastpow[op_, 0] := 1;
  fastpow[op_, 1] := op;
  fastpow[op_, i_Integer] /; i > 1 := 
    fastpow[op, i] = Expand[fastpow[op, -1 + i] \[CenterDot] op];

  If[VARIANT == "COMn", (* Centre-of-mass coupling to higher-order *)
    MyPrint["SIAM/COMn"];
    nph = ToExpression @ optionvalue["Nph"];
    If[!paramexists["order", "extra"],
      MyError["Define the expansion order for the exponential function!"];
    ];
    exporder = ToExpression @ param["order", "extra"];
    MyPrint["nph=", nph, " exporder=", exporder];

    With[{displ = phononplus[nph] + phononminus[nph]},
      expg = Series[Exp[-ggg x], {x, 0, exporder}] // Normal;
      MyPrint @ Timing[expg = expg /. 
        {x -> displ, x^n_ -> fastpow[displ, n]};];
      MyPrint @ Timing[expg = Expand[expg];];
      MyPrint["Length[expg]=", Length[expg]];
      rhs[i_] := rhs[i] = ket[i] + 
        Plus @@ Cases[expg, x_ nc[ket[a_], bra[i]] -> x ket[a]];
      
      nc[a___, EXPG, ket[i_]] := nc[a, rhs[i]];
      nc[a___, EXPG, vc[j___, ket[i_]]] := nc[a, rhs[i], vc[j]];
      ap[a___, EXPG, vc[j___, ket[i_]]] := ap[a, rhs[i], vc[j]];
      operatorQ[EXPG] ^= True;
      
      Hc = gammaPol nc[hop[f[0], d[]], EXPG];
      Print["Hc=", Hc];
    ];
    H1 = H1 + omega phononnumber[nph];
    MAKEPHONON = 1; (* One phonon mode *)
  ];

  (* Field affects spin down electrons only *)
  If[VARIANT == "MAGFIELDDOWN",
    checkQSZ[];
    H1 = H1 - B 1/2 number[d[], DO];
  ];
    
  If[VARIANT == "MAJORANA",  
    Hc = Hc + gammaPol (anodelta anomaloushop[f[0], d[]]
                     +  majdelta 
               (nc[f[CR,0,UP],d[CR,UP]] + nc[d[AN,UP],f[AN,0,UP]] +
                nc[f[CR,0,DO],d[CR,DO]] + nc[d[AN,DO],f[AN,0,DO]]) 
                );
  ];   

  If[VARIANT == "MAJORANAUPDO",  
    Hc = Hc + gammaPol (anodelta anomaloushop[f[0], d[]]
                     +  majdeltaup 
                     (nc[f[CR,0,UP],d[CR,UP]] + nc[d[AN,UP],f[AN,0,UP]])
                     + majdeltado
                     (nc[f[CR,0,DO],d[CR,DO]] + nc[d[AN,DO],f[AN,0,DO]])
                );
  ];   

  If[VARIANT == "MAJORANAUPDOB",
    Hc = Hc + gammaPol (anodelta anomaloushop[f[0], d[]]
                     +  majdeltaup 
                     (nc[f[CR,0,UP],d[CR,UP]] + nc[d[AN,UP],f[AN,0,UP]])
                     + majdeltado
                     (nc[f[CR,0,DO],d[CR,DO]] + nc[d[AN,DO],f[AN,0,DO]])
                   );
                     
    H1 = H1 + B spinz[d[]];
  ];   

  (* majdeltaup=0 corresponds to maj1up=maj2up=1. *)
  (* majdeltaup=1 corresponds to maj1up=2 maj2up=0. *)
  If[VARIANT == "MAJORANANEWMODEL",
    Hc = gammaPol (
      (nc[f[CR,0,UP],d[AN,UP]]+nc[d[CR,UP],f[AN,0,UP]]) (maj1up+maj2up)/2 
     +(nc[f[CR,0,UP],d[CR,UP]]+nc[d[AN,UP],f[AN,0,UP]]) (maj1up-maj2up)/2 
     +(nc[f[CR,0,DO],d[AN,DO]]+nc[d[CR,DO],f[AN,0,DO]]) (maj1do+maj2do)/2 
     +(nc[f[CR,0,DO],d[CR,DO]]+nc[d[AN,DO],f[AN,0,DO]]) (maj1do-maj2do)/2
    );
                     
    H1 = H1 + B spinz[d[]];
  ];   

  If[VARIANT == "MAJORANANEWMODELX",
    Hc = gammaPol (
      (nc[f[CR,0,UP],d[AN,UP]]+nc[d[CR,UP],f[AN,0,UP]]) (maj1up+maj2up)/2 
     +(nc[f[CR,0,UP],d[CR,UP]]+nc[d[AN,UP],f[AN,0,UP]]) (maj1up-maj2up)/2 
     +(nc[f[CR,0,DO],d[AN,DO]]+nc[d[CR,DO],f[AN,0,DO]]) (maj1do+maj2do)/2 
     +(nc[f[CR,0,DO],d[CR,DO]]+nc[d[AN,DO],f[AN,0,DO]]) (maj1do-maj2do)/2
    );
                     
    H1 = H1 + B spinz[d[]] + Bx spinx[d[]];
  ];   

  If[VARIANT == "ASSISTED", (* SIAM with the assisted hopping term *)
    Hc = gammaPol (hop[d[], f[0]] + 
      dd Sum[nc[number[d[], 1-SP], hop[d[], f[0], SP]], {SP, 0, 1}]);
  ];

  H = H0 + H1 + Hc;
];

(* Kondo model with arbitrary spin *)
If[ MODEL == "KONDO",
  def1ch[0];

  If[!paramexists["spin", "extra"],
    MyError["Define the spin of the impurity!"];
  ];
  SPIN = ToExpression @ param["spin", "extra"];
  MyPrint["SPIN=", SPIN];

  Module[{sz, sp, sm, sx, sy, oz, op, om, ss},
    sz = spinketbraZ[SPIN];
    sp = spinketbraP[SPIN];
    sm = spinketbraM[SPIN];
    sx = spinketbraX[SPIN];
    sy = spinketbraY[SPIN];    

    oz = nc[ sz, spinz[ f[0] ] ];
    op = nc[ sp, spinminus[ f[0] ] ];
    om = nc[ sm, spinplus[ f[0] ] ];

    ss = oz + 1/2 (op + om) // Expand;
    MyVPrint[2, "KONDO ss=", ss];
    
    Hc0 = Jkondo ss;
    Hc = Hc0 + B sz;

    If[VARIANT == "POTENTIAL",
      Hc = Hc + eps0 (number[f[0]]-1);
    ];
      
    If[VARIANT == "ANISO",
      sp2 = pow[sp, 2];
      sm2 = pow[sm, 2];
      Hc = Hc - anD pow[sz, 2] 
              - anB2/2 (sp2 + sm2)
              - anB4/2 (pow[sp2,2] + pow[sm2,2]);
    ];
    If[VARIANT == "ANISONEW",
      Hc = Hc + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]);
    ];

    If[VARIANT == "XXZANISONEW",
      Hc = JkondoZ oz + JkondoP 1/2 (op+om) +
           B sz // Expand;
      Hc = Hc + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]);
    ];
    
    If[VARIANT == "FIELDX",
      Hc = Hc + Bx sx;
    ];

    If[VARIANT == "FIELDXYZ",
      Hc = Hc0 + Bx sx + By sy + Bz sz;
    ];
    
    If[VARIANT == "GENERAL",
      Hc0 = JkondoZ oz + 1/2 JkondoP (op + om) // Expand;
      Haniso = anD pow[sz,2] + anE (pow[sx, 2] - pow[sy, 2]);
      Hc = Hc0 + Haniso + Bx sx + By sy + Bz sz;
    ];
      
    If[VARIANT == "DM",  
      Hc = Hc + JkondoDM (nc[sx, spiny[f[0]]] - nc[sy, spinx[f[0]]]);
    ];
  ];
    
  H = H0 + Hc;
   
  MAKESPINKET = SPIN; (* See QSZaddspinket[] function in SNEG *)
];

(* Combination of SIAM and KONDO models: impurity spin is side-coupled to the d-level. *)
If[ MODEL == "KLM",
  def1ch[1];

  If[!paramexists["spin", "extra"],
    MyError["Define the spin of the impurity!"];
  ];
  SPIN = ToExpression @ param["spin", "extra"];
  MyPrint["SPIN=", SPIN];

  Module[{sz, sp, sm, sx, sy, oz, op, om, ss},
    sz = spinketbraZ[SPIN];
    sp = spinketbraP[SPIN];
    sm = spinketbraM[SPIN];
    sx = spinketbraX[SPIN];
    sy = spinketbraY[SPIN];    

    oz = nc[ sz, spinz[ d[] ] ];
    op = nc[ sp, spinminus[ d[] ] ];
    om = nc[ sm, spinplus[ d[] ] ];

    ss = oz + 1/2 (op + om) // Expand;
    MyVPrint[2, "KONDO ss=", ss];

    (* B also on impurity site! *)
    Hc = Hc + B spinz[d[]];
    
    Hk = Jkondo ss;
    Hk = Hk + B sz + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]);
  ];
    
  H = H0 + H1 + Hc + Hk;
   
  MAKESPINKET = SPIN; (* See QSZaddspinket[] function in SNEG *)
];

(* Two-impurity one-channel Kondo model *)
If[ MODEL == "TWOKONDO",
  def1ch[0];

  If[!paramexists["spin1", "extra"],
    MyError["Define the spin of impurity 1!"];
  ];
  If[!paramexists["spin2", "extra"],
    MyError["Define the spin of impurity 2!"];
  ];
  SPIN1 = ToExpression @ param["spin1", "extra"];  
  SPIN2 = ToExpression @ param["spin2", "extra"];
  MyPrint["SPIN1=", SPIN1, "  SPIN2=", SPIN2];
    
  (* Identity operators *)
  id1 = spinketbraI[SPIN1, {1,0}];
  id2 = spinketbraI[SPIN2, {0,1}];
  
  (* First impurity *)
  sz1 = spinketbraZ[SPIN1, {1,0}] ~ nc ~ id2;
  sp1 = spinketbraP[SPIN1, {1,0}] ~ nc ~ id2;
  sm1 = spinketbraM[SPIN1, {1,0}] ~ nc ~ id2;
  sx1 = spinketbraX[SPIN1, {1,0}] ~ nc ~ id2;
  sy1 = spinketbraY[SPIN1, {1,0}] ~ nc ~ id2;

  (* Second impurity *)
  sz2 = spinketbraZ[SPIN2, {0,1}] ~ nc ~ id1;
  sp2 = spinketbraP[SPIN2, {0,1}] ~ nc ~ id1;
  sm2 = spinketbraM[SPIN2, {0,1}] ~ nc ~ id1;
  sx2 = spinketbraX[SPIN2, {0,1}] ~ nc ~ id1;
  sy2 = spinketbraY[SPIN2, {0,1}] ~ nc ~ id1;
  
  oz1 = nc[ sz1, spinz[ f[0] ] ];
  op1 = nc[ sp1, spinminus[ f[0] ] ];
  om1 = nc[ sm1, spinplus[ f[0] ] ];

  oz2 = nc[ sz2, spinz[ f[0] ] ];
  op2 = nc[ sp2, spinminus[ f[0] ] ];
  om2 = nc[ sm2, spinplus[ f[0] ] ];
  
  ss1 = oz1 + 1/2 (op1 + om1) // Expand;
  ss2 = oz2 + 1/2 (op2 + om2) // Expand;
  
  (* Basic model *)
  Hc = Jkondo1 ss1 + Jkondo2 ss2 + B (sz1+sz2); 
  
  If[VARIANT == "SS",
    S1S2 = spinketbraspinspin[{SPIN1, SPIN2}];
    Hc += J12 S1S2;
  ];

  If[VARIANT == "SSANISO",
    S1S2 = spinketbraspinspin[{SPIN1, SPIN2}];
    Hc += J12 S1S2;
    (* Anisotropy tensor. Usually anz=anD, anx=-any=anE. *)
    Hc += an1z pow[sz1, 2] + an1x pow[sx1, 2] + an1y pow[sy1,2];
    Hc += an2z pow[sz2, 2] + an2x pow[sx2, 2] + an2y pow[sy2,2];
  ];

  If[VARIANT == "FULL" || VARIANT == "FULLBY",
    Hc = Jkondo1 ss1 + Jkondo2 ss2; 
    Hc += B  (g1 sz1 + g2 sz2);
    Hc += Bx (g1 sx1 + g2 sx2);
    If[VARIANT == "FULLBY",
      Hc += By (g1 sy1 + g2 sy2);
    ];

    S1S2 = spinketbraspinspin[{SPIN1, SPIN2}];
    Hc += J12 S1S2;
    (* Anisotropy tensor. Usually anz=anD, anx=-any=anE. *)
    Hc += an1z pow[sz1, 2] + an1x pow[sx1, 2] + an1y pow[sy1,2];
    Hc += an2z pow[sz2, 2] + an2x pow[sx2, 2] + an2y pow[sy2,2];
  ];

  H = H0 + Hc;  
  
  MAKESPINKET = {SPIN1, SPIN2}; (* See QSZaddspinket[] function in SNEG *)
];

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

(* Single impurity Anderson model with a centre-of-mass phonon mode
modulating the hybridisation strength. *)
If[ MODEL == "ONE" && VARIANT == "COM",
  def2ch[1];
  nph = ToExpression @ optionvalue["Nph"];

  (* In single-impurity models the Sqrt[2] factors are indeed required for
  Gamma to be the ** total ** hybridisation strength! *)
  
  With[{displ = phononplus[nph] + phononminus[nph]},
    Hc = gammaPol/Sqrt[2] nc[1 - ggg displ, hop[f[0], d[]] ] +
         gammaPol/Sqrt[2] nc[1 + ggg displ, hop[f[1], d[]] ];
  ];
  H1 = H1 + omega phononnumber[nph];
  MAKEPHONON = 1; (* One phonon mode *)

  H = H0 + H1 + Hc;
  nnop[d[]] = -1; (* Override back *)
  If[isLR[],
  
    (* In addition, a -> -a and a^\dag -> -a^\dag, i.e. both bras and kets
    change sign! i.e. transform ket[i] into (-1)^i ket[i], since a ket is a
    (a^\dag)^i/i! \ket[0] state! This implies that only P=-1 states are
    coupled to the phonon mode, while all P=1 states are decoupled. *)
  
    lrchain = {f[0], d[], f[1]};
    lrextrarule = {ket[i__] -> (-1)^i ket[i]};
  ];
];

(* Same as ONE/KONDO, but without projecting out the excited states. 
   This version works fine with ISO. For Jkondo1==Jkondo2 we get a 2CK NFL GS.*)

If[ MODEL == "ONE" && VARIANT == "JKONDO",
  def2ch[1];

  Hc = Jkondo1 spinspin[f[0], d[]] + Jkondo2 spinspin[f[1], d[]];

  (* Strictly speaking, we should keep the U term. However, the 0 and doubly
     occupied states are decoupled from leads. So this is not a bug per se. *)

  H1 = 0;
  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)

  If[isLR[],
  (* Note: this only makes sense if Jkondo1 == Jkondo2. *)
    If[!(Jkondo1 == Jkondo2),
      MyError["Jkondo1 != Jkondo2"];
    ];
    lrchain = {f[0], d[], f[1]};
  ]; 
];

If[ MODEL == "ONE" && VARIANT == "SIAMKONDO",
  def2ch[1];

  Hc = gammaPol hop[f[0], d[]] + Jkondo2 spinspin[f[1], d[]];
  H = H0 + H1 + Hc;

  nnop[d[]] = -1; (* Override back *)
];

(* Superconductor-dot-superconductor system. cf. M.S. Choi, PRL 2004 *)
If[ MODEL == "ONE" && VARIANT == "SDS",  
  def2ch[1];

  (* For phi=0, Delta=0, we recover a single-channel problem with the usual hybridisation strength. *)
  Hc = gammaPol hop[f[0], d[]] Cos[phi/4] +
       gammaPol hop[f[1], d[]] Sin[phi/4];
  
  H = H0 + H1 + Hc;
];

(* DQD in ring geometry with a flux *)
(* Cf. AF Feiguin CA Busser, arxiv:1106.1464 *)
If[ MODEL == "ONE" && VARIANT == "FLUX",
  def2ch[1];

  (* For phi=0, we recover standard SIAM with hybridization Gamma *)
  Hc = gammaPol hop[f[0], d[]] Cos[phi] +
       gammaPol hop[f[1], d[]] Sin[phi];
  
  H = H0 + H1 + Hc;
       
  lrchain = {f[0], d[], f[1]};
];

(* The following model corrects a sign error in VARIANT == "SDS". *)
(* In practice, the difference is, however, very small. *)
If[ MODEL == "ONE" && VARIANT == "SDSCORRECT",  
  def2ch[1];

  (* For phi=0, Delta=0, we recover a single-channel problem
     with the usual hybridisation strength. *)
  Hc = gammaPol hop[f[0], d[]] Cos[phi/4] +
       gammaPol hop[f[1], d[]] Sin[phi/4];
       
  Print["H0=", H0];
       
  H = -H0 + H1 + Hc; (* Different sign!! *)
];

(* Fully general SD-QD-SD model. *)
If[ MODEL == "ONE" && VARIANT == "COMPLEXSDS",  
  def2ch[1];

  (* For phi=0, Delta=0, we recover a single-channel problem
     with the usual hybridisation strength. *)
  (* Mind the 1/Sqrt[2] factors!!! *)
  Hc = gammaPolch1/Sqrt[2] hopphi[f[0], d[], phi/2] +
       gammaPolch2/Sqrt[2] hopphi[f[1], d[], 0];
  
  H = H0 + H1 + Hc;
];

If[ MODEL == "ONE" && VARIANT == "SDSMAGFIELD",
  def2ch[1];

  (* For phi=0, Delta=0, we recoved a single-channel problem
     with the usual hybridisation strength. *)
  Hc = gammaPol hop[f[0], d[]] Cos[phi/4] +
       gammaPol hop[f[1], d[]] Sin[phi/4];

  H0 = H0 + B spinz[d[]];
  
  H = H0 + H1 + Hc;
];

If[ MODEL == "TWOKONDO2CH",
  def2ch[0];

  If[!paramexists["spin1", "extra"],
    MyError["Define the spin of impurity 1!"];
  ];
  If[!paramexists["spin2", "extra"],
    MyError["Define the spin of impurity 2!"];
  ];
  SPIN1 = ToExpression @ param["spin1", "extra"];  
  SPIN2 = ToExpression @ param["spin2", "extra"];
  MyPrint["SPIN1=", SPIN1, "  SPIN2=", SPIN2];
    
  (* Identity operators *)
  id1 = spinketbraI[SPIN1, {1,0}];
  id2 = spinketbraI[SPIN2, {0,1}];
  
  (* First impurity *)
  sz1 = spinketbraZ[SPIN1, {1,0}] ~ nc ~ id2;
  sp1 = spinketbraP[SPIN1, {1,0}] ~ nc ~ id2;
  sm1 = spinketbraM[SPIN1, {1,0}] ~ nc ~ id2;
  sx1 = spinketbraX[SPIN1, {1,0}] ~ nc ~ id2;
  sy1 = spinketbraY[SPIN1, {1,0}] ~ nc ~ id2;

  (* Second impurity *)
  sz2 = spinketbraZ[SPIN2, {0,1}] ~ nc ~ id1;
  sp2 = spinketbraP[SPIN2, {0,1}] ~ nc ~ id1;
  sm2 = spinketbraM[SPIN2, {0,1}] ~ nc ~ id1;
  sx2 = spinketbraX[SPIN2, {0,1}] ~ nc ~ id1;
  sy2 = spinketbraY[SPIN2, {0,1}] ~ nc ~ id1;
  
  (* Coupling of SPIN1 and SPIN2 to ch1: Jkondo1, Jkondo2 *)
  oz1 = nc[ sz1, spinz[ f[0] ] ];
  op1 = nc[ sp1, spinminus[ f[0] ] ];
  om1 = nc[ sm1, spinplus[ f[0] ] ];

  oz2 = nc[ sz2, spinz[ f[0] ] ];
  op2 = nc[ sp2, spinminus[ f[0] ] ];
  om2 = nc[ sm2, spinplus[ f[0] ] ];
  
  ss1 = oz1 + 1/2 (op1 + om1) // Expand;
  ss2 = oz2 + 1/2 (op2 + om2) // Expand;

  (* Coupling of SPIN1 and SPIN2 to ch2: Jkondo1ch2, Jkondo2ch2 *)
  oz1ch2 = nc[ sz1, spinz[ f[1] ] ];
  op1ch2 = nc[ sp1, spinminus[ f[1] ] ];
  om1ch2 = nc[ sm1, spinplus[ f[1] ] ];

  oz2ch2 = nc[ sz2, spinz[ f[1] ] ];
  op2ch2 = nc[ sp2, spinminus[ f[1] ] ];
  om2ch2 = nc[ sm2, spinplus[ f[1] ] ];
  
  ss1ch2 = oz1ch2 + 1/2 (op1ch2 + om1ch2) // Expand;
  ss2ch2 = oz2ch2 + 1/2 (op2ch2 + om2ch2) // Expand;
  
  (* Basic model *)
  Hc = Jkondo1 ss1 + Jkondo2 ss2 + B (sz1+sz2);
  
  If[VARIANT == "FULL" || VARIANT == "FULLBY",
    Hc =  Jkondo1    ss1    + Jkondo2    ss2;
    Hc += Jkondo1ch2 ss1ch2 + Jkondo2ch2 ss2ch2;
    Hc += B  (g1 sz1 + g2 sz2);
    Hc += Bx (g1 sx1 + g2 sx2);
    If[VARIANT == "FULLBY",
      Hc += By (g1 sy1 + g2 sy2);
    ];

    S1S2 = spinketbraspinspin[{SPIN1, SPIN2}];
    Hc += J12 S1S2;
    (* Anisotropy tensor. Usually anz=anD, anx=-any=anE. *)
    Hc += an1z pow[sz1, 2] + an1x pow[sx1, 2] + an1y pow[sy1,2];
    Hc += an2z pow[sz2, 2] + an2x pow[sx2, 2] + an2y pow[sy2,2];
  ];

  H = H0 + Hc;  
  
  MAKESPINKET = {SPIN1, SPIN2}; (* See QSZaddspinket[] function in SNEG *)
  
  hybopf0 = ( Chop @ Expand @ komutator[Hhyb /. params, f[#1, 0, #2]] )&;
  hybopf1 = ( Chop @ Expand @ komutator[Hhyb /. params, f[#1, 1, #2]] )&;
  
  If[isLR[],
    lrchain = {f[0], f[1]};
    lrextrarule = { ket[a_, b_] :> ket[b,a] };
  ];
  
];

(* As MODEL=TWOKONDO2CH, but for a single impurity ! *)

If[ MODEL == "ONEKONDO2CH",
  def2ch[0];

  If[!paramexists["spin1", "extra"],
    MyError["Define the spin of impurity 1!"];
  ];
  SPIN1 = ToExpression @ param["spin1", "extra"];  
  MyPrint["SPIN1=", SPIN1];
    
  (* First impurity *)
  sz1 = spinketbraZ[SPIN1];
  sp1 = spinketbraP[SPIN1];
  sm1 = spinketbraM[SPIN1];
  sx1 = spinketbraX[SPIN1];
  sy1 = spinketbraY[SPIN1];

  oz1 = nc[ sz1, spinz[ f[0] ] ];
  op1 = nc[ sp1, spinminus[ f[0] ] ];
  om1 = nc[ sm1, spinplus[ f[0] ] ];
  
  ss1 = oz1 + 1/2 (op1 + om1) // Expand;

  oz1ch2 = nc[ sz1, spinz[ f[1] ] ];
  op1ch2 = nc[ sp1, spinminus[ f[1] ] ];
  om1ch2 = nc[ sm1, spinplus[ f[1] ] ];
  
  ss1ch2 = oz1ch2 + 1/2 (op1ch2 + om1ch2) // Expand;
  
  If[VARIANT == "FULL" || VARIANT == "FULLBY",
    Hc =  Jkondo1    ss1;
    Hc += Jkondo1ch2 ss1ch2;
    Hc += B  (g1 sz1);
    Hc += Bx (g1 sx1);
    If[VARIANT == "FULLBY",
      Hc += By (g1 sy1);
    ];

    (* Anisotropy tensor. Usually anz=anD, anx=-any=anE. *)
    Hc += an1z pow[sz1, 2] + an1x pow[sx1, 2] + an1y pow[sy1,2];
  ];

  H = H0 + Hc;  
  
  MAKESPINKET = SPIN1; (* See QSZaddspinket[] function in SNEG *)

  lrchain = {f[0], f[1]};
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

(* Kondo model with arbitrary spin *)
If[ MODEL == "KONDO2CH",
  def2ch[0];

  If[!paramexists["spin", "extra"],
    MyError["Define the spin of the impurity!"];
  ];
  SPIN = ToExpression @ param["spin", "extra"];
  MyPrint["SPIN=", SPIN];

  Module[{sz, sp, sm, oz1, op1, om1, oz2, op2, om2, ss1, ss2},
    sx = spinketbraX[SPIN];
    sy = spinketbraY[SPIN];
    sz = spinketbraZ[SPIN];
    sp = spinketbraP[SPIN];
    sm = spinketbraM[SPIN];

    oz1 = nc[sz, 1/2(number[f[0], UP] - number[f[0], DO])];
    op1 = nc[sp, f[CR, 0, DO], f[AN, 0, UP]];
    om1 = nc[sm, f[CR, 0, UP], f[AN, 0, DO]];

    oz2 = nc[sz, 1/2(number[f[1], UP] - number[f[1], DO])];
    op2 = nc[sp, f[CR, 1, DO], f[AN, 1, UP]];
    om2 = nc[sm, f[CR, 1, UP], f[AN, 1, DO]];

    ss1 = oz1 + 1/2 (op1 + om1) // Expand;
    ss2 = oz2 + 1/2 (op2 + om2) // Expand;
    
    Hc = Jkondo1 ss1 + Jkondo2 ss2 + B sz;

    If[VARIANT == "XXZ",
      Hc = Jkondo1Z oz1 + Jkondo1P 1/2 (op1+om1) +
           Jkondo2Z oz2 + Jkondo2P 1/2 (op2+om2) +
           B sz // Expand;
    ];

    If[VARIANT == "ANISO",
      sp2 = pow[sp, 2];
      sm2 = pow[sm, 2];
      MyPrint["sp2=", sp2];
      MyPrint["sm2=", sm2];
      Hc = Hc - anD pow[sz, 2] 
              - anB2/2 (sp2 + sm2)
              - anB4/2 (pow[sp2,2] + pow[sm2,2]);
    ];
    If[VARIANT == "ANISONEW",
      sp2 = pow[sp, 2];
      sm2 = pow[sm, 2];
      Hc = Hc + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]);
    ];

    If[VARIANT == "XXZANISONEW",
      Hc = Jkondo1Z oz1 + Jkondo1P 1/2 (op1+om1) +
           Jkondo2Z oz2 + Jkondo2P 1/2 (op2+om2) +
           B sz // Expand;
      sp2 = pow[sp, 2];
      sm2 = pow[sm, 2];
      Hc = Hc + anD pow[sz, 2] + anE (pow[sx, 2] - pow[sy,2]);
    ];
          
  ];
    
  H = H0 + Hc;
   
  MAKESPINKET = SPIN; (* See QSZaddspinket[] function in SNEG *)

  If[isLR[],
    lrchain = {f[0], f[1]};
  ]; 
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
