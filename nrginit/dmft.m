(*
    NRG Ljubljana -- dmft.m -- Support for arbitrary energy-dependent
    hybridization functions and for Dynamic Mean Field Theory
*)
    
(*    

NOTES
=====

The hybridization function (dmftgamma) must be given as a function of k, k
being the normalized energy in the interval [-1:1].
Unless an external module is called (parameter "run"), we use the expression
defined by parameter "gamma" in block [dmft] of the parameter file.
IMPROTANT: All occurencies of 'eps' are replaced by 'k', while 'ch' is
the channel number 1,...,COEFCHANNELS.

Alternatively, dmftgamma may be given as a list of COEFCHANNELS arrays of
tabulated hybridization functions.

FORMERLY, the convention was that gamma at k=0 is equal to 1; the overall
hybridization strength scale was then set by the factor realGamma. In other
words, function gamma defined the "form" of the hybridization function, not
its strength (Example: gamma[a_] = 1 for flat bands). For DMFT applications,
it is more practical to simply use arbitrary hybridization function; the
overall strength of the hybridization is then set by the integrals
thetaCh[a]. IMPORTANT: in this case, don't forget to either set Gamma to 1 
or to override the definition of theta0Ch[a].

Recall that, by definition (a form factor!), theta0Ch[a_] = 2 for a flat
band which has rho=1/2. Thus the relation between theta0Ch and thetaCh
(integral of rho over [-1:1]) should then be theta0Ch[a_] = 2 thetaCh[a_].

NOTE: Campo's approach works very well for analytical expressions, but it is
less numerically stable for tabulated data due to the (1/k) integration
kernel as opposed to (k) in Yoshida case. In particular, it is significantly
more difficult to maintain the particle-hole symmetry due to its sensitivity
to small differences in negative and positive part of the hybridization
function. 
*)

MyPrint["dmft.m started"];

(* NOTE: the precision of the results is expanded to PREC, since 
   high precision is required in the Lanczos tridiagonalization. 
   The precision must be likewise expanded in tabulated hybridization
   functions to ensure high precision integration. *)

setprec[z_] := SetPrecision[z, PREC];

(* Check if x is a positive real number, abort if it is not. *)
dumpcheckpos[msg_, x_] := Module[{},
  MyPrintForm[msg <> "[``][``]=``", aa, m, c10 @ x];
  If[Element[x, Reals] =!= True,
    MyError["Coefficient is not a real number."];
  ];
  If[x < 0,
    MyError["Coefficient is negative."];
  ];
];
dumpcheck[msg_, x_] := Module[{},
  MyPrintForm[msg <> "[``][``]=``", aa, m, c10 @ x];
  If[Element[x, Reals] =!= True,
    MyError["Coefficient is not a real number."];
  ];
];

(* Interpolation order for tabulated hybridization functions. Mathematica
default is 3. A more conservative choice would be to use ORDER = 1 in order
to avoid possible surprises with high-order polynomial approximations. *)

ORDER = 3;
Interp[l_] := Interpolation[l, InterpolationOrder -> ORDER];

(* Process the hybridization function in the form of a tabulated list to
enable efficient and accurate evaluation of the integrals. All results are
returned as function of parameter 'kk'. *)
  
dmftinterpolating[l_List] := Module[{lprec, l3, x3, Tneg, Tpos, 
  F0, Fpos, Fneg, Fposmul, Fnegmul, Fposdiv, Fnegdiv, xmin, xmax,
  Tnegmul, Tposmul, Tnegdiv, Tposdiv},
  
  MyPrint["dmftinterpolating[], len=", Length[l]];  

  lprec = setprec[l];

  (* Interpolate around 0 *)
  (* ATTENTION: the polynomial function should give strictly positive values in
     the relevant interval !!!! *)
  l3 = Take[Sort[lprec, Abs[First[#1]] < Abs[First[#2]] &], 3];
  x3 = l3[[All,1]];
  xmin = Min[x3];
  xmax = Max[x3];
  F0 = InterpolatingPolynomial[l3, kk];

    (* xmin, xmax determine the regions where the hybridisation function is
       determined by tabulated values: [-1:xmin] for N, [xmax:+1] for P. *)
  
  (* Positive and negative parts of the hybridization function *)
  Tneg = Select[lprec, First[#] <= xmin &];
  Tpos = Select[lprec, First[#] >= xmax &];

  Fneg = Interp[Tneg][kk];
  Fpos = Interp[Tpos][kk];
  
  (* kk Delta(kk) *)
  Tnegmul = Map[{#[[1]], #[[2]] #[[1]]}&, Tneg];
  Tposmul = Map[{#[[1]], #[[2]] #[[1]]}&, Tpos];
  Fnegmul = Interp[Tnegmul][kk];
  Fposmul = Interp[Tposmul][kk];

  (* Delta(kk)/kk *)
  Tnegdiv = Map[{#[[1]], #[[2]]/#[[1]]}&, Tneg];
  Tposdiv = Map[{#[[1]], #[[2]]/#[[1]]}&, Tpos];
  Fnegdiv = Interp[Tnegdiv][kk];
  Fposdiv = Interp[Tposdiv][kk];
  
  setprec @ {F0, Fpos, Fneg, Fposmul, Fnegmul, Fposdiv, Fnegdiv, xmin, xmax}
];

getdmftgamma[] := Module[{},
  dmftgamma = 0; (* GLOBAL! *)
  dmftscdelta = 0; (* GLOBAL! *)

  (* A) Define hybridization using parameter gamma in [dmft] block ... *)
  If[paramexists["gamma", "dmft"],
    dmftgamma = param["gamma", "dmft"];
  ];

  (* B) ... or load an external module to determine "dmftgamma" and "dmftscdelta". *)
  If[paramexists["run", "dmft"],
    fn = param["run", "dmft"];
    Print["Setting hybridisation using script ", fn];
    loadmodule[fn];
  ];
  
  If[dmftgamma == 0,
    MyError["dmftgamma not defined. Check [dmft] block in param file!"];
  ];
];

(* Reality check :-) *)
checkisgammareal[] := Module[{},
  For[aa = 1, aa <= COEFCHANNELS, aa++,
    (* Avoid extremal points -1 and 1, where gamma may diverge. *)
    TBLSTEP = 0.05;
    tblN = Table[gammaN[aa], {kk, -0.99, 0, TBLSTEP}];
    tblP = Table[gammaP[aa], {kk, 0, 0.99, TBLSTEP}];
    
    tbl = Join[tblN, tblP];
    If[ (And @@ Map[Element[#, Reals]&, tbl]) =!= True,
      MyError["Check definition for gamma! tbl=", tbl];
    ];
  ];
];
  
getgamma[] := Module[{},
  (* A) dmftgamma is an expression given as a String *)
  If[Head[dmftgamma] == String,
    MyPrint["dmftgamma is string ", dmftgamma];
    TABULATED = False;
    dmftgammaexpr = ToExpression[dmftgamma];
    For[aa = 1, aa <= COEFCHANNELS, aa++,
      gammafn[aa] = setprec @ Simplify[ dmftgammaexpr /. {eps -> kk, ch -> aa} ];
      gammaP[aa] = gammafn[aa];
      gammaN[aa] = gammafn[aa];
    ];
  ];
    
  (* B) dmftgamma is a List of tables, one table per channel. *)
  If[Head[dmftgamma] == List && Length[dmftgamma] == COEFCHANNELS,
    MyPrint["dmftgamma is tabulated"];
    TABULATED = True;
    For[aa = 1, aa <= COEFCHANNELS, aa++,
      gammaTAB[aa] = dmftinterpolating @ dmftgamma[[aa]];
      gammaP[aa] = gammaTAB[aa] [[2]];
      gammaN[aa] = gammaTAB[aa] [[3]];
    ];
  ];

  checkisgammareal[];
];

checkisxxintnumeric[] := Module[{},
    If[!NumericQ[xxintN[aa] /. kk->-0.1], MyError["xxintN not numeric."]];
    If[!NumericQ[xxintP[aa] /. kk->0.1],  MyError["xxintP not numeric."]];
];

checkisxxint0numeric[] := Module[{},
  If[!NumericQ[xxint0[aa] /. kk -> 10^-16], MyError["xxint0 not numeric."]];
];

xxint[] := Module[{},
  xxintN[aa] = Integrate[xxN, kk]; (* Indefinite integrals! *)
  xxintP[aa] = Integrate[xxP, kk];
  checkisxxintnumeric[];
  If[TABULATED,
    xxint0[aa] = Integrate[xx0, kk];
    checkisxxint0numeric[];
  ];
];

(* Setup integration routines. Helper functions xxint* need to be defined! 
See xxint[]. *)
integrationroutines[] := Module[{},
    intN[aa][a_, b_] /; a <= xmin && b <= xmin := 
      (xxintN[aa] /. kk->b) - (xxintN[aa] /. kk->a);
    intP[aa][a_, b_] /; a >= xmax && b >= xmax := 
      (xxintP[aa] /. kk->b) - (xxintP[aa] /. kk->a);

    If[TABULATED,
      intN[aa][a_, b_] /; a <= xmin && b > xmin :=  
        (xxintN[aa] /. kk->xmin) - (xxintN[aa] /. kk->a) +
        (xxint0[aa] /. kk->b)    - (xxint0[aa] /. kk->xmin);
      intN[aa][a_, b_] /; a > xmin && b > xmin :=  
        (xxint0[aa] /. kk->b)    - (xxint0[aa] /. kk->a);
      intP[aa][a_, b_] /; a < xmax && b >= xmax :=  
        (xxint0[aa] /. kk->xmax) - (xxint0[aa] /. kk->a) +
        (xxintP[aa] /. kk->b)    - (xxintP[aa] /. kk->xmax);
      intP[aa][a_, b_] /; a < xmax && b < xmax :=  
        (xxint0[aa] /. kk->b)    - (xxint0[aa] /. kk->a);
    ];
];
    
(* Conventional discretization points. *)
km[m_] = LAMBDA^(1-Z-m);
km[0] = 1;

Module[{},
  getdmftgamma[];
  getgamma[];

  (* Perform (indefinite) integrations and calculate the overall integral theta *)
  For[aa = 1, aa <= COEFCHANNELS, aa++,
    If[TABULATED,
      {xx0, xxP, xxN, xmin, xmax} = gammaTAB[aa] [[{1,2,3,8,9}]],
    (* else *)
      xxP = xxN = gammafn[aa];
      xmin = xmax = 0;
    ];
    
    xxint[];
    integrationroutines[];

    (* Take real parts to avoid any possible spurious imaginary
    parts arising from the integration. *)
    thetaCh[aa] = Re @ setprec[ intN[aa][-1, 0]  + intP[aa][0, 1] ];
    MyPrintForm["thetaCh[``]=``", aa, c10 @ thetaCh[aa]];

    (* Recall: df[], dfminus[] are integrals of the hybridisation function. *)
    For[m = 0, m <= mMAX, m++,
      df[aa, m] = Re @ setprec @ intP[aa][km[m+1], km[m]];
      dfminus[aa, m] = Re @ setprec @ intN[aa][-km[m], -km[m+1]];
      dumpcheckpos["df+", df[aa, m]];
      dumpcheckpos["df-", dfminus[aa, m]];
    ];

    (* de[], deminus[] are the energies of representative states of the
    conduction band in the star representation. The definition depends on
    the discretization scheme used. For DY, it is (\int x Gamma(x) dx)/(\int
    Gamma(x) dx), while for DC it is (\int Gamma(x) dx)/(\int Gamma(x)/x
    dx). *)
      
    If[DY,
      If[TABULATED,
        {xxP, xxN} = gammaTAB[aa] [[{4,5}]];
        xx0 = gammaTAB[aa] [[1]] kk, (* Multiplied by kk! *)
      (* else *)
        xxP = xxN = gammafn[aa] kk
      ];

      xxint[];

      (* NOTE THE MINUS SIGN in the definition of deminus[]! *)
      For[m = 0, m <= mMAX, m++,
        de[aa, m] = Re @ setprec[ intP[aa][km[m+1], km[m]]/df[aa, m] ];
        deminus[aa, m] = 
          - Re @ setprec[ intN[aa][-km[m], -km[m+1]]/dfminus[aa, m] ];
      ];
    ];

    If[DC,
      If[TABULATED,
        {xxP, xxN} = gammaTAB[aa] [[{6,7}]];
        xx0 = gammaTAB[aa] [[1]]/kk, (* Divided by kk! *)
      (* else *)
        xxP = xxN = gammafn[aa]/kk
      ];
        
      xxint[];  

      For[m = 0, m <= mMAX, m++,
        de[aa, m] = Re @ setprec[ df[aa, m] / intP[aa][km[m+1], km[m]] ];
        deminus[aa, m] = 
          -Re @ setprec[ dfminus[aa, m] / intN[aa][-km[m], -km[m+1]] ];
      ];
    ];
      
    For[m = 0, m <= mMAX, m++,
      dumpcheckpos["de+", de[aa, m]];
      dumpcheckpos["de-", deminus[aa, m]];
    ];
  ]; (* LOOP over COEFCHANNELS *)
    
  (* bug trap *)
  de[___] := MyError["oops dmft e"];
  df[___] := MyError["oops dmft f"];
];

(* Superconductivity. *)
dosc[] := Module[{},
  If[!(Head[dmftscdelta] == List && Length[dmftscdelta] == COEFCHANNELS),
    MyError["dmft.m error: dmftscdelta must be tabulated!"];
  ];
  If[!DY, MyError["dmft.m error: Only discretization=Y is supported!"]];

  For[aa = 1, aa <= COEFCHANNELS, aa++,
    l1 = dmftgamma[[aa]];
    l2 = dmftscdelta[[aa]];
    If[Length[l1] != Length[l2],
      MyError["dmft.m error: gamma and scdelta should be defined on the same grid!"];
    ];
    l3 = l1 * l2;
    l3[[All, 1]] = l1[[All,1]];
    
    scdeltaTAB[aa] = dmftinterpolating @ l3;
    scdeltaP[aa] = scdeltaTAB[aa] [[2]];
    scdeltaN[aa] = scdeltaTAB[aa] [[3]];

    {xx0, xxP, xxN} = scdeltaTAB[aa] [[{1,2,3}]];

    xxint[];

    For[m = 0, m <= mMAX, m++,
      dg[aa, m] = Re @ setprec[ intP[aa][km[m+1], km[m]]/df[aa, m] ];
      dgminus[aa, m] = Re @ setprec[ intN[aa][-km[m], -km[m+1]]/dfminus[aa, m] ];
      dumpcheck["dg+", dg[aa, m]];
      dumpcheck["dg-", dgminus[aa, m]];
    ];
  ];
];

If[isSPSU2[], 
  dosc[], 
(* else *)
  dg[___] := MyError["oops dmft g"];
];

MyPrint["dmft.m done"];
