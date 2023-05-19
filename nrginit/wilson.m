MyPrint["Wilson chain"];

ClearAll[thetaCh]; (* Bug honey-pot *)

(* ---- Tridiagonalisation approach (parameter "tri" in param file):
old - direct use of the recursion relations
sc - as above, but extended for superconducting hosts with arbitrary DOS
     and even-frequency pairing function
sc2 - as above, but for fully general pairing function
cpp - tridiagonalisation performed in the C++ part of the code
orth - recursion relations + orthogonality requirements
none - don't output the coeffcient table
manual - load the discretization tables from a file
manual_namnu - as above, for the superconducting case
*)

TRI = paramdefault["tri", "old"];
If[TRI == "old",
  defaultprec = 1000;
  dothelanczos = dothelanczosold;
];
If[TRI == "sc",
  defaultprec = 1000;
  dothelanczos = dothelanczossc;
];
If[TRI == "sc2",
  defaultprec = 1000;
  dothelanczos = dothelanczossc2;
];
If[TRI == "orth",
  defaultprec = 50;
  dothelanczos = dothelanczosorth;
];
If[TRI == "cpp",
  defaultprec = 30; (* Should do... *)
  dothelanczos = dothelanczosold;
];
If[TRI == "none",
  defaultprec = 30;
  dothelanczos = dothelanczosold;
];
If[TRI == "nambu",
  defaultprec = 30; (* ?? *)
  dothelanczos = dothelanczosnambu;
];
If[TRI == "manual",
  defaultprec = 30; (* Should be enough *)
  dothelanczos = loaddiscretizationtables;
];
If[TRI == "manual" || TRI == "manual_nambu" || TRI == "manual_nambu_new",
  defaultprec = 30; (* Should be enough *)
  dothelanczos = loaddiscretizationtables;
];
If[option["GENERATE_TEMPLATE"],
  dothelanczos = None;
  Nmax = 0;
];

(* If not specified: determine the interface type according to the current implementation in the C++ part of the code,
as encoded in defaultchaintype[] *)
(* legacy = enforce legacy interface *)
(* matrix = enforce matrix interface *)
defaultchaintype[_] := "legacy";
WILSONCHAIN = paramdefault["wilsonchain", defaultchaintype[SYMTYPE]];

(* Use arbitrary precision arithmetics *)
PREC = paramdefaultnum["prec", defaultprec];
MyVPrint[2,"PREC=", PREC];
setpr[expr_] := SetPrecision[expr, PREC];

LAMBDA = setpr @ lambda;
Z = setpr @ z;

(* DISCNMAX: Number of the discretization intervals. If nrxi/DISCNMAX is
negative, compute as many xi constants as there will be iterations. *)
DISCNMAX = paramdefaultnum["nrxi", -1];

If[paramexists["Nmax"] && paramexists["Tmin"],
  MyError["Specify either Nmax or Tmin, not both!"]; ];

If[paramexists["Nmax"],
  Nmax = ToExpression @ param["Nmax"];
];

Tmin=0;

If[paramexists["T"] && paramexists["Tmin_ratio"],
  If[paramnum["T"] > 0 && paramnum["Tmin_ratio"] > 0,
    Tmin = paramnum["T"] * paramnum["Tmin_ratio"];
    MyPrint["Tmin_ratio ==> Tmin=", Tmin];
  ];
];

If[paramexists["Tmin"],
  Tmin = paramnum["Tmin"];
  MyPrint["Tmin=", Tmin];
];

If[Tmin > 0,
  Nmax = 0;
  While[SCALE[Nmax+1] >= Tmin, Nmax++];
  MyPrint["Tmin=", Tmin, " ==> Nmax=", Nmax];
];

If[DISCNMAX < 0,
  DISCNMAX = Nmax;
];

MyPrint["DISCNMAX=", DISCNMAX];
If[DISCNMAX < 0 || DISCNMAX >= 999, MyError["Error"] ];

(* mMAX is the number of u(..., m) values used in the Lanczos calculation.
   Note that an equal number of the hybridization function integrals must
   be computed! *)
mMAX = Max[{80, 2 DISCNMAX}]; (* **** WAS: DISCNMAX + 40 **** *)
If[paramexists["mMAX"], mMAX = ToExpression @ param["mMAX"]];
MyVPrint[2,"mMAX=", mMAX];
If[mMAX <= 0 || mMAX >= 999, MyError["Error."]];

(* Override DISCNMAX in the case where the full tridiagonalisation
   is performed in the C++ part of the code. *)
If[TRI == "cpp" || TRI == "none",
  (* NOTE: we need to compute Ninit diagonalisation coefficients. *)
  TRUEDISCNMAX = DISCNMAX;
  DISCNMAX = Ninit;
];

(**********************************)

(* df[a, m] are the integrals over energy of the hybridisation form
function. They are also often denoted as \gamma_m (see Eq. (60) in the
Campo-Oliveira paper). a is the channel index: 1, ..., COEFCHANNELS and m=0,...
the interval number. The coefficients depend only on the dispersion and on
the discretization points. For constant density of states,
df[m]=eps[m]-eps[m+1], i.e. the interval width. *)

(* theta is the energy integral of the hybridisation form function over all
discretization intervals, i.e. a sum over all df and dfminus. *)

If[BAND == "flat",
  eps[a_, 0] = 1;
  eps[a_, m_] = LAMBDA^(-Z-m+1);
  df[a_, m_] := eps[a,m] - eps[a,m+1];
  dfminus = df; (* p-h symmetric *)
  thetaCh[a_] = 2; (* Denoted as {\bar \gamma} in Campo,Oliveira paper. *)
];

If[BAND == "cosine",
  df[a_, 0] = (-(Sqrt[-1 + LAMBDA^(2*Z)]/LAMBDA^(2*Z)) +
    ArcSec[LAMBDA^Z])/2;

  df[a_, m_] = (LAMBDA^(-m - Z)*(LAMBDA*Sqrt[1 - LAMBDA^(-2*(-1 + m +
    Z))] - Sqrt[1 - LAMBDA^(-2*(m + Z))] + LAMBDA^(m +
    Z)*(ArcCsc[LAMBDA^(-1 + m + Z)] - ArcCsc[LAMBDA^(m + Z)])))/2;

  dfminus = df; (* p-h symmetric *)
  thetaCh[a_] = Pi/2;

  (* NOTE: what enters the expression for the operator f_0 are the
  combinations df[]/thetaCh. Thus the normalization of the density of states
  is automatically taken care of. *)
];

(* xx *)
MyImport[fn_String, opts___] := Module[{l},
  l = Import[fn, opts];
  If[l === $Failed, MyError["Failed importing file ", fn]];
  l
];

ImportTable[fn_String] := Module[{l, dim},
  l = MyImport[fn];
  dim = Dimensions[l];
  If[Length[dim] == 2 && dim[[2]] == 2,
    (* Array of complex numbers masquerading as a matrix *)
    l = Map[First[#]+I Last[#]&, l];
  ];
  l
];

If[BAND == "dmft",   (* Run an external module! *)
  timestart["dmft"];
  loadmodule["dmft.m", True]; (* Must supply df, dfminus and thetaCh. *)
  timeadd["dmft"];
];

If[BAND == "nambu",
  0; (* Do nothing here *)
];

(*
   If band=manual, we handle the discretization manually outside initial.m.
   This should be performed before running nrginit, since some of the
   Wilson chain parameters are needed to define the initial NRG cluster.
   To be used with tri=none.
*)
If[BAND == "manual",
  Print["band=manual, importing theta, COEFCHANNELS=", COEFCHANNELS];
  Module[{fn, th},
    If[COEFCHANNELS == 1,
      fn = "theta.dat";
      th = Flatten[ImportTable[fn]] [[1]];
      thetaCh[a_] = th;
      Print["thetaCh[1]=", thetaCh[1]];
    ];
    If[COEFCHANNELS > 1,
      Do[
        fn = "theta" <> ToString[nrch] <> ".dat";
        th = Flatten[ImportTable[fn]] [[1]];
        thetaCh[nrch] = th;
        Print["thetaCh["<>ToString[nrch]<>"]=", thetaCh[nrch]],
          {nrch, COEFCHANNELS}];
      If[CHANNELS == 1 && POL2x2,
        Module[{th, th2},
          (* Special handling for matrix strucure in spin space! *)
          (* We need to take a matrix square root! *)
          th2 = {{thetaCh[1], thetaCh[3]}, {thetaCh[4], thetaCh[2]}};
          th = MatrixPower[th2, 1/2];
          thetaCh[1] = th[[1,1]]^2; (* Raise to power 2, because Sqrt is taken! *)
          thetaCh[2] = th[[2,2]]^2;
          thetaCh[3] = th[[1,2]]^2;
          thetaCh[4] = th[[2,1]]^2;
          Do[
            Print["thetaCh_fixed["<>ToString[nrch]<>"]=", thetaCh[nrch]],
            {nrch, COEFCHANNELS}];
        ];
      ];
    ];
  ];
];

(* A new version, where we deal directly with the hoppings between the impurity
   cluster and the first shell of the Wilson chain(s). *)
If[BAND == "manual_V",
  Print["band=manual_V, importing V, VDIM=", VDIM];
  Module[{fn, Vnn},
    If[VDIM == 1,
      fn = "V.dat";
      Vnn = Flatten[ImportTable[fn]] [[1]];
      V[_,_] = Chop[Vnn]; (* Chop[] is safe here, since Vnn is order 1 *)
      thetaCh[_] = Vnn^2; (* For backward compatibility - DO NOT USE! *)
      Print["V11=", V[1,1], " thetaCh[1]=", thetaCh[1]];
    ];
    If[VDIM > 1,
      Do[
        fn = "V" <> ToString[nrch1] <> ToString[nrch2] <> ".dat";
        Vnn = Flatten[ImportTable[fn]] [[1]];
        V[nrch1, nrch2] = Chop[Vnn];
        If[nrch1 == nrch2,
          thetaCh[nrch1] = Vnn^2; (* For backward compatibility - DO NOT USE! *)
        ];
        Print["V["<>ToString[nrch1]<>","<>ToString[nrch2]<>"]=", V[nrch1,nrch2]],
          {nrch1, VDIM},{nrch2, VDIM}];
    ];
  ];
];

hook[BAND];

(* This is called from maketable[]. *)
inittheta0ch[] := Module[{},
  (* Here we take into account the overall hybridization strength. *)
  If[DownValues[theta0Ch] === {} && realGamma >= 0,
    theta0Ch[a_] = bandrescale thetaCh[a] realGamma; (* YYY *)
    (* A single factor of bandrescale here! *)
  ];
  If[DownValues[theta0Ch] === {} && realGamma < 0,
    theta0Ch[a_] = bandrescale thetaCh[a]; (* YYY *)
  ];

  MyPrintForm["thetaCh=``", cstr10 /@ Array[thetaCh, COEFCHANNELS] ];
  MyPrintForm["theta0Ch=``", cstr10 /@ Array[theta0Ch, COEFCHANNELS] ];
  MyPrintForm["gammaPolCh=``", cstr10 /@ (Array[gammaPolCh, COEFCHANNELS] /. params) ];

  (* Backward compatibility for single-channel models. *)
  theta0 = theta0Ch[1];

  (* Bug trap *)
  Scan[ If[Negative @ theta0Ch[#], MyError["thetaCh negative"]]&, Range[COEFCHANNELS]];
];

(* COMMON CODE: the Lanczos diagonalization proper. de[] must be defined on
input. de[a, i] is the integral of the hybridisation function of channel 'a'
multiplied by energy, divided by the integral of the hybridization function
only: see Eq. (21) in Bulla et al. JPCM 9 10463 (1997). (Also denoted as
\Epsilon_n!) See also M. Sindel, PhD dissertation (2004), Appendix A:
"Derivation of the NRG-Equations" and Eqs. (15)-(17) in R. Bulla, Th.
Pruschke, A. C. Hewson, "Anderson impurity in pseudo-gap Fermi systems",
JPCM 9 10463 (1997). *)

(* Epsilon_n are sort of average energies in individual discretization
intervals!! *)

lanczosinit[] := Module[{},
  (* The algorithm used is a transcription of the procedure described by Kan
  Chen and C. Jayaprakash in "X-ray-edge singularities with nonconstant
  density of states: A renormalization-group approach", Phys. Rev. B 52,
  14436 (1995). Hereafter referenced as CJ. *)

  (* We are considering H_c=\sum_{m=0}^{\infty} \Lambda^{-m} (s_m a_m^\dag
  a_m-t_m b_m^\dag d_b). Since \Lambda^{-m} has been factored out, we divide
  coefficients de[], deminus[] by this factor. ds, dt are thus order 1 ! *)

  (* The minus sign in the definition of H_c is to be taken into account
  in the definition of deminus[] !! cf. dmft.m *)

  demem[a_, m_]      := demem[a, m]      = setpr @ de[a, m];
  deminusmem[a_, m_] := deminusmem[a, m] = setpr @ deminus[a, m];

  (* Diagonal matrix A in the CJ paper. *)
  diagA[a_] := diagA[a] =
    setpr @ Join[Table[de[a, m],       {m, 0, mMAX}],
                 Table[-deminus[a, m], {m, 0, mMAX}]];

  ClearAll[dzeta, xi2, xi, du, dv];

  (* Definition: f_n=\sum_m (u_{nm} a_m + v_{nm} b_m). Eq. (12) in CJ. *)

  (* The following lines are the same for all discretization schemes, since
  they are determined by the discretization points only. *)

  du[a_][-1, m_] = 0;
  dv[a_][-1, m_] = 0;

  du0[a_][0, m_] := Module[{fval},
    fval = df[a, m];
    If[fval < 0, MyError["negative hyb. integral df ", a, " ", m]];
    Sqrt[fval]/Sqrt[thetaCh[a]]
  ];
  dv0[a_][0, m_] := Module[{fval},
    fval = dfminus[a, m];
    If[fval < 0, MyError["negative hyb. integral dfminus ", a, " ", m]];
    Sqrt[fval]/Sqrt[thetaCh[a]]
  ];

  (* RECALL: du[0], dv[0] are essentially square roots of the hybridization
  function integrals. *)

  (* Useful test to see if df[], dfminus[] and thetaCh[] are self-consistently
     defined. We also rescale exactly to one! *)
  Module[{a, checksum},
    For[a = 1, a <= COEFCHANNELS, a++,
      checksum = Sum[du0[a][0,m]^2+dv0[a][0,m]^2, {m, 0, mMAX}];
      MyPrintForm["Discretization checksum [-1] (channel ``): ``",
            a, N[1-checksum, 10]];
      uvrescalefactor[a] = checksum;

      du[a_][0, m_] := du[a][0,m] = du0[a][0,m]/Sqrt[uvrescalefactor[a]];
      dv[a_][0, m_] := dv[a][0,m] = dv0[a][0,m]/Sqrt[uvrescalefactor[a]];

      (* Increase parameter DISCCHECKSUM in DMFT, since the effective
      hybridization will be somewhat oscillatory! *)
      DISCCHECKSUM = 10^-10;
      If[paramexists["discchecksum", "dmft"],
        DISCCHECKSUM = ToExpression @ param["discchecksum", "dmft"];
      ];
      If[Abs[checksum-1] > DISCCHECKSUM,
        MyError["Oops. Check your discretization code."];
      ];

    ];
    uvrescalefactor[a_] := MyError["Unknown band index."]; (* bug trap *)
  ];
];

dothelanczosnambu[] := Module[{},
  MyAssert[COEFCHANNELS == 2];
  (* For single-channel problems: *)
  (* a=1 spin up *)
  (* a=2 spin down *)

  dzeta[a_][m_] := 000;
  xi[a_][m_] := 000;

  sckappa[a_][n_] := 0;
  scdelta[a_][n_] := 000;

  (* Read thetaCh from file *)
];

dothelanczosold[] := Module[{},
  (* Eq. (18) in CJ. *)
  dzeta[a_][n_] :=
    dzeta[a][n] = Sum[(demem[a, m] du[a][n, m]^2 - deminusmem[a, m] dv[a][n, m]^2),
      {m, 0, mMAX}];

  (* Eq. (17) in CJ. *)
  xi[a_][n_] := xi[a][n] = Sqrt[xi2[a][n]];
  xi2[a_][-1] := 0;
  xi2[a_][n_] :=
    xi2[a][n] = Sum[(demem[a, m]^2 du[a][n, m]^2 + deminusmem[a, m]^2 dv[a][n, m]^2),
    {m, 0, mMAX}] - xi2[a][n - 1] - dzeta[a][n]^2;

  (* Eq. (15) in CJ. *)
  du[a_][n_, m_] :=
    du[a][n, m] = ((demem[a, m] - dzeta[a][n - 1]) du[a][n - 1, m] -
    xi[a][n - 2] du[a][n - 2, m])/xi[a][n - 1];

  (* Eq. (16) in CJ. *)
  dv[a_][n_, m_] :=
  dv[a][n, m] = ((-deminusmem[a, m] - dzeta[a][n - 1]) dv[a][n - 1, m] -
    xi[a][n - 2] dv[a][n - 2, m])/xi[a][n - 1];
];

(* Calculation of Wilson chain coefficients: ** superconducting host with
EVEN-frequency pairing function ** *)
(* Based on the work by Oliver Bodensiek, 2008 *)
dothelanczossc[] := Module[{},
  (* CONVENTION:
  f_{n,alpha} = U_{n alpha, m beta} a_{m,beta} + V_{n alpha, m beta} b_{m,beta},
  a_{m,beta} = U_{n alpha, m beta} f_{n,alpha},
  b_{m,beta} = U_{n alpha, m beta} f_{n,alpha}.
  *)

  duMAT[_, _, _][-1, _] = 0;
  dvMAT[_, _, _][-1, _] = 0;

  (* f_{0,alpha} = 1/\sqrt{theta_0} \sum_n (\gamma_n^+ a_{n,alpha}
                                           +\gamma_n^- b_{n,alpha}).
     Recall that du[], dv[] are defined to be gamma^{+-}_n/sqrt{theta_0},
     and that gamma^+=df, gamma^-=dfminus. See lanczsosinit[].
  *)

  duMAT[a_, alpha_, beta_][0, m_] := du[a][0, m] KroneckerDelta[alpha, beta];
  dvMAT[a_, alpha_, beta_][0, m_] := dv[a][0, m] KroneckerDelta[alpha, beta];

  (* deMAT are hybridisation matrices: diagonal contains hybridisation
  coefficients (de/deminus), out-of-diagonal coefficients are pairing
  coefficients. (2,2) component has inverted sign, since f f^\dag = -f^\dag
  f + 1. *)

  deplusMAT[a_, 1, 1][m_]  := deplusMAT[a, 1, 1][m]  = +setpr @ de[a, m];
  deplusMAT[a_, 2, 2][m_]  := deplusMAT[a, 2, 2][m]  = -setpr @ de[a, m];
  deplusMAT[a_, 1, 2][m_]  := deplusMAT[a, 1, 2][m]  = +setpr @ dg[a, m];
  deplusMAT[a_, 2, 1][m_]  := deplusMAT[a, 2, 1][m]  = +setpr @ dg[a, m];

  (* We define deminus[]/dgminus[] to be positive quantities, while in the
  Nambu matrix formalism we use absolute quantities, therefore in deminusMAT
  matrices we need to flip the sign of 1,1 and 2,2 components as compared to
  deplusMAT! *)

  deminusMAT[a_, 1, 1][m_] := deminusMAT[a, 1, 1][m] = -setpr @ deminus[a, m];
  deminusMAT[a_, 2, 2][m_] := deminusMAT[a, 2, 2][m] = +setpr @ deminus[a, m];
  deminusMAT[a_, 1, 2][m_] := deminusMAT[a, 1, 2][m] = +setpr @ dgminus[a, m];
  deminusMAT[a_, 2, 1][m_] := deminusMAT[a, 2, 1][m] = +setpr @ dgminus[a, m];

  (* dzeta is the on-site energy matrix; diagonal components are the
    on-site energy (note that the sign of 2,2 component is inverted), while
    the off-diagonal component is the on-site electron pairing. *)

  dzetaMAT[a_, alpha_, beta_][n_] := dzetaMAT[a, alpha, beta][n] =
    Sum[
      (duMAT[a, alpha, mu][n,m] deplusMAT [a, mu, nu][m] duMAT[a, beta, nu][n, m] +
       dvMAT[a, alpha, mu][n,m] deminusMAT[a, mu, nu][m] dvMAT[a, beta, nu][n, m]),
    {mu, 2}, {nu, 2}, {m, 0, mMAX}];

  (* Auxiliary matrix M *)

  dmMAT[a_, alpha_, beta_][n_, m_] := dmMAT[a, alpha, beta][n, m] =
    Sum[
      duMAT[a, beta, mu][n, m] deplusMAT[a, alpha, mu][m] -
      dzetaMAT[a, mu, beta][n] duMAT[a, mu, alpha][n, m],
    {mu, 1, 2}] - xiMAT[a, beta][n-1] duMAT[a, beta, alpha][n-1,m];

  (* Auxiliary matrix N *)

  dnMAT[a_, alpha_, beta_][n_, m_] := dnMAT[a, alpha, beta][n, m] =
    Sum[
      dvMAT[a, beta, mu][n, m] deminusMAT[a, alpha, mu][m] -
      dzetaMAT[a, mu, beta][n] dvMAT[a, mu, alpha][n, m],
    {mu, 1, 2}] - xiMAT[a, beta][n-1] dvMAT[a, beta, alpha][n-1,m];

  (* Hopping matrix, assumed to be diagonal, thus we only keep track of a
     single index alpha. *)

  xiMAT[a_, 1][n_] := xiMAT[a, alpha][n] = +Sqrt[xi2MAT[a, 1][n]];
  xiMAT[a_, 2][n_] := xiMAT[a, alpha][n] = -Sqrt[xi2MAT[a, 2][n]];
  xi2MAT[a_, alpha_][-1] := 0;
  xi2MAT[a_, alpha_][n_] := xi2MAT[a, alpha][n] =
    Sum[
      Sum[(dmMAT[a, mu, alpha][n, m])^2 + (dnMAT[a, mu, alpha][n, m])^2, {mu, 1, 2}],
    {m, 0, mMAX}];

  (* Recursion for U and V matrixes *)

  duMAT[a_, alpha_, beta_][n_, m_] := duMAT[a, alpha, beta][n, m] =
   (
    Sum[
       deplusMAT[a, beta, mu][m]      duMAT[a, alpha, mu  ][n - 1, m] -
       dzetaMAT[a, mu, alpha][n - 1]  duMAT[a, mu,    beta][n - 1, m],
    {mu, 2}] - xiMAT[a, alpha][n - 2] duMAT[a, alpha, beta][n - 2, m]
   )/xiMAT[a, alpha][n - 1];

 dvMAT[a_, alpha_, beta_][n_, m_] := dvMAT[a, alpha, beta][n, m] =
   (
    Sum[
       deminusMAT[a, beta, mu][m]     dvMAT[a, alpha, mu  ][n - 1, m] -
       dzetaMAT[a, mu, alpha][n - 1]  dvMAT[a, mu,    beta][n - 1, m],
    {mu, 2}] - xiMAT[a, alpha][n - 2] dvMAT[a, alpha, beta][n - 2, m]
   )/xiMAT[a, alpha][n - 1];

  (* Extract required components *)
  xi[a_][n_] := xi[a][n] = xiMAT[a, 1][n];
  dzeta[a_][n_] := dzeta[a][n] = dzetaMAT[a, 1, 1][n];
  scdelta[a_][n_] := scdelta[a][n] = dzetaMAT[a, 1, 2][n];
  sckappa[a_][n_] := sckappa[a][n] = 0;
];

(* Calculation of Wilson chain coefficients: ** superconducting host
with arbitrary (even/odd) pairing function frequency dependence ** *)
dothelanczossc2[] := Module[{},
  (* CONVENTION:
  f_{n,alpha} = U_{n alpha, m beta} a_{m,beta} + V_{n alpha, m beta} b_{m,beta},
  a_{m,beta} = U_{n alpha, m beta} f_{n,alpha},
  b_{m,beta} = U_{n alpha, m beta} f_{n,alpha}.
  *)

  duMAT[_, _, _][-1, _] = 0;
  dvMAT[_, _, _][-1, _] = 0;

  (* f_{0,alpha} = 1/\sqrt{theta_0} \sum_n (\gamma_n^+ a_{n,alpha}
                                           +\gamma_n^- b_{n,alpha}).
     Recall that du[], dv[] are defined to be gamma^{+-}_n/sqrt{theta_0},
     and that gamma^+=df, gamma^-=dfminus. See lanczsosinit[].
  *)

  duMAT[a_, alpha_, beta_][0, m_] := du[a][0, m] KroneckerDelta[alpha, beta];
  dvMAT[a_, alpha_, beta_][0, m_] := dv[a][0, m] KroneckerDelta[alpha, beta];

  (* deMAT are hybridisation matrices: diagonal contains hybridisation
  coefficients (de/deminus), out-of-diagonal coefficients are pairing
  coefficients. (2,2) component has inverted sign, since f f^\dag = -f^\dag
  f + 1. *)

  deplusMAT[a_, 1, 1][m_]  := deplusMAT[a, 1, 1][m]  = +setpr @ de[a, m];
  deplusMAT[a_, 2, 2][m_]  := deplusMAT[a, 2, 2][m]  = -setpr @ de[a, m];
  deplusMAT[a_, 1, 2][m_]  := deplusMAT[a, 1, 2][m]  = +setpr @ dg[a, m];
  deplusMAT[a_, 2, 1][m_]  := deplusMAT[a, 2, 1][m]  = +setpr @ dg[a, m];

  (* We define deminus[]/dgminus[] to be positive quantities, while in the
  Nambu matrix formalism we use absolute quantities, therefore in deminusMAT
  matrices we need to flip the sign of 1,1 and 2,2 components as compared to
  deplusMAT! *)

  deminusMAT[a_, 1, 1][m_] := deminusMAT[a, 1, 1][m] = -setpr @ deminus[a, m];
  deminusMAT[a_, 2, 2][m_] := deminusMAT[a, 2, 2][m] = +setpr @ deminus[a, m];
  deminusMAT[a_, 1, 2][m_] := deminusMAT[a, 1, 2][m] = +setpr @ dgminus[a, m];
  deminusMAT[a_, 2, 1][m_] := deminusMAT[a, 2, 1][m] = +setpr @ dgminus[a, m];

  (* dzeta is the on-site energy matrix; diagonal components are the
    on-site energy (note that the sign of 2,2 component is inverted), while
    the off-diagonal component is the on-site electron pairing. *)

  dzetaMAT[a_, alpha_, beta_][n_] := dzetaMAT[a, alpha, beta][n] =
    Sum[
      (duMAT[a, alpha, mu][n,m] deplusMAT [a, mu, nu][m] duMAT[a, beta, nu][n, m] +
       dvMAT[a, alpha, mu][n,m] deminusMAT[a, mu, nu][m] dvMAT[a, beta, nu][n, m]),
    {mu, 2}, {nu, 2}, {m, 0, mMAX}];

  (* Auxiliary matrix M *)

  dmMAT[a_, alpha_, beta_][n_, m_] := dmMAT[a, alpha, beta][n, m] =
    Sum[
      duMAT[a, beta, mu][n, m]  deplusMAT[a, alpha, mu][m] -
      dzetaMAT[a, mu, beta][n]  duMAT[a, mu, alpha][n, m] -
      xiMAT[a, mu, beta][n-1]   duMAT[a, mu, alpha][n-1,m],
    {mu, 1, 2}] ;

  (* Auxiliary matrix N *)

  dnMAT[a_, alpha_, beta_][n_, m_] := dnMAT[a, alpha, beta][n, m] =
    Sum[
      dvMAT[a, beta, mu][n, m] deminusMAT[a, alpha, mu][m] -
      dzetaMAT[a, mu, beta][n] dvMAT[a, mu, alpha][n, m] -
      xiMAT[a, mu, beta][n-1]  dvMAT[a, mu, alpha][n-1,m],
    {mu, 1, 2}] ;

  (* Hopping matrix, assumed to be diagonal, thus we only keep track of a
     single index alpha. *)

  xiMAT[_, _, _][-1] := 0;
  xiMAT[a_, 1, 2][n_] := xiMAT[a, 1, 2][n] = Sqrt[s2[a][n]];
  xiMAT[a_, 2, 1][n_] := xiMAT[a, 2, 1][n] = Sqrt[s2[a][n]];
  xiMAT[a_, 1, 1][n_] := xiMAT[a, 1, 1][n] = +Sqrt[t2[a][n]];
  xiMAT[a_, 2, 2][n_] := xiMAT[a, 2, 2][n] = -Sqrt[t2[a][n]];

  t2[a_][n_] := (t2s2[a][n] - s2mt2[a][n])/2;
  s2[a_][n_] := (t2s2[a][n] + s2mt2[a][n])/2;

  xiINV[a_][n_] := xiINV[a][n] = Inverse[
    {{xiMAT[a, 1, 1][n], xiMAT[a, 1, 2][n]},
     {xiMAT[a, 2, 1][n], xiMAT[a, 2, 2][n]}}
  ];

  xiINVMAT[a_, alpha_, beta_][n_] := xiINVMAT[a, alpha, beta][n] = xiINV[a][n] [[alpha, beta]];

  (* t2s2 = t^2+s^2, i.e. diagonal element of the hopping matrix squared *)
  t2s2[a_][n_] := t2s2[a][n] =
    Sum[
      Sum[(dmMAT[a, mu, 1][n, m])^2 +
          (dnMAT[a, mu, 1][n, m])^2,    {mu, 1, 2}],
      {m, 0, mMAX}];

  (* s2mt2 = s^2-t^2, i.e. the out-of-diagonal element of the (t.sigma_x.t) matrix,
      where t is the hopping matrix and sigma_x is the Pauli matrix x. *)
  s2mt2[a_][n_] := s2mt2[a][n] =
    Sum[dmMAT[a, 1, 1][n, m] dmMAT[a, 2, 2][n, m] +
        dmMAT[a, 2, 1][n, m] dmMAT[a, 1, 2][n, m] +
        dnMAT[a, 1, 1][n, m] dnMAT[a, 2, 2][n, m] +
        dnMAT[a, 2, 1][n, m] dnMAT[a, 1, 2][n, m],
    {m, 0, mMAX}];

  (* Recursion for U and V matrixes *)

  duMAT[a_, tau_, beta_][n_, m_] := duMAT[a, tau, beta][n, m] =
   Sum[ xiINVMAT[a, tau, alpha][n-1] *
    Sum[
       deplusMAT[a, beta, mu][m]      duMAT[a, alpha, mu  ][n - 1, m] -
       dzetaMAT[a, mu, alpha][n - 1]  duMAT[a, mu,    beta][n - 1, m] -
       xiMAT[a, mu, alpha][n - 2]     duMAT[a, mu,    beta][n - 2, m],
    {mu, 2}],
   {alpha, 2}];

  dvMAT[a_, tau_, beta_][n_, m_] := dvMAT[a, tau, beta][n, m] =
   Sum[ xiINVMAT[a, tau, alpha][n-1] *
    Sum[
       deminusMAT[a, beta, mu][m]     dvMAT[a, alpha, mu  ][n - 1, m] -
       dzetaMAT[a, mu, alpha][n - 1]  dvMAT[a, mu,    beta][n - 1, m] -
       xiMAT[a, mu, alpha][n - 2]     dvMAT[a, mu,    beta][n - 2, m],
    {mu, 2}],
   {alpha, 2}];

  (* Extract required components *)
  xi[a_][n_] := xi[a][n] = xiMAT[a, 1, 1][n];
  sckappa[a_][n_] := sckappa[a][n] = xiMAT[a, 1, 2][n];
  dzeta[a_][n_] := dzeta[a][n] = dzetaMAT[a, 1, 1][n];
  scdelta[a_][n_] := scdelta[a][n] = dzetaMAT[a, 1, 2][n];
];

dothelanczosorth[] := Module[{},
  (* Eq. (18) in CJ. *)
  dzeta[a_][n_] :=
  dzeta[a][n] = Sum[(demem[a, m] du[a][n, m]^2 - deminusmem[a, m] dv[a][n, m]^2),
    {m, 0, mMAX}];

  (* Eq. (17) in CJ. *)
  xi[a_][n_] := xi[a][n] = Sqrt[xi2[a][n]];
  xi2[a_][-1] := 0;
  xi2[a_][n_] := xi2[a][n] = Module[{},
    MyPrint[n];
    (* This is an appropriate point to perform the rescaling. See revise[] below. *)
    revise[a][n];
    Sum[(demem[a, m]^2 du[a][n, m]^2 + deminusmem[a, m]^2 dv[a][n, m]^2),
      {m, 0, mMAX}] - xi2[a][n - 1] - dzeta[a][n]^2
  ];

  (* Eq. (15) in CJ. *)
  duexpr[a_][n_, m_] := ((demem[a, m] - dzeta[a][n - 1]) du[a][n - 1, m] -
      xi[a][n - 2] du[a][n - 2, m])/xi[a][n - 1];
  du[a_][n_, m_] /; m >= Quotient[n, 2] := du[a][n, m] = duexpr[a][n, m];

  (* Eq. (16) in CJ. *)
  dvexpr[a_][n_, m_] := ((-deminusmem[a, m] - dzeta[a][n - 1]) dv[a][n - 1, m] -
    xi[a][n - 2] dv[a][n - 2, m])/xi[a][n - 1];
  dv[a_][n_, m_] /; m >= Quotient[n, 2] := dv[a][n, m] = dvexpr[a][n, m];

  (* Orthogonality equations, Eq.(22) in CJ *)
  eqlhs[a_][n_, j_] := Join[ Table[du[a][j, m], {m, 0, Quotient[n, 2]-1}],
                             Table[dv[a][j, m], {m, 0, Quotient[n, 2]-1}] ];
  eqrhs[a_][n_, j_] := Sum[-du[a][n, m] du[a][j, m], {m, Quotient[n, 2], mMAX}] +
                       Sum[-dv[a][n, m] dv[a][j, m], {m, Quotient[n, 2], mMAX}];

  solsys2[a_][n_] := solsys2[a][n] = Module[{lhs, rhs, sol},
    lhs = Table[eqlhs[a][n, j], {j, 0, 2*Quotient[n, 2]-1}];
    rhs = Table[eqrhs[a][n, j], {j, 0, 2*Quotient[n, 2]-1}];
    sol = LinearSolve[lhs][rhs];
    sol
  ];

  du[a_][n_, m_] /; m < Quotient[n, 2] :=
    du[a][n, m] = solsys2[a][n] [[ 1+m ]] ;
  dv[a_][n_, m_] /; m < Quotient[n, 2] :=
    dv[a][n, m] = solsys2[a][n] [[ 1+m+Quotient[n,2] ]];

  (* Vector U in the CJ paper. *)
  Uvec[a_][n_] := Join[Table[du[a][n, m], {m, 0, mMAX}],
                       Table[dv[a][n, m], {m, 0, mMAX}]];
  normn[a_][n_] := Uvec[a][n] . Uvec[a][n];

  (* Reduce roundoff error by renormalizing du and dv. *)
  revise[a_][n_] := Module[{alpha, m},
    alpha = 1/Sqrt[normn[a][n]];
    For[m = 0, m <= mMAX, m++,
      du[a][n,m] = setpr[ alpha du[a][n,m] ];
      dv[a][n,m] = setpr[ alpha dv[a][n,m] ];
    ];
    xi[a][n-1] = (diagA[a] Uvec[a][n-1]).Uvec[a][n];
    xi2[a][n-1] = (xi[a][n-1])^2;
  ];

  (* Trigger calculation at this point! This is important to ensure that
     the coefficient renormalizations are done properly. *)
  For[a = 1, a <= COEFCHANNELS, a++,
    Table[xi[a][i], {i, 0, DISCNMAX}];
  ];
  MyPrint["Lanczos done."];
];

(* Added 12.9.2012 *)
(* Removed discfaktor[n] on 21 Sep 2016 *)
loaddiscretizationtables[] := Module[{imp1,imp2},
  MyPrint["Loading discretization data from files."];
  If[COEFCHANNELS == 1,
    (* Load xi[a], dzeta[a]. *)
    imp1=Flatten[ImportTable["xi.dat"]];
    Print["xi=", Table[ xi[_][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
    (* Attention: zeta vs. dzeta ! *)
    imp2=Flatten[ImportTable["zeta.dat"]];
    Print["zeta=", Table[ dzeta[_][n] = imp2[[n+1]], {n,0,DISCNMAX} ]];
  ];
  If[COEFCHANNELS > 1,
    Do[
      Print["nrch=", nrch];
      imp1=Flatten[ImportTable["xi" <> ToString[nrch] <> ".dat"]];
      Print["xi=", Table[ xi[nrch][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
      imp2=Flatten[ImportTable["zeta" <> ToString[nrch] <> ".dat"]];
      Print["zeta=", Table[ dzeta[nrch][n] = imp2[[n+1]], {n,0,DISCNMAX} ]],
      {nrch, COEFCHANNELS}];
    If[RUNGS,
      Do[
        Print["nrch=", nrch];
        imp1=Flatten[ImportTable["xiR" <> ToString[nrch] <> ".dat"]];
        Print["xiR=", Table[ xiR[nrch][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
        imp2=Flatten[ImportTable["zetaR" <> ToString[nrch] <> ".dat"]];
        Print["zetaR=", Table[ zetaR[nrch][n] = imp2[[n+1]], {n,0,DISCNMAX} ]],
        {nrch, COEFCHANNELS}];
    ];
  ];
];

loadtablescdelta[] := Module[{imp1,imp2},
  MyPrint["Loading scdelta from a file."];
  If[COEFCHANNELS == 1,
    imp1=Flatten[ImportTable["scdelta.dat"]];
    Print["scdelta=", Table[ scdelta[_][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
  ];
  If[COEFCHANNELS > 1,
    Do[
      imp1=Flatten[ImportTable["scdelta" <> ToString[nrch] <> ".dat"]];
      Print["scdelta=", Table[ scdelta[nrch][n] = imp1[[n+1]], {n,0,DISCNMAX} ]],
      {nrch, COEFCHANNELS}];
  ];
];

loadtablesckappa[] := Module[{imp1,imp2},
  MyPrint["Loading sckappa from a file."];
  If[COEFCHANNELS == 1,
    imp1=Flatten[ImportTable["sckappa.dat"]];
    Print["sckappa=", Table[ sckappa[_][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
  ];
  If[COEFCHANNELS > 1,
    Do[
      imp1=Flatten[ImportTable["sckappa" <> ToString[nrch] <> ".dat"]];
      Print["sckappa=", Table[ sckappa[nrch][n] = imp1[[n+1]], {n,0,DISCNMAX} ]],
      {nrch, COEFCHANNELS}];
  ];
];

(* Use the following for debugging purposes. *)
discretizationChecks[] := Module[{},
  normalizationcheck[a_] := Module[{tab},
    tab = Table[1-Sum[du[a][n,m]^2, {n, 0, DISCNMAX}], {m, 0, mMAX}];
    MyPrintForm["Normalization: 1-sums_n u_{nm}^2 (channel ``)", a];
    MyPrint[cstr10[#]] & /@ tab;
  ];

  (* This orthogonality relation is non-trivial, in particular for
     the cases with a gap in the DOS !! *)
  orthscalar[a_][n_, j_] := Sum[du[a][n, m] du[a][j, m], {m, 0, mMAX}]  +
                            Sum[dv[a][n, m] dv[a][j, m], {m, 0, mMAX}];

  orthogonalitycheck[a_] := Module[{eq, tab, unknowns},
    MyPrint["Orthogonality (channel ``)", a];
    tab = Table[orthscalar[a][10, j], {j, 0, 10}];
    MyPrintForm["Orthogonality: sums_m u_{nm} u_{jm}+v_{nm} v_{jm} (channel ``)", a];
    MyPrint[cstr10[#]] & /@ tab;
  ];

  normalizationcheck[1];
  orthogonalitycheck[1];
];


(* See Yoshida, Whitaker, Oliveira, "Renormalization-group calculation of
excitation properties for impurity models", PRB 41 9403 (1990). Hereafter
referenced as YWO. *)

If[DY,
  (* Evaluated in lanczos-Yoshida.nb *)
  If[BAND == "flat",
    de[a_, 0] = (1 + LAMBDA^-Z)/2;
    de[a_, m_] = (1 + LAMBDA^-1)/2 LAMBDA^(1-Z-m);
    deminus = de; (* p-h symmetric *)
    dg[_, _] = dgminus[_, _] = paramdefaultnum["bcsgap2", "0"];
  ];

  (* Evaluated in lanczos-Yoshida.nb *)
  If[BAND == "cosine",
     de[a_, 0] = (-2*(-1 + LAMBDA^(2*Z))^(3/2))/(3*
       LAMBDA^Z*(Sqrt[-1 + LAMBDA^(2*Z)] - LAMBDA^(2*Z)*ArcSec[LAMBDA^Z]));

     de[a_, m_] = (2*LAMBDA^(m + Z)*((1 - LAMBDA^(-2*(-1 + m + Z)))^(3/2)
       - (1 - LAMBDA^(-2*(m + Z)))^(3/2)))/(3*(-(LAMBDA* Sqrt[1 -
       LAMBDA^(-2*(-1 + m + Z))]) + Sqrt[1 - LAMBDA^(-2*(m + Z))] +
       LAMBDA^(m + Z)*(-ArcCsc[LAMBDA^(-1 + m + Z)] + ArcCsc[LAMBDA^(m +
       Z)])));

     deminus = de; (* p-h symmetric *)
  ];
];

(* See Campo, Oliveira, "Alternative discretization in the numerical
renormalization-group method", PRB 72 104432 (2005). Hereafter referenced as
CO. *)

If[DC,
  (* Evaluated in lanczos-CampoOliveira.nb *)
  If[BAND == "flat",
    de[a_, m_] := (eps[a, m] - eps[a, m+1]) / Log[eps[a, m]/eps[a, m+1]];
    deminus = de; (* p-h symmetric *)
  ];

  (* Evaluated in lanczos-CampoOliveira.nb *)
  If[BAND == "cosine",
    de[a_, 0] = (Sqrt[-1 + LAMBDA^(2*Z)]/LAMBDA^(2*Z) -
     ArcSec[LAMBDA^Z])/(2*Sqrt[1 - LAMBDA^(-2*Z)] -
      2*(Z*Log[LAMBDA] + Log[1 + Sqrt[1 - LAMBDA^(-2*Z)]]));

    de[a_, m_] = -(LAMBDA^(-m - Z)*(-(LAMBDA*Sqrt[1 - LAMBDA^(-2*(-1 + m
      + Z))]) + Sqrt[1 - LAMBDA^(-2*(m + Z))] + LAMBDA^(m +
      Z)*(-ArcCsc[LAMBDA^(-1 + m + Z)] + ArcCsc[LAMBDA^(m + Z)])))/(2*(Sqrt[
      1 - LAMBDA^(-2*(-1 + m + Z))] - Sqrt[1 - LAMBDA^(-2*(m + Z))] +
      Log[LAMBDA] - Log[1 + Sqrt[1 - LAMBDA^(-2*(-1 + m + Z))]] + Log[1 +
      Sqrt[1 - LAMBDA^(-2*(m + Z))]]));

    deminus = de; (* p-h symmetric *)
  ];
];

(* New discretization scheme. See Zitko, Pruschke: "Energy resolution
and discretization artefacts in numerical renormalization group", (2008). *)
If[DZ,
  (* The case of constant density of states can be handled analytically. *)
  If[BAND == "flat",
    de[a_, 0] = (1-LAMBDA^-Z+Log[LAMBDA]-Z Log[LAMBDA])/Log[LAMBDA];
    de[a_, m_] = (eps[a, m] - eps[a, m+1]) / Log[eps[a, m]/eps[a, m+1]];
    deminus = de; (* p-h symmetric *)
  ]; (* BAND == "flat" *)

  (* Coside band, i.e. semi-circular DOS. Calculation handled by initial.m
     by numerically solving the discretization ODE using NDSolve. *)
  If[BAND == "cosine",
    Module[{z, rho0, norm, rho, in, omega,
            zmax, zfaktor, eqs, sol, fsol, g, tab},

      rho0[omega_] = Sqrt[1-omega^2];
      norm = NIntegrate[rho0[omega], {omega, -1, 1}];
      rho[omega_] = rho0[omega]/norm;
      in[omega_] = Integrate[rho[omega], omega];
      e[z_] := If[z<2, 1, lambda^(2-z)];
      alpha[z_] = -(in[e[z]]-in[e[z+1]]);

      rho0[omega_] := Sqrt[Abs[1-omega^2]] + 10^-5;
      (* THE NORMALIZATION IS IMPORTANT HERE! *)
      norm = NIntegrate[rho0[omega], {omega, -1, 1}];
      (* MyPrint["norm=", norm]; *)
      rho[omega_] := rho0[omega]/norm;

      zmax = 30; (* ideally mMAX+2 *)
      zfaktor = (1-lambda^-1)/Log[lambda];

      diff[z_] := -alpha[z]/lambda^(2-z) / (rho[f[z] lambda^(2-z)]);
      eqs = {f'[z] == Log[lambda] f[z] - diff[z], f[1] == 1/lambda};
      eqs = SetPrecision[eqs, 32];
      (* MyPrint[eqs]; *)

      sol = NDSolve[eqs, f, {z, 1, zmax}, WorkingPrecision -> 32];
      sol = sol[[1]];
      fsol = f /. sol;

      tab = Table[fsol[j], {j, 1, zmax}];
      (* Scan[MyPrint, tab]; *)

      (* Perform an extrapolation to smaller values! *)
      Eps[z_] := If[z<=zmax, fsol[z], fsol[zmax]] lambda^(2-z);

      MyPrint["Eps(zmax)/SCALE=", Eps[zmax]/(zfaktor lambda^(2-zmax))];
      tab = Table[Eps[j+1.0]/(zfaktor lambda^(2-j-1)), {j, 1, zmax}];
      (* Scan[MyPrint, tab]; *)
    ];
    de[a_, m_] := de[a, m] = Eps[1+m+z];
    deminus = de;
  ]; (* BAND == "cosine" *)

  (* Arbitrary asymmetric DOS, coefficients obtained by solving the discretization ODE
     using an external tool. This one supports multiple different conduction channels. *)
  If[BAND == "asymode",
    Module[{a, lfn, l, omega,
            rho0, xmax, tab, tabneg, rho, rhoneg,
            zfaktor, fsol, fsolneg},
     (* GLOBAL: intrho, intrhoneg *)

     For[a = 1, a <= COEFCHANNELS, a++,

      lfn = paramdefault["dos", "Delta.dat"];
      (* Multichannel support: append a channel number! *)
      If[a > 1, lfn = lfn <> ToString[a]];
      l = MyImport[lfn, "Table"];
      If[!MatrixQ[l], MyError["Loading DOS failed: not a matrix. ", lfn]];
      (* Rescaling [21 Sep 2016] *)
      l = Map[{First[#]/bandrescale, bandrescale Last[#]}&, l]; (* YYY *)
      (* Linear interpolation!! *)
      rho0 = Interpolation[l, InterpolationOrder -> 1];
      MyPrint["rho[0]=", rho0[0]];

      xmax = N @ paramdefaultnum["xmax", 30];

      tab = Select[l, Positive[ #[[1]] ]& ];
      tab = Sort[tab];
      tab = Prepend[tab, {0, tab[[1, 2]]}];
      tab = setpr @ tab;
      If[tab[[-1,1]] < 1.0, tab = Append[tab, {1, tab[[-1,2]]}]; ]; (* fix boundary *)
      rho = Interpolation[tab, InterpolationOrder -> 1]; (* Linear interpolation!! *)
      intrho[a][omega_] = Integrate[rho[omega], omega];

      eps[_, 0] = 1;
      (* Rescale here! *)
      If[!parambool["hardgap"],
        eps[_, m_] = LAMBDA^(-Z-m+1),
      (* else *)
        boundary = SetPrecision[paramnum["boundary"], PREC];
        eps[_, m_] = (1-boundary) LAMBDA^(-Z-m+1) + boundary;
      ];

      df[a_, m_] := df[a, m] =
        setpr[ intrho[a][eps[a,m]] - intrho[a][eps[a,m+1]] ];

      tabneg = Select[l, Negative[ #[[1]] ]& ];
      tabneg[[All,1]] = -tabneg[[All,1]]; (* Change sign! *)
      tabneg = Sort[tabneg];
      tabneg = Prepend[tabneg, {0, tab[[1,2]]}];
      tabneg = setpr @ tabneg;
      If[tabneg[[-1,1]] < 1.0, tabneg = Append[tabneg, {1, tabneg[[-1,2]]}]; ]; (* fix boundary *)
      rhoneg = Interpolation[tabneg, InterpolationOrder -> 1]; (* Linear interpolation!! *)
      intrhoneg[a][omega_] = Integrate[rhoneg[omega], omega];

      MyPrint["pos=", intrho[a][0.1]];
      MyPrint["neg=", intrhoneg[a][0.1]];

      dfminus[a_, m_] := dfminus[a, m] =
        setpr[ intrhoneg[a][eps[a,m]] - intrhoneg[a][eps[a,m+1]] ];

      thetaCh[a] = setpr[Integrate[   rho[x], {x, 0, 1}] +
                         Integrate[rhoneg[x], {x, 0, 1}]];
      MyPrint["theta=", thetaCh[a]];

      theta0Ch[a] = N @ thetaCh[a]; (* For DMFT only! *)

      solpath = paramdefault["solpath", ".."];
      lfn = solpath <> "/FSOL.dat";
      If[a > 1, lfn = lfn <> ToString[a]];
      fsol[a] = MyImport[lfn, "Table"];
      fsol[a] = Interpolation[fsol[a], InterpolationOrder -> 1];

      (* Perform an extrapolation to larger arguments x! *)
      Eps[a_, x_] := If[x <= xmax-1, fsol[a][x],
                                     fsol[a][xmax-1]] lambda^(2-x);

      zfaktor = (1-lambda^-1)/Log[lambda];

      (* Show results for error checking. *)
      tab = Table[{j, Eps[a, j+z]/(zfaktor lambda^(2-j-z))}, {j, 1, xmax}];
      Scan[MyPrint, tab];

      lfn = solpath <> "/FSOLNEG.dat";
      If[a > 1, lfn = lfn <> ToString[a]];
      fsolneg[a] = MyImport[lfn, "Table"];
      fsolneg[a] = Interpolation[fsolneg[a], InterpolationOrder -> 1];

      (* Perform an extrapolation to larger arguments x! *)
      Epsneg[a_, x_] := If[x <= xmax-1, fsolneg[a][x],
                                        fsolneg[a][xmax-1]] lambda^(2-x);

      (* Show results for error checking. *)
      tabneg = Table[{j, Epsneg[a, j+z]/(zfaktor lambda^(2-j-z))}, {j, 1, xmax}];
      Scan[MyPrint, tabneg];

      de[a_, m_] := de[a, m] = setpr @ Eps[a, 1+m+z];
      deminus[a_, m_] := deminus[a, m] = setpr @ Epsneg[a, 1+m+z];
    ]; (* loop over channels *)
   ]; (* Module *)
  ]; (* BAND == "asymode" *)

  (* Adaptable mesh approach. FSOL.dat and GSOL.dat need to be generated using
     an external tool. *)
  If[BAND == "adapt",
    Module[{lfn, l,
            rho0, xmax, tab, rho, intrho,
            gsol, gsolneg,
            tabneg, rhoneg, intrhoneg,
            zfaktor, fsol, fsolneg},

      lfn = paramdefault["dos", "Delta.dat"];
      l = MyImport[lfn];
      If[!MatrixQ[l], MyError["Loading DOS failed: not a matrix. ", lfn]];
      (* Rescaling [21 Sep 2016] *)
      l = Map[{First[#]/bandrescale, bandrescale Last[#]}&, l]; (* YYY *)
      (* Linear interpolation!! *)
      rho0 = Interpolation[l, InterpolationOrder -> 1];
      MyPrint["rho[0]=", rho0[0]];

      xmax = N @ paramdefaultnum["xmax", 30];

      tab = Select[l, Positive[ #[[1]] ]& ];
      tab = setpr @ tab;
      rho = Interpolation[tab, InterpolationOrder -> 1]; (* Linear interpolation!! *)
      intrho[omega_] = Integrate[rho[omega], omega];

      solpath = paramdefault["solpath", ".."];
      gsol = MyImport[solpath <> "/GSOL.dat"];
      gsol = setpr @ gsol;
      gsol = Interpolation[gsol, InterpolationOrder -> 1];

      eps[a_, 0] = 1;
      eps[a_, m_] := gsol[Z+m+1] LAMBDA^(-Z-m+1);
      df[a_, m_] := df[a, m] = setpr[ intrho[eps[a,m]] - intrho[eps[a,m+1]] ];

      tabneg = Select[l, Negative[ #[[1]] ]& ];
      tabneg[[All,1]] = -tabneg[[All,1]]; (* Change sign! *)
      tabneg = setpr @ tabneg;
      rhoneg = Interpolation[tabneg, InterpolationOrder -> 1]; (* Linear interpolation!! *)
      intrhoneg[omega_] = Integrate[rhoneg[omega], omega];

      MyPrint["pos=", intrho[0.1]];
      MyPrint["neg=", intrhoneg[0.1]];

      gsolneg = MyImport[solpath <> "/GSOLNEG.dat"];
      gsolneg = setpr @ gsolneg;
      gsolneg = Interpolation[gsolneg, InterpolationOrder -> 1];

      epsminus[a_, 0] = 1;
      epsminus[a_, m_] := gsolneg[Z+m+1] LAMBDA^(-Z-m+1);
      dfminus[a_, m_] := dfminus[a, m] =
        setpr[ intrhoneg[epsminus[a,m]] - intrhoneg[epsminus[a,m+1]] ];

      thetaCh[a_] = setpr[Integrate[rho[x], {x, 0, 1}] + Integrate[rhoneg[x], {x, 0, 1}]];
      MyPrint["theta=", thetaCh[1]];

      theta0Ch[a_] = N @ thetaCh[a]; (* For DMFT only! *)

      fsol = MyImport[solpath <> "/FSOL.dat"];
      fsol = Interpolation[fsol, InterpolationOrder -> 1];

      (* Perform an extrapolation to larger arguments x! *)
      Eps[x_] := If[x <= xmax-1, fsol[x], fsol[xmax-1]] lambda^(2-x);

      zfaktor = (1-lambda^-1)/Log[lambda];

      (* Show results for error checking. *)
      tab = Table[{j, Eps[j+z]/(zfaktor lambda^(2-j-z))}, {j, 1, xmax}];
      Scan[MyPrint, tab];

      fsolneg = MyImport[solpath <> "/FSOLNEG.dat"];
      fsolneg = Interpolation[fsolneg, InterpolationOrder -> 1];

      (* Perform an extrapolation to larger arguments x! *)
      Epsneg[x_] := If[x <= xmax-1, fsolneg[x], fsolneg[xmax-1]] lambda^(2-x);

      (* Show results for error checking. *)
      tabneg = Table[{j, Epsneg[j+z]/(zfaktor lambda^(2-j-z))}, {j, 1, xmax}];
      Scan[MyPrint, tabneg];
    ];
    de[a_, m_] := de[a, m] = setpr @ Eps[1+m+z];
    deminus[a_, m_] := deminus[a, m] = setpr @ Epsneg[1+m+z];
  ]; (* BAND == "adapt" *)


  (* See Hoeck and Schnack, Phys. Rev. B 87, 184408 (2014), Appendix A. *)
  (* See Eqs. (A12) and (A13). Added 15 Jun 2016. *)

  (* When we add local field in NRG Ljubljana, we typically write B spinz[d].
     The Hamiltonian term is g_imp mu_B B S_z = g_imp mu_B B (1/2) sigma,
     with sigma = +-1.
     Factor (1/2) is taken into account in spinz[]. Thus B corresponds to
     g_imp mu_B B. *)

  (* In Hoeck & Schnack, \epsilon_{k\sigma} = \epsilon_k + \sigma g_bulk mu_B B
     where, importantly, sigma = +-(1/2). Also mu = +-(1/2) in their equations.
     h is defined as h=g_bulk mu_B B, which is consistent with our convention for
     impurity B. *)

  (* In other words: both B and bulkh in NRG Ljubljana are defined as Zeeman energy.
     bulkh is thus the splitting between the bottom band edges for spin up and
     spin down electrons. For decoupled band, computing SZf0 we should thus find
     S_z(f0) = 1/2 (n_up - n_down). n_up = rho (W-h/2), n_down = rho (W+h/2).
     Thus S_z(f0) = -1/2 rho W h = -h/4, since rho=1/(2W). *)

  If[BAND == "flat_with_bulk_field",
    MyAssert[POLARIZED == True];
    MyAssert[SYMTYPE == "SPU1" || SYMTYPE == "QSZ"];

    (* Suffix denotes the h=0 values. *)
    eps0[a_, 0] = 1;
    eps0[a_, m_] = LAMBDA^(-Z-m+1);
    df0[a_, m_] := eps0[a,m] - eps0[a,m+1];
    dfminus0 = df0; (* p-h symmetric *)
    de0[a_, 0] = (1-LAMBDA^-Z+Log[LAMBDA]-Z Log[LAMBDA])/Log[LAMBDA];
    de0[a_, m_] = (eps0[a, m] - eps0[a, m+1]) / Log[eps0[a, m]/eps0[a, m+1]];
    deminus0 = de0; (* p-h symmetric *)

    bulkh = paramdefaultnum["bulkh", 0]; (* Bulk magnetic field in units of W. *)
    MyPrint["bulkh=", bulkh];
    bulkh = setpr @ bulkh; (* Important! *)

    rescalep[1] := 1 + bulkh/2; (* Important: we divide by 2! *)
    rescalep[2] := 1 - bulkh/2;
    rescalep[_] := MyError["oops"];
    rescalem[1] := 1 - bulkh/2;
    rescalem[2] := 1 + bulkh/2;
    rescalem[_] := MyError["oops"];

    df[a_, m_]      := setpr[ rescalep[a] df0[a,m] ]; (* \gamma^2 *)
    dfminus[a_, m_] := setpr[ rescalem[a] dfminus0[a,m] ];
    de[a_, m_]      := setpr[ rescalep[a] de0[a,m] ]; (* \Epsilon *)
    deminus[a_, m_] := setpr[ rescalem[a] deminus0[a,m] ];

    (* This remains the same, because the band only shifts, while its total
       spectral weight remains the same. *)
    thetaCh[a_] = 2; (* Denoted as {\bar \gamma} in Campo,Oliveira paper. *)
  ]; (* BAND == "flat" *)

];

timestart["xi"];

hookfile["hook_pre_lanczosinit"];
If[TRI != "manual" && TRI != "manual_nambu" && TRI != "manual_nambu_new" && !option["GENERATE_TEMPLATE"],
  lanczosinit[];
];
hookfile["hook_post_lanczosinit"];

dothelanczos[];

(* Handle exceptions *)
If[TRI == "manual_nambu",
  If[SYMTYPE == "SPU1",
    sckappa[a_][n_] := 0;
    scdelta[a_][n_] := 1(zetaR[a][n]); (* on-site 'rung' term is the same as pairing *)
    (* zetaR is loaded only if RUNGS=True. Beware: not tested, implementation not complete. *)
  ];
  If[SYMTYPE == "SPSU2",
    sckappa[a_][n_] := 0; (* zero for problems without energy-dependence of phases *)
    loadtablescdelta[]; (* loaded from scdelta.dat *)
  ];
  If[SYMTYPE == "P" || SYMTYPE == "PP",
    sckappa[a_][n_] := 0;
    loadtablescdelta[];
  ];
];

If[TRI == "manual_nambu_new",
  If[SYMTYPE == "SPSU2",
    loadtablescdelta[]; (* loaded from scdelta.dat *)
    loadtablesckappa[]; (* loaded from sckappa.dat *)
  ];
];

If[paramexists["disccheck"],
  discretizationChecks[];
];

(* TRAP: At this point, thetaCh should be known. *)
If[!ValueQ[ thetaCh[1] ],
  MyError["Unknown BAND type."];
];

MyPrintForm["BAND=`` thetaCh=``", BAND,
  cstr10 /@ Array[thetaCh, COEFCHANNELS] ];

(* Staggered potential in the Wilson chain for problems with hard gap or
with BCS gap (after suitable p-h and Bogoljubov transformations). NOTE: we
add it to dzeta[] determined by the energy-dependent hybridisation function!
See Yoshioka, Ohashi, J. Phys. Soc. Japan, 69 1812 (2002), Eq. (2.19). Note
that zeta[a][0] = -Delta. This is multiplied by faktor to obtain
-Deltatilde.  See Eq. (2.19) with N=-1 and factor out Lambda^(-1/2) to
obtain H=Himp + Sqrt[Gamma] hop - Delta (n-1). *)

zeta[a_][n_] := dzeta[a][n] +
  If[paramexists["gap"], - (-1)^n paramnum["gap"], 0] +
  If[paramexists["shift0"] && n == 0, paramnum["shift0"], 0];

gapdefined = False;

hookfile["hook_bcs"];

If[!gapdefined && paramexists["bcsgap"],
  gapdefined = True;
  scdelta[_][_] = paramnum["bcsgap"];
  sckappa[_][_] = 0;
];

If[!gapdefined && paramexists["bcsgap1"] && paramexists["bcsgap2"] && !paramexists["bcsgap3"],
  gapdefined = True;
  scdelta[1][_] = paramnum["bcsgap1"];
  scdelta[2][_] = paramnum["bcsgap2"];
  sckappa[_][_] = 0;
];

If[!gapdefined && paramexists["bcsgap1"] && paramexists["bcsgap2"] && paramexists["bcsgap3"],
  gapdefined = True;
  scdelta[1][_] = paramnum["bcsgap1"];
  scdelta[2][_] = paramnum["bcsgap2"];
  scdelta[3][_] = paramnum["bcsgap3"];
  sckappa[_][_] = 0;
];

(* *** Bulk magnetic field support. *** *)

If[POLARIZED && isQSZ[] && paramexists["globalB"],
  For[i = 1, i <= CHANNELS, i++,
    zeta[i][0] = zeta[i][0] + (1/2) paramnum["globalB"];
    zeta[i+CHANNELS][0] = zeta[i+CHANNELS][0] - (1/2) paramnum["globalB"];
  ];
];

If[POLARIZED && isU1[] && paramexists["globalB"],
  For[i = 1, i <= CHANNELS, i++,
    zeta[i][0] = zeta[i][0] - (1/2) paramnum["globalB"];
    zeta[i+CHANNELS][0] = zeta[i+CHANNELS][0] + (1/2) paramnum["globalB"];
  ];
];

(* Compare with globalB! This one modifies all coefficients. *)
If[POLARIZED && (isSPU1[] || isP[] || isPP[] || isNONE[]) && paramexists["globalh"],
  For[i = 1, i <= CHANNELS, i++,
    For[m = 0, m <= mMAX, m++,
      zeta[i][m] = zeta[i][m] + (1/2) paramnum["globalh"];
      zeta[i+CHANNELS][m] = zeta[i+CHANNELS][m] - (1/2) paramnum["globalh"];
    ];
  ];
];


showtable[name_, channel_, table_] := Module[{},
  MyPrintForm["`` (channel ``)", name, channel];
  MyPrint[cstr10[First[#]]]& /@ table;
  (* XXX
  If[(And @@ Map[Element[#, Reals]&, table]) != True,
    MyError["showtable: Coefficients must be real numbers."];
  ];
  *)
];

OUTPREC = 20;

(* Legacy interface for building Wilson chain coefficient tables *)
Module[{a, precxi, preczeta},
  For[a = 1, a <= COEFCHANNELS, a++,
      MyPrintForm["Discretization (channel ``)", a];

      xitable[a]   = Table[{N[bandrescale xi[a][i],   OUTPREC]}, {i, 0, DISCNMAX}];
      zetatable[a] = Table[{N[bandrescale zeta[a][i], OUTPREC]}, {i, 0, DISCNMAX}];

      showtable["xitable", a, xitable[a]];
      showtable["zetatable", a, zetatable[a]];
      MyPrint["Precision last xi:", precxi = Precision[xi[a][DISCNMAX]] ];
      MyPrint["Precision last zeta: ", preczeta = Precision[zeta[a][DISCNMAX]] ];

      If[RUNGS,
        xiRtable[a]   = Table[{N[bandrescale xiR[a][i],   OUTPREC]}, {i, 0, DISCNMAX}];
        zetaRtable[a] = Table[{N[bandrescale zetaR[a][i], OUTPREC]}, {i, 0, DISCNMAX}];

        showtable["xiRtable", a, xiRtable[a]];
        showtable["zetaRtable", a, zetaRtable[a]];
      ];

      If[isSC[],  (* Superconducting host *)
        scdeltatable[a] = Table[{N[bandrescale scdelta[a][i], OUTPREC]}, {i, 0, DISCNMAX}];
        sckappatable[a] = Table[{N[bandrescale sckappa[a][i], OUTPREC]}, {i, 0, DISCNMAX}];

        showtable["scdeltatable", a, scdeltatable[a]];
        showtable["sckappatable", a, sckappatable[a]];
      ];

      eptable[a]  = Table[{N[de[a, m], OUTPREC]},      {m, 0, mMAX}];
      emtable[a]  = Table[{N[deminus[a, m], OUTPREC]}, {m, 0, mMAX}];
      u0ptable[a] = Table[{N[du[a][0, m], OUTPREC]},   {m, 0, mMAX}];
      u0mtable[a] = Table[{N[dv[a][0, m], OUTPREC]},   {m, 0, mMAX}];
    ];
];

MyPrint["Discretization done."];
timeadd["xi"];
