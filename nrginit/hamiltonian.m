(***************************** HAMILTONIAN *******************************)

(* Define default Anderson-like impurity Hamiltonians and provide default
values for nnop[]. *)

adddots[nrdots_] := Module[{},
  H1 = 0; (* Zero-impurity cases and Kondo-like models with spin
             operators only. *)

  If[nrdots >= 1,
    H1 = delta number[d[]] + U/2 pow[number[d[]]-1, 2];
  ];
  If[nrdots >= 2,
    Ha = delta number[a[]] + U/2 pow[number[a[]]-1, 2];
    nnop[a[]] = -1;
  ];
  If[nrdots >= 3,
    Hb = delta number[b[]] + U/2 pow[number[b[]]-1, 2];
    nnop[b[]] = -1;
  ];
  If[nrdots >= 4,
    He = delta number[e[]] + U/2 pow[number[e[]]-1, 2];
    nnop[e[]] = -1;
  ];
  If[nrdots >= 5,
    Hg = delta number[g[]] + U/2 pow[number[g[]]-1, 2];
    nnop[g[]] = -1;
  ];
];

(* Average number of electrons per Wilson chain site. Depends
on the degeneracy, i.e. on the electron spin. *)
AVGOCCUP = (2 BANDSPIN + 1)/2;

(* c is the channel index, i is the chain site index *)
fop[c_Integer, 0, x___]  := f[c, x];
fop[c_Integer, i_Integer, x___] := f[c, i, x];
fopCR[c_Integer, 0, x___]  := f[CR, c, x];
fopCR[c_Integer, i_Integer, x___] := f[CR, c, i, x];
fopAN[c_Integer, 0, x___]  := f[AN, c, x];
fopAN[c_Integer, i_Integer, x___] := f[AN, c, i, x];

(* On-site Hamiltonian on i-th site (0,1,2) of the Wilson chain
for channel ch (1,2,3). *)
HBANDonsite[ch_, i_] := Module[{reg, ireg},
  If[!POLARIZED && !POL2x2,
    reg = coefzeta[ch, i] (number[fop[ch-1, i]] - AVGOCCUP);
  ];
  If[POLARIZED && !POL2x2,
    reg = coefzeta[ch, i]          (number[fop[ch-1, i], UP] - AVGOCCUP/2) +
          coefzeta[ch+CHANNELS, i] (number[fop[ch-1, i], DO] - AVGOCCUP/2);
  ];
  (* Added 10.9.2012 *)
  If[POL2x2,
    reg = coefzeta[ch, i]            (number[fop[ch-1, i], UP] - AVGOCCUP/2) +
          coefzeta[ch+CHANNELS, i]   (number[fop[ch-1, i], DO] - AVGOCCUP/2) +
          coefzeta[ch+2*CHANNELS, i] (fopCR[ch-1, i, UP] ~ nc ~ fopAN[ch-1, i, DO]) +   (* fixed 21.10.2019 *)
          coefzeta[ch+3*CHANNELS, i] (fopCR[ch-1, i, DO] ~ nc ~ fopAN[ch-1, i, UP]); (* Hermitian! *)
  ];

  (* Superconducting pairing contribution. This is isospinx[fop, n=0]. *)
  ireg = coefdelta[ch, i]            nc[fopCR[ch-1, i, UP], fopCR[ch-1, i, DO]] +
         Conjugate[coefdelta[ch, i]] nc[fopAN[ch-1, i, DO], fopAN[ch-1, i, UP]];

  reg + ireg
];

(* Hop with spin change. Two such terms are required to form a Hermitian
   Hamiltonian. *)
hop[op1_?fermionQ[j1___], op2_?fermionQ[j2___], sigma1_, sigma2_] :=
  op1[CR, j1, sigma1] ~ nc ~ op2[AN, j2, sigma2] +
  op2[CR, j2, sigma2] ~ nc ~ op1[AN, j1, sigma1];

(* Hopping Hamiltonian term for hopping from i-1-th to i-th site. Note that
the order of arguments to anomaloushop[] is important to get the sign right. *)

HBANDhop[ch_, 0]  := 0;

HBANDhop[ch_, i_] := Module[{reg, ireg},
  If[!POLARIZED && !POL2x2,
    reg = coefxi[ch, i-1] hop[fop[ch-1, i-1], fop[ch-1, i]];
  ];
  If[POLARIZED && !POL2x2,
    reg = coefxi[ch, i-1]          hop[fop[ch-1, i-1], fop[ch-1, i], UP] +
          coefxi[ch+CHANNELS, i-1] hop[fop[ch-1, i-1], fop[ch-1, i], DO];
  ];
  (* Added 10.9.2012 *)
  If[POL2x2,
    reg = coefxi[ch, i-1]          hop[fop[ch-1, i-1], fop[ch-1, i], UP] +
          coefxi[ch+CHANNELS, i-1] hop[fop[ch-1, i-1], fop[ch-1, i], DO] +
          coefxi[ch+2*CHANNELS, i-1] hop[fop[ch-1, i-1], fop[ch-1, i], UP, DO] +
          coefxi[ch+3*CHANNELS, i-1] hop[fop[ch-1, i-1], fop[ch-1, i], DO, UP];
  ];

  (* Support for anomalous hopping terms in the presence of superconductivity:
     these are the f^\dag_{i-1} f^\dag_i terms. *)
  ireg = coefkappa[ch, i-1] anhop[fop[ch-1, i-1], fop[ch-1, i]];
  reg + ireg
];

(* NOTE: H0 is the chain Hamiltonian for problems where the band is not
particle-hole symmetric, or for problems where we take into account more
than the size 0 of the Wilson chain (see parameter "Ninit"). *)

HBAND[] := Module[{H0onsite, H0hop, H0bcs},
  H0onsite = Sum[HBANDonsite[ch, i], {ch, CHANNELS}, {i, 0, Ninit}];
  H0hop =    Sum[HBANDhop[ch, i],    {ch, CHANNELS}, {i, 0, Ninit}];
  H0 = H0onsite + H0hop;

  (* If we don't work with SPSU2 symmetry type, we drop all anomalous terms. *)
  If[!isSC[],
    H0 = H0 /. { coefdelta[___] -> 0, coefkappa[___] -> 0};
  ];

  If[RUNGS && CHANNELS == 2,
    H0r = Sum[coefrung[1,i] hop[fop[0, i], fop[1, i]], {i, 0, Ninit}];
    MyPrint["H0r=", H0r];
    H0 = H0 + H0r;
  ];

  MyPrint["H0=", H0];
];

(* The default band-impurity coupling Hamiltonian. *)
HC[] := Module[{},
  (* Note: gammaPolCh/zeta/... are indexed as ch=1, 2, ..., while f[]
  operators are indexed as ch=0, 1, ... *)

  If[!POLARIZED && !POL2x2,
    Hc = Sum[gammaPolCh[ch] hop[f[ch-1], d[]], {ch, CHANNELS}];
  ];
  If[POLARIZED && !POL2x2,
    Hc = Sum[gammaPolCh[ch]          hop[f[ch-1], d[], UP], {ch, CHANNELS}] +
         Sum[gammaPolCh[ch+CHANNELS] hop[f[ch-1], d[], DO], {ch, CHANNELS}];
  ];
  (* Added 10.9.2012 *)
  If[POL2x2,
    Hc = Sum[gammaPolCh[ch]            hop[f[ch-1], d[], UP],     {ch, CHANNELS}] +
         Sum[gammaPolCh[ch+CHANNELS]   hop[f[ch-1], d[], DO],     {ch, CHANNELS}] +
         Sum[gammaPolCh[ch+2*CHANNELS] hop[f[ch-1], d[], UP, DO], {ch, CHANNELS}] +
         Sum[gammaPolCh[ch+3*CHANNELS] hop[f[ch-1], d[], DO, UP], {ch, CHANNELS}];
  ];
];

(* Called from def1ch[] and def2ch[] to set up the support for spin-polarized
conduction-band calculations. *)
InitPolarized[] := Module[{},
  (* Defaults. *)
  COEFCHANNELS = CHANNELS;
  VDIM = CHANNELS; (* NEW 2019: dimension of hybridisation matrix. Used for BAND=manual_V *)

  (* Changes in the case of spin-polarized conduction bands. *)
  If[POLARIZED,
    If[!MemberQ[{"QSZ", "U1", "SPU1", "P", "PP", "NONE"}, SYMTYPE],
      MyError["Spin polarized calculation not supported with chosen SYMTYPE."]
    ];
    COEFCHANNELS = 2 CHANNELS; (* Double the number of coefficient sets for channels. *)
    VDIM = 2 * CHANNELS;
  ];
  If[POL2x2,
    If[!MemberQ[{"U1"}, SYMTYPE],
      MyError["Spin 2x2 structure not supported with chosen SYMTYPE."];
    ];
    COEFCHANNELS = 4 CHANNELS; (* Quadruple the number of coefficient sets for channels. *)
    VDIM = 2 * CHANNELS; (* !! *)
  ];
  If[POLARIZE && POL2x2,
    MyError["POLARIZED or POL2x2? Choose one!"];
  ];
];

(* Default settings for 1-channel models *)
def1ch[nrdots_:1] := Module[{},
  CHANNELS = 1;
  NRDOTS = nrdots;
  MyPrint["def1ch, NRDOTS=", NRDOTS];

  InitPolarized[];

  HC[];
  HBAND[];

  (* Default operator numbering for isospin symmetry generation: d[]
  has an inverted sign, but the first site in the Wilson's chain does
  not! *)
  nnop[d[] ] = -1;
  nnop[f[0]] = 0;
  nnop[f[0, 1]] = -1;

  adddots[NRDOTS];
];

(* Default settings for 2-channel models *)
def2ch[nrdots_:1] := Module[{},
  CHANNELS = 2;
  NRDOTS = nrdots;
  MyPrint["def2ch, NRDOTS=", NRDOTS];

  InitPolarized[];

  HC[];
  HBAND[];

  nnop[f[0]] = 0;
  nnop[f[1]] = 0;

  (* For SYMTYPE=ISO, we must have nnop[f[0]] = nnop[f[1]] = 0. *)
  (* For SYMTYPE=ISO2, we must have nnop[f[0]] = 0, nnop[f[1]] = 1. *)
  If[SYMTYPE == "ISO2" || SYMTYPE == "ISO2LR",
    nnop[f[0]] = 0;
    nnop[f[1]] = 1;
  ];

  (* WARNING: in TQD, d is the impurity in the middle, so the
  default rule for d[] has to be overriden!! *)

  nnop[d[] ] = -2;

  adddots[nrdots];
];

(* Default settings for 3-channel models *)
def3ch[nrdots_:1] := Module[{},
  CHANNELS = 3;
  NRDOTS = nrdots;
  MyPrint["def3ch, NRDOTS=", NRDOTS];

  InitPolarized[];

  HC[];
  HBAND[];

  (* The current implementation of SYMTYPE=ISO is such that nnop[f[i]]
     are all equal! *)
  nnop[f[0]] = 0;
  nnop[f[1]] = 0;
  nnop[f[2]] = 0;

  nnop[d[] ] = -1;

  adddots[nrdots];
];

(* Default settings for 4-channel models *)
def4ch[nrdots_:1] := Module[{},
  CHANNELS = 4;
  NRDOTS = nrdots;
  MyPrint["def4ch, NRDOTS=", NRDOTS];

  InitPolarized[];

  HC[];
  HBAND[];

  (* The current implementation of SYMTYPE=ISO is such that nnop[f[i]]
     are all equal! *)
  nnop[f[0]] = 0;
  nnop[f[1]] = 0;
  nnop[f[2]] = 0;
  nnop[f[3]] = 0;

  nnop[d[] ] = -1;

  adddots[nrdots];
];

CHANNELS = -1; (* Bug trap *)
NRDOTS = -1; (* Bug trap *)
MAKESPINKET = Null; (* Operator(s?) that is converted to spin kets *)
MAKEORBKET = Null;
MAKEPHONON = Null; (* 1 = one phonon mode, etc. *)
BZSPIN = Null;

(** Reflection symmetry (parity) **)

(* List 'lrchain' holds the operators that describe the ordering of sites in
the real-space configuration of the impurities, for example: {f[0], d[],
f[1]} for two-channel single-impurity Kondo problem. This list is used to
generate parity-adapted basis states. List 'lrextrarule' contains additional
transformation rules that need to be applied to the basis when constructing
the mirror-symmetric states. For example: in the case of antisymmetric
coupling of a phonon mode to a hopping term, all phonon kets need to be
multiplied by -1. *)

lrchain = {};
lrextrarule = {};


(* Some useful functions for Hamiltonian construction and error checking. *)
checkQSZ[] := If[Nor[isQSZ[], isP[], isPP[], isNONE[], isDBLQSZ[], isDBLSU2[], isDBLISOSZ[], isSU2[], isU1[],
                     isSPU1[], isSPU1LR, isISOSZ[]],
 MyError["SYMTYPE==QSZ (and similar) only"]
];

(**************************** MODEL DEPENDENT ******************************)

(* Trivial examples: no impurity and single impurity Anderson model *)

(* No impurity: clean conduction band. WARNING: thermodynamic properties
   of the band still depend on the "hybridization function"! *)
If[ MODEL == "CLEAN",
  def1ch[0];
  H = H0;
];

(* Single Impurity Anderson Model and its variants *)
If[ MODEL == "SIAM" && (VARIANT == "" || VARIANT == "MAGFIELD"),
  def1ch[1];

  (* SIAM in magnetic field. Only for SYMTYPE=QSZ. *)
  If[VARIANT == "MAGFIELD",
    checkQSZ[];
    H1 = H1 + B spinz[d[]]; (* Zeeman term *)
  ];

  H = H0 + H1 + Hc;
];

(* Load additional model definitions *)
loadmodule["models.m"];
loadmodule["custommodels.m", False];

(* If a model name ends in .m, we load a package file of this name!
This is an alternative to adding new model definitions to custommodels.m *)

SCRIPTMODEL = False;
If[StringLength[MODEL] >= 2 && StringTake[MODEL, -2] == ".m",
  loadmodule[MODEL];
  SCRIPTMODEL = True; (* Use to generate suitable filename *)
];
