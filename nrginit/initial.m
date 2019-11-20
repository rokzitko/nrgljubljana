
(* 
   NRG Ljubljana -- initial.m -- Basis construction, initial Hamiltonian
   diagonalisation and calculation of irreducible matrix elements 
  
   Copyright (C) 2005-2019 Rok Zitko

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.
   
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
   
   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

   Contact information:
   Rok Zitko
   F1 - Theoretical physics
   "Jozef Stefan" Institute
   Jamova 39
   SI-1000 Ljubljana
   Slovenia

   rok.zitko@ijs.si
*)

VERSION = "2019.11";

(* Logging of Mathematica output: this is useful for bug hunting *)
If[!ValueQ[mmalog],
  mmalog = OpenWrite["mmalog"];
  AppendTo[$Output, mmalog];
  AppendTo[$Messages, mmalog]; (* Important! *)
];

SetOptions[$Output, PageWidth -> 240];
SetOptions[$Messages, PageWidth -> 240];

(* Generalization of Print[] which wraps all strings in StandardForm[]
format directives to remove the quotation marks. *) 
Print2[l__] := Print @@ ({l} /. x_String -> StandardForm[x]);

Print2["NRG Ljubljana ", VERSION, " (c) Rok Zitko, rok.zitko@ijs.si, 2005-2019"];

(* Print a warning/error message at specified verbosity level *)
DEBUG = 1;
MyPrint[msg__] := If[DEBUG >= 1, Print2[msg]];
MyPrintForm[form_, msg___] := If[DEBUG >= 1, 
  Print2 @ ToString @ StringForm[form, Sequence @@ (InputForm /@ {msg})]
];
MyVPrint[verbosity_, msg__] := If[DEBUG >= verbosity, Print2[msg]];

(* Prints an error message and exits *)
If[!ValueQ[ExitOnError], ExitOnError = True];
MyError[msg__] := (Print2[msg]; Print2["Aborting.\n"]; 
  If[ExitOnError, Exit[1] ]);

(* Prints a warning message -- does not exit *)
MyWarning[msg__] := Print2[msg];

SetAttributes[MyAssert, HoldFirst];
MyAssert[expr_] := If[Evaluate[expr] =!= True, MyError["Assertion failed: ", Hold[expr] ]];

Off[InterpolatingFunction::dmvali];
Off[InterpolatingFunction::dmval];

(* Convert to a C-language compatible string, truncate to 10 significant
digits *)
cstr10[x_] := ToString[N[x, 10], CForm];

c10[x_] := CForm[N[x,10]]; (* used in dmft.m *)

(* Check version of sneg *)
majorversion[x_] := ToExpression @ First @ StringSplit[x, "."];
minorversion[x_] := ToExpression @ Last @ StringSplit[x, "."];

VERrequired = "1.188";
VERsneg = $SnegVersion;
If[!( (majorversion[VERsneg] > majorversion[VERrequired]) ||
     ((majorversion[VERsneg] == majorversion[VERrequired]) &&
      (minorversion[VERsneg] >= minorversion[VERrequired])) ),
   MyError["Lanski sneg: ", VERsneg,
           " Required version of SNEG is ", VERrequired];
];

MyPrint["Mathematica version: ", $Version];
MyPrint["sneg version: ", $SnegVersion];

(* No recursion limit *)
$RecursionLimit = Infinity;

(* Timing code. Used for profiling and benchmarking. *)
ClearAll[timingdata, time0];
timingdata[_] = 0;
timestart[name_] := (time0[name] = AbsoluteTime[]);
timeadd[name_] := Module[{time1, timediff},
  time1 = AbsoluteTime[];
  timediff = time1-time0[name];
  timingdata[name] += timediff;
];
timereport[] := Module[{t},
  MyPrint["Timing report"];
  t = Map[{#[[1,1,1]], #[[2]]}&, DownValues[timingdata]];
  Scan[MyPrint, t];
];

(* Load a module, Get::noopen message is suppressed. *)
silentGet[x__] := Module[{ret},
  Off[Get::noopen];
  ret = Get[x];
  On[Get::noopen];
  ret
];

(* Load an external module (another .m) file. *)
loadmodule[filename_String, exitonfailure_:True] := Module[{ret},
  MyPrint["Loading module ", filename];
  ret = silentGet[filename, Path -> modulespath];
  If[ret === $Failed,
    MyWarning["Can't load " <> filename <> ". " <>
      If[exitonfailure, "Aborting.", "Continuing."]];
    If[exitonfailure == True, Exit[1]];
  ];
  ret (* The result is returned! *)
];

(* External hook: if the variable 'var' ends with a ".m" suffix, 
   the corresponding external module is called. *)
hook[var_] := If[StringTake[var,-2] == ".m", loadmodule[var, True]];

(* List of directories in which we look for modules. *)
modulespath = {".", ".."};
If[ValueQ[PACKAGEPATH], modulespath = Join[modulespath, PACKAGEPATH]];

(* List of directories in which we look for dump files of basis vectors and
Hamiltonian matrices. In addition to the current directory, we also look
in the parent; this is useful for parameter sweeps. *)
dumppath = {".", ".."};

(* Parse parameter file, if not already performed by nrginit *)
loadmodule["initialparse.m"];
If[PARSED =!= True,
  MyPrint["Parsing parameters"];
  parse["param"];
];

(* External hook with a named external module *)
hookfile[keyword_] := Module[{filename},
  filename = paramdefault[keyword, ""];
  If[filename =!= "",
    loadmodule[filename, True]
  ];
];

(***************************************************************)

If[!ValueQ[SYMTYPE],  SYMTYPE = "QS"];

(* Implemented symmetry types are:
 "QS", ==>                 U(1)_charge x SU(2)_spin
 "QST", ==>                U(1)_charge x SU(2)_spin x SO(3)_orbital [three orbitals]
 "QSTZ", ==>               U(1)_charge x SU(2)_spin x U(1)_orbital [three orbitals]
 "QSZTZ", ==>              U(1)_charge x U(1)_spin x U(1)_orbital [three orbitals]
 "QJ", ==>                 U(1)_charge x SU(2)_total_momentum [three orbitals]
 "ISO" and "ISO2", ==>     SU(2)_isospin x SU(2)_spin
 "QSZ", ==>                U(1)_charge x U(1)_spin
 "ISOSZ", ==>              SU(2)_isospin x U(1)_spin
 "ISOSZLR", ==>            SU(2)_isospin x U(1)_spin x Z_2
 "ISOLR" and "ISO2LR", ==> SU(2)_isospin x SU(2)_spin x Z_2
 "QSLR", ==>               U(1)_charge x SU(2)_spin x Z_2
 "QSC3", ==>               U(1)_charge x SU(2)_spin x Z_3 [3 channels]
 "QSZLR", ==>              U(1)_charge x U(1)_spin x Z_2
 "SU2", ==>                SU(2)_charge
 "DBLSU2", ==>             SU(2)_charge1 x SU(2)_charge2
 "DBLISOSZ", ==>           SU(2)_charge1 x SU(2)_charge2 x U(1)_spin
 "U1", ==>                 U(1)_charge
 "SPSU2", ==>              SU(2)_spin
 "SPU1", ==>               U(1)_spin
 "SPU1LR", ==>             U(1)_spin x Z_2
 "SPSU2LR", ==>            SU(2)_spin x Z_2
 "SPSU2C3", ==>            SU(2)_spin x Z_3 [3 channels]
 "SPSU2T", ==>             SU(2)_spin x SO(3)_orbital [three orbitals]
 "SL", ==>                 spinless fermions [U(1)_charge symmetry]
 "SL3", ==>                spinless fermions [threefold U(1)_charge]
 "ANYJ", ==>               U(1) x U(1), arbitrary spin of the conduction-band electrons
 "P", ==>                  Z_2 fermion parity,
 "PP", ==>                 (Z_2)^2 fermion parity,
 "NONE", ==>               no symmetry
*)

If[SYMTYPE == "runtime",
  If[!paramexists["symtype"],
    MyError["Cannot determine the symmetry type."];
  ];
  SYMTYPE = param["symtype"];
];

knownsymtypes = 
{"QS", "QST", "QSTZ", "QSZTZ", "QJ", "ISO", "ISO2", "QSZ", "ISOLR", "ISO2LR", "QSLR", "QSC3",
"QSZLR", "DBLSU2", "DBLISOSZ",
"SU2", "U1", "SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2C3", "SPSU2T",
"SL", "SL3", "ANYJ", "P", "PP", "NONE", "ISOSZ", "ISOSZLR"};
If[!(MemberQ[knownsymtypes, SYMTYPE]),
  MyError["Unknown SYMTYPE."];
];

(* isLR[] returns True if the problem has Z_2 parity symmetry. For symmetry 
types in this list, we create a parity-adapted set of basis states as 
step number 5 in the basis-state-generation part of the code. *)
lrsymtypes = {"QSLR", "QSZLR", "ISOLR", "ISO2LR", "ISOSZLR", "SPU1LR", "SPSU2LR"};
isLR[] := MemberQ[lrsymtypes, SYMTYPE]; 

qstypes = {"QS", "QSLR"};
qsztypes = {"QSZ", "QSZLR"};
isotypes = {"ISO", "ISO2", "ISOLR", "ISO2LR"};
isosztypes = { "ISOSZ", "ISOSZLR" };
sctypes = {"SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2T", "P", "PP", "NONE", "SPSU2C3"};
orbtypes = {"QST", "SPSU2T", "QSTZ", "QSZTZ"}; (* quantum number T *)
isQS[] := MemberQ[qstypes, SYMTYPE];
isQSZ[] := MemberQ[qsztypes, SYMTYPE];
isISO[] := MemberQ[isotypes, SYMTYPE];
isISOSZ[] := MemberQ[isosztypes, SYMTYPE];
isSC[] := MemberQ[sctypes, SYMTYPE];
isORB[] := MemberQ[orbtypes, SYMTYPE];

isSU2[] := (SYMTYPE === "SU2");
isDBLSU2[] := (SYMTYPE === "DBLSU2");
isDBLISOSZ[] := (SYMTYPE === "DBLISOSZ");
isU1[] := (SYMTYPE === "U1");
isSPSU2[] := (SYMTYPE === "SPSU2");
isSPU1[] := (SYMTYPE === "SPU1");
isSPU1LR[] := (SYMTYPE === "SPU1LR");
isSPSU2LR[] := (SYMTYPE === "SPSU2LR");
isSL[] := (SYMTYPE === "SL");
isSL3[] := (SYMTYPE === "SL3");
isANYJ[] := (SYMTYPE === "ANYJ");
isP[] := (SYMTYPE === "P");
isPP[] := (SYMTYPE === "PP");
isNONE[] := (SYMTYPE === "NONE");
isQST[] := (SYMTYPE === "QST");
isQSTZ[] := (SYMTYPE === "QSTZ");
isQSZTZ[] := (SYMTYPE === "QSZTZ");
isQJ[] := (SYMTYPE === "QJ");
isSPSU2T[] := (SYMTYPE === "SPSU2T");
isQSC3[] := (SYMTYPE === "QSC3");
isSPSU2C3[] := (SYMTYPE === "SPSU2C3");
(* TODO: remove the above! *)
is[sym_] := SYMTYPE === sym;

(********)

MODEL          = paramdefault["model", "SIAM"];
VARIANT        = paramdefault["variant", ""];
OPTIONS        = paramdefault["options", ""];
PERTURB        = paramdefault["perturb", ""];
OPS            = paramdefault["ops", ""];
lambda         = N @ paramdefaultnum["Lambda", 2.0];
z              = N @ paramdefaultnum["z", 1.0];
BAND           = paramdefault["band", "flat"];
DEBUG          = paramdefaultnum["mmadebug", 1];
bandrescale    = paramdefaultnum["bandrescale", 1];
DISCRETIZATION = paramdefault["discretization", "Z"];

If[StringLength[DISCRETIZATION] >= 1,
  DISCRETIZATION = StringTake[DISCRETIZATION, {1}] ];

knowndiscretizationtypes = {"C", "Y", "Z"};
If[!(MemberQ[knowndiscretizationtypes, DISCRETIZATION]),
  MyError["Unknown DISCRETIZATION."];
];

DY = DISCRETIZATION == "Y"; (* Original logarithmic mesh approach *)
DC = DISCRETIZATION == "C"; (* Campo and Oliveira's approach *)
DZ = DISCRETIZATION == "Z"; (* My little modification of "C" *)

(* Are the discretization coefficients for spin up and spin down different? *)
(* For QSZ and U1, this effectively doubles the number of coefficient sets,
   one number of sets for spin-up electrons, the other for spin-down electrons.
   For SPU1, the number is also doubled, but the off-diagonal elements are
   identical in both sets. *)
POLARIZED = paramdefaultbool["polarized", False];

(* For U1, this effectively quadruples the number of coefficient sets. This
   allows to describe the full 2x2 structure in the spin space. *)
POL2x2 = paramdefaultbool["pol2x2", False];

(* Allow channel-mixing terms in the Wilson chain. *)   
RUNGS = paramdefaultbool["rungs", False];

(* It is possible to take into account more than a single site of the Wilson
chain in the initial.m part of the code. This is controlled by parameter
"Ninit"; its value can be interpreted as the highest index of the f-orbital
that is still retained in the initial Hamiltonian. *)

Ninit = paramdefaultnum["Ninit", 0];

(* Commonly used model parameters. These four are special, 
since they can be defined in either [extra] or [param] blocks,
and they have default values if not defined anywhere. *)

getmodelparam[key_, default_] := Module[{},
  If[paramexists[key, "extra"], Return[paramnum[key, "extra"]]];
  If[paramexists[key, "param"], Return[paramnum[key, "param"]]];
  Return[default];
];

realU     = getmodelparam["U", 0.1];
realGamma = getmodelparam["Gamma", 0.1];
realdelta = getmodelparam["delta", 0.];
realt     = getmodelparam["t", 0.];

parsevalue[str_String] /;  klicaj[str] := ToExpression[StringDrop[str, 1]];
parsevalue[str_String] /; !klicaj[str] := importnum[str];        

(* All PARAM=VALUE lines in [extra] block of the input file get transformed 
   into PARAM->VALUE rules in params. Since rules are applied in the order 
   of their appearance in the list, this implies that the preexisting rules
   take precendence over these automatically appended ones. *)
addextraparams[] := Module[{params2},
  If[ValueQ[listdata["extra"]],
    params2 = Map[ToExpression[First[#]] -> parsevalue @ Last[#] &, 
                  listdata["extra"]];
    params = Join[params, params2];
  ];                 
];

(* Define extraPARAM=VALUE for all PARAM=VALUE lines in [extra]
   block of the input file (for backward compatibility). *)
If[ValueQ[listdata["extra"]],
   exmap = Map[ {ToExpression["extra" <> First[#]],     
          parsevalue @ Last[#] }&, listdata["extra"] ];
   MapThread[Set, Transpose[exmap]];
];

If[!ValueQ[BASISRULE],BASISRULE = ""]; (* Transformation rule for basis states *)
(* This can be used, for example, to project out doubly occupied impurity
state to simulate U=infinity Anderson model, etc. *)

(* Workaround: StringSplit is only available in Mathematica starting in 5.1.
The following routine thus extends the compatibility of "NRG Ljubljana" to
Mathematica 5.0. *)
MyStringSplit[l_]:=Module[{c, p},
      c=Characters[l];
      p=Position[c," "]//Flatten;
      p=Join[{0},p,{Length[c]+1}];
      Table[StringTake[l,{p[[i-1]]+1,p[[i]]-1}],{i,2,Length[p]}]
];
MyStringSplit[""] := {};      

(* OPS is a space delimited list of operators that we request to compute
during the NRG iteration. *)
(* 7.3.2016: sort alphabetically *)
lops = Sort @ MyStringSplit[OPS];

(* calcopq[op] returns True if we requested calculation of operator 'op'. 
This is used in the generation of the input file for NRG iteration. *)
calcopq[op_] := MemberQ[lops, op];

(* Generate a list of all operators with a given prefix. *)
(* Returns a list of {string, suffix, expr}, where 'string' is the matching
   string, while 'suffix' is the string without the prefix, while 'expr'
   is 'suffix' with eventual leading ( and trailing ) stripped. *)
calcoplist[prefix_] := Module[{len, l, stringstrip},
  len = StringLength[prefix];
  l = Select[lops, StringLength[#] >= len &];
  l = Select[l, StringTake[#, len] == prefix &];

  (* stringstrip[] is ugly, but works with Mathematica 5.0 *)
  stringstrip[str_] := 
    stringstrip @ StringDrop[str, 1] /; 
    StringLength[str] >= 1 && StringTake[str, 1] == "(";
  stringstrip[str_] := 
    stringstrip @ StringDrop[str, -1] /; 
    StringLength[str] >= 1 && StringTake[str, -1] == ")";
  stringstrip[str_] := str;

  If[l === {},
    {},
    Map[{#, tempstr=StringDrop[#, len], stringstrip @ tempstr}&, l]
  ]
];

(* OPTIONS is a space delimited list of additional options that are
defined in the parameters file. *)
loptions = MyStringSplit[OPTIONS];
MyPrint["Options: ", loptions];

(* option[keyword] returns True if option 'keyword' is specified *)
option[keyword_] := MemberQ[locateoption[keyword], True];

(* Careful: locateoption[] is not fool-proof. This should better be done
with regex searches, but this only appeared in Mathematica starting in
version 5.1 *)

locateoption[keyword_] := Map[StringPosition[#, keyword] != {} &, loptions];

(* Options can have values: key=value *)
optionvalue[keyword_] := Module[{loc, pos, pair},
  loc = locateoption[keyword];
  If[Count[loc, True] != 1,
    MyError["Error parsing keyword:", keyword];
  ];
  pos = Position[loc, True] [[1,1]];
  pair = loptions [[pos]];
  If[StringLength[pair] <= StringLength[keyword],
    MyError["Keyword found, no value:", pair];
  ];
  StringDrop[pair, StringLength[keyword]+1] (* A string is returned! *)
];

(* We allow for dynamic adding of keywords. For some models, specific
workarounds are required which can be enabled at run-time. *)

addoption[keyword_] := AppendTo[loptions, keyword];

(* 
KNOWN OPTIONS 
============= 
PARAMPRE - apply parameters before calling matrixrepresntationvc[]
WRITE - write basis and Hamiltonian matrix files
READBASIS - read basis from a file
READHAM - read the Hamiltonian matrices from a file
EPSCLIP - clip nonrepresentably small (in double type) floating points to zero
TEMPLATE - create a template for the output file 'data' rather than an actual output file.
           This is used to create the 'data' file on computers without Mathematica.
LRTRICK - rewrite the basis states so that they are even viz. odd wrt parity
COMPLEX - enforce the generation of a complex-numbers version of the data file
NOSHUR - do not use the Shur decomposition to diagonalize matrices
MPVCSLOW - set MPVCFAST=False
GENERATE_TEMPLATE - generate a template file with no Wilson coefficients, to be used
     in the NRG Ljubljana <-> TRIQS interface
*)

(* Spin of conduction-band electrons. *)
BANDSPIN = 1/2; (* Default for all codes except ANYJ. *)

If[SYMTYPE == "ANYJ",
  If[!paramexists["spin"],
    MyError["Define the spin of the conduction band electrons."];
  ];
  BANDSPIN = ((ToExpression @ param["spin"])-1)/2;
  MyPrint["BAND SPIN=", BANDSPIN];
];

(***************************** GENERAL STUFF *******************************)

(* 
 The following parameters are common to all models:
 delta = deviation from p-h symmetry, delta=epsilon+U/2
 U = e-e repulsion parameter
 gamma = pi rho |V|^2, hybridisation strength, gamma/D = (t'/t)^2 for 
   single embedded dot (model SIAM). It is not used directly.
 gammaPol = sqrt{gamma} \sim V, i.e. proportional to hopping 
 t = coupling of the side dot
 For other parameters, consult the Hamiltonian definitions!

 NOTE: all dimensionfull parameters are expressed in units of the   
 bandwidth D. 
*)

snegrealconstants[delta, U, gammaPol, t];

SetAttributes[gammaPolCh, NumericFunction];
SetAttributes[hybV, NumericFunction];
SetAttributes[coefzeta,  NumericFunction];
SetAttributes[coefrung, NumericFunction];
SetAttributes[coefxi, NumericFunction];
SetAttributes[coefdelta, NumericFunction];
SetAttributes[coefkappa, NumericFunction];

(* Assumption: Wilson chain coefficients are real. This is not always the case! *)
Conjugate[gammaPolCh[i__]] ^= gammaPolCh[i];
Conjugate[hybV[i__]] ^= hybV[i];
Conjugate[coefzeta[i__]] ^= coefzeta[i];
Conjugate[coefrung[i__]] ^= coefrung[i];
Conjugate[coefxi[i__]] ^= coefxi[i];
Conjugate[coefdelta[i__]] ^= coefdelta[i];
Conjugate[coefkappa[i__]] ^= coefkappa[i];

(* Quantities 'theta0' and 'gammaA' are defined below, when
discretization is set up. *)

params = { 
    gammaPol        -> Sqrt[(1/Pi) theta0       gammaA], (* GP *)
    gammaPolCh[ch_] :> Sqrt[(1/Pi) theta0Ch[ch] gammaA],
    hybV[i_,j_]     :> Sqrt[1/Pi] V[i,j],    
    coefzeta[ch_, j__] :> N[ bandrescale zeta[ch][j] ], (* channel index = 1,2,3 *)
    coefxi[ch_, j__]   :> N[ bandrescale xi[ch][j] ], (* N[] -> MachinePrecision! *)
    coefrung[ch_, j__]  :> N[ bandrescale zetaR[ch][j] ],
    coefdelta[ch_, j__] :> N[ bandrescale scdelta[ch][j] ],
    coefkappa[ch_, j__] :> N[ bandrescale sckappa[ch][j] ], (* YYY *)
    U              -> realU,
    delta          -> realdelta,
    t              -> realt
};

(* See M. Sindel, PhD dissertation, Appendix A.1 *)
(* \sqrt{ (\Gamma/\pi) \theta } *)
(* \theta_0 = \int_{-1}^{1} \Gamma(\epsilon) d\epsilon *)
(* \theta_0 = \theta \Gamma, i.e. \Gamma defines the overall scale, 
   \theta is the integral of the frequency dependant part. *)

(* CONVENTION: f[0], f[1] are zero sites of the Wilson chain for the first
and the second channel, d[] is the impurity orbital, a[], b[], e[], g[] are
additional impurity orbitals. *)

snegfermionoperators[{f, BANDSPIN}, a, b, d, e, g];

(* Declare additional parameters using snegrealconstants[]. addexnames[] is
called from maketable[] or manually when debugging! *) 

addexnames[] := Module[{exnames},
  exnames = Map[StringDrop[#,5]&, Names["extra*"]];
  MyPrint["exnames=", exnames];
  snegrealconstants @@ Map[Symbol, exnames];
];


(***************************** HAMILTONIAN *******************************)

(* Define default Anderson-like impurity Hamiltonians and provide default
values for nnop[]. *)

adddots[nrdots_] := Module[{},
  MyPrint["adddots, nrdots=", nrdots];

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
  ireg = coefdelta[ch, i] (nc[fopCR[ch-1, i, UP], fopCR[ch-1, i, DO]] + 
                           nc[fopAN[ch-1, i, DO], fopAN[ch-1, i, UP]]);

  reg + ireg
];

(* Anomalous hopping operator *)
SetAttributes[anomaloushop, Listable];
anomaloushop[op1_?fermionQ[j1___], op2_?fermionQ[j2___], sigma_] :=
  op1[CR, j1, sigma] ~ nc ~ op2[CR, j2, 1-sigma] +
  op2[AN, j2, 1-sigma] ~ nc ~ op1[AN, j1, sigma];

(* NOTE THE MINUS SIGN! *)
anomaloushop[op1_?fermionQ[j1___], op2_?fermionQ[j2___]] /;
  (spinof[op1] == spinof[op2] == 1/2) :=
  anomaloushop[op1[j1], op2[j2], UP] - anomaloushop[op1[j1], op2[j2], DO];

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

  (* Add support for RUNGS *)

  ireg = coefkappa[ch, i-1] anomaloushop[fop[ch-1, i-1], fop[ch-1, i]];
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
  
  MyPrint["COEFCHANNELS:", COEFCHANNELS];  
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
checkQSZ[] := If[Nor[isQSZ[], isP[], isPP[], isNONE[], isDBLSU2[], isDBLISOSZ[], isSU2[], isU1[], 
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

(* Important: all left sides of the table params should be defined as sneg
parameters. Currently they are all real constants! Call addparamnames[]
after adding parameters to params list. *)
(* TODO: only if symbols? *)

addparamnames[] := Module[{paramnames},   
   paramnames = params[[All, 1]];
   paramnames = Select[paramnames, AtomQ]; (* Drop parametrized rules *)
   MyVPrint[3, "Declaring as constants: ", paramnames];
   snegrealconstants @@ paramnames;
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

addextraparams[];
addparamnames[];
MyPrint["params=", params];

If[NRDOTS == -1, MyError["NRDOTS = -1"]]; (* Bug trap *)
If[CHANNELS == -1, MyError["CHANNELS = -1"]]; (* Bug trap *)

MyPrint["NRDOTS:", NRDOTS];
MyPrint["CHANNELS:", CHANNELS];

basopsch = Table[f[n], {n, 0, CHANNELS-1}];
Do[basopsch = Join[basopsch, Table[f[n, i], {n, 0, CHANNELS-1}]], {i, 1, Ninit}];
basopsdot = Take[{d[], a[], b[], e[], g[]}, NRDOTS];
basisops = basops = Sort[Join[basopsch, basopsdot]];

MyPrint["basis:", basisops];
MyPrint["lrchain:", lrchain];
MyPrint["lrextrarule:", lrextrarule];

(* Check for consistency between the number of operators and the numbers of
   dots and channels *)
NROPS = Length[basisops];
If[NROPS != NRDOTS + CHANNELS (1+Ninit),
  MyError["Number of operators does not match what we were expecting."];
];

MyPrint["NROPS:", NROPS];


If[NRDOTS >= 1,
(** We define some useful operators, such as total spin and total charge. **)

  opstot = Plus @@ Map[spinxyz, basopsdot];
  opq = Plus @@ Map[(number[#]-1)&, basopsdot];
  ntot = Plus @@ Map[number, basopsdot];

  (* Square of total occupancy *)
  ntot2 = pow[ntot, 2];

  (* Total S^2 operator. *)
  ops2 = inner[opstot, opstot] //Expand;

  (* Total S_Z, S_+ and S_- operators. I define them here mostly for
  convenience and for testing purposes. *)
  opsz = opstot[[3]];
  opsplus =  opstot[[1]] + I opstot[[2]];
  opsminus = opstot[[1]] - I opstot[[2]],
(* else *)
  opstot = opq = ntot = ntot2 = ops2 = opsz = opsplus = opsminus = 0;
];

MyVPrint[3, "opstot:", opstot];
MyVPrint[3, "opq:", opq];
MyVPrint[3, "ntot:", ntot];

(* Even combination of (all) creation operators *)
If[NRDOTS >= 1,
  creven = (Plus @@ Map[# /. op_[] -> op[CR, sigma] &, basopsdot]) / Sqrt[NRDOTS];
  aneven = conj[creven];
  neven = Sum[ nc[creven, aneven], {sigma, 0, 1}] // Expand;
  sxyzeven = spinxyzgen[(creven/.sigma->#)&];
  s2even = inner[sxyzeven, sxyzeven] // Expand;
  MyVPrint[3, "neven:", neven]; (* Report level 3 = nonessntial stuff *)
  MyVPrint[3, "s2even:", s2even];
];

(* Odd combination of creation operators for TWO DOTS. This might be
   useful even for more than two dots, as one of the orthogonal
   combinations of states. *) 
If[NRDOTS >= 2,
  crodd = 1/Sqrt[2] ( d[CR, sigma] - a[CR, sigma] );
  anodd = conj[crodd];
  nodd = Sum[ nc[crodd, anodd], {sigma, 0, 1}] // Expand;
  sxyzodd = spinxyzgen[(crodd/.sigma->#)&];
  s2odd = inner[sxyzodd, sxyzodd] // Expand;
  MyVPrint[3, "nodd:", nodd];
  MyVPrint[3, "s2odd:", s2odd];
  seso = inner[sxyzeven, sxyzodd] // Expand;
  MyVPrint[3, "seso:", seso];
];


(* Expand will factor common operators. We also collect common factors in
advance. *)

H = Expand[H];
MyPrint["Hamiltonian generated. ", H];

(* Is Hamiltonian Hermitian? *)
If[paramdefaultbool["checkHc", True],
  Hc = conj[H];
  Hdiff = Simplify[Expand[H-Hc]];
  MyPrint["H-conj[H]=", Hdiff];
  hookfile["Hcsimpl"];
  If[POL2x2, (* hack *)
    Hdiff = Hdiff /. coefzeta[4,i_]->coefzeta[3,i];
  ];
  Hdiff = Hdiff /. {0. -> 0, Complex[0.,0.] -> 0};
  If[Hdiff =!= 0,
    MyError["Non-Hermitian Hamiltonian!"]
  ];
];

If[DEBUG >= 4,
(* Rewrite the Hamiltonian in terms of high-level functions. This
   makes the debugging of Hamiltonian definitions significantly easier. 
   Note, however, that this operation may take a lot of time! *)
  MyPrint["Hamiltonian -> ", SnegSimplify[H]];
];

(*************************)

(* In NRG, it is customary to define a dimensionless Hamiltonian H_N whose
smallest coefficients are of order 1. To be definite, we rescale the
Hamiltonian so that the hopping matrix elements along the Wilson chain
become precisely 1 in the asymptotic N->inf limit.  This is achieved by
dividing the truncated Hamiltonian with a suitable prefactor which can be
further decomposed into a N-dependant part LAMBDA^(-(N-1)/2), and a constant
prefactor. *)

(* See Yoshida, Whitaker, Oliveira, PRB 41 9403, Eqs. (23), (29), (39-40). *)
If[DY, sc0 = (1+1/lambda)/2 lambda^-(z-1) ]; (* FK *)

(* See Campo, Oliveira, PRB 72, 104432, Eq. (46). *)
If[DC || DZ, sc0 = (1-1/lambda)/Log[lambda] lambda^-(z-1) ];

(* (SC): This corresponds to the function of the same name in nrg.cc ! *)
(* Here we use the rescale factor to account for different bandwidths. *)
SCALE[n_] = sc0 lambda^(-(n-1)/2) bandrescale; (* YYY *)

MyPrint["SCALE[0]=", SCALE[0]];

faktor = 1/SCALE[1];
MyPrint["faktor=", faktor];

(* Lambda correction. 'gammaA' is denoted as A_\Lambda in Krishnamurhty et
al. and other papers. \Gamma^C=A_\Lambda \Gamma and, by analogy,
J_K^C=A_\Lambda J_K, where the superscript C denotes corrected quantity.
This corrections improves convergence to the continuum limit in the original
logarithmic discretization approach (Y). *)

gammaA = If[DY, 1/2 (1+1/lambda)/(1-1/lambda) Log[lambda], 1];


(* Additional error check: are all required parameters defined? *)
(* Strategy: apply parameters, replace operators by 1. Is the result
   a number? *)
checkdefinitions[] := Module[{Htest},
  Htest = H /. params;
  Htest = Htest //. { op_?operatorQ -> 1, bra[__] -> 1, ket[__] -> 1 };
  Htest = N[Htest];
  MyPrint["checkdefinitions[] -> ", Htest];   

  If[!NumberQ[Htest],
    MyError["Undefined paramters detected: ", Htest];
  ];
];
  
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
  needqszlist = {"QSZ", "QSZLR", "ANYJ", "ISOSZ", "ISOSZLR"};
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

  MyPrint["Basis states generated."];
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
  
  If[isANYJ[],
    (* Keep using half-integer quantum numbers! *)
    bvc = bzop2bzvc[bz, vak];
    degnr[{___}] := 1;
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
  MyPrint["BASIS NR=", nrstates1];

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
      (* Maximum spin state *)
      bzspin = {First @ spinbasis[MAKESPINKET]};
      MyPrint["spin basis=", bzspin];

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

    DBLISOSZaddspinketNRG[bvc_List, SP_] := Module[{spinbz},
      spinbz = spinbasis[SP];
      basistensorproduct[spinbz, bvc, Function[{qn1, qn2},
        Module[{x=qn2}, x[[3]] += 2*First[qn1]; x] ]  (* Note the factor 2! *)
      ]
    ];

    If[isDBLISOSZ[],
      bvc = DBLISOSZaddspinketNRG[bvc, MAKESPINKET];
      bz = bzvc2bzop[bvc];
    ];

    If[isANYJ[],
      bvc = QSZaddspinket[bvc, MAKESPINKET];
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

  (**** NEW: Basis states generally do not change for a given (MODEL,
  VARIANT) pair, therefore we can store them to a file and reuse them
  later. This can save considerable amount of time. ****)

  MyPrint["Basis: " <> basisfilename];
  MyPut[bvc, basisfilename];
  
  timeadd["basis"];
];

(* At this point bvc contains the basis in the occupation number
representation. bz is not used directly past this point. *)

bz =.; (* Safety measure! *)


(********************** EXACT DIAGONALIZATION ************************)

(* At this point the basis had been generated (or it had been read
from a file) and we are ready to perform all the necessary
calculations: diagonalization of the Hamiltonian and unitary
transformations of all other operators that we want to evaluate during
the NRG iteration. *)

Scan[(bazavc[ First[#] ] = Last[#])&, bvc];

(* Prepare for optimizations in op2matrix[]: *)

(* getdiffvc[] returns all the vc[] elements which appear in the expression
'm'. This is used to obtain the action of an operator on these elements,
thus enabling the reuse of the results and avoiding recalculations. *)

getdiffvc[m_] := Union[Flatten[(m /. {Plus -> List}) //. x_ v_vc :> v]];

(* transformrule[] creates a mapping between the occupation-number vectors
and a set of orthonomal vectors in the R^n space. This is used to transform
operator expressions for operators into matrix representations. *)

transformrule[vecs_] := Thread[vecs -> IdentityMatrix[Length[vecs]]];

(* bazavcdiffvc[] returns the vc[] elements which appear in the basis
vectors in each invariant subspace. *)

Scan[(bazavcdiffvc[First[#]] = getdiffvc[Last[#]]) &, bvc];

(* bazavctransf[] returns the R^n vector representation of the basis. This
is a (dim x len) matrix, where 'dim' is the number of elements in the basis
and 'len' the number of vc[] elements which are present. *)

Scan[(bazavctransf[First[#]] = 
  Last[#] /. transformrule[bazavcdiffvc[First[#]]]) &, bvc];

(* List of all invariant subspaces, i.e. all possible combinations of the
  quantum numbers such as (Q,2S+1) or (2I+1,2S+1) *)

subspaces = bvc[[All, 1]];
nrsub = Length[subspaces];
ndx2subspace[i_] := Position[subspaces, i][[1,1]];

(* Integer or half-integer spin? Used for "QS". *)
integerspin[ss_] := OddQ[ss];
halfintegerspin[ss_] := EvenQ[ss];

(* We hereby assume that the representative members in each
subspace have maximal spin (and isospin) projection quantum
number, i.e. S_z=S, I_z=I. *)
subspacesz[ss_] := SS2S[ss];
subspaceiz[ii_] := II2I[ii];
subspacejz[jj_] := JJ2J[jj];

(* Note: for orbital moments we do not use the 2T+1 convention, 
thus this one is simply equal to T. *)
subspacetz[t_] := t;

(* Dimension of the invariant subspace *)
dim[inv_] := Length[ bazavc[inv] ];

(* StringJoinSep[] joins strings with a separator between each element. *)
StringJoinSep[sep_, strlist_] := StringJoin @@ 
  (Most @ Flatten @ ({#, sep}& /@ strlist));
StringJoinSep[sep_, {}] := "";

(* Ensure that rational numbers are converted to ASCII fractions with 'f' separators. *)
myToString[n_Rational] := ToString[Numerator[n]] <> "f" <> ToString[Denominator[n]];
myToString[n_] := ToString[n];

(* Convert Invar to a string that is safe to be used as a filename.
Furthermore, the mapping from Invars to strings should be 1-1. *)

Invar2String[inv_] := StringJoinSep[".", myToString /@ inv];

MPVCFAST = If[option["MPVCSLOW"], False, True];

(* Create a matrix representation of the operator 'op'. *)

op2matrix[op_, inv1_, inv2_] := Module[{bz1, bz2, mm, b1, b2},
  (* Step 1: calculate the matrix elements between all the vc[] vectors which
     appear in inv1 and inv2. *)
  bz1 = bazavcdiffvc[inv1];
  bz2 = bazavcdiffvc[inv2];
  
    (* We use ..vcfast[] here. This function is only applicable to simple
    monomial vectors, such as the ones which appear in bz1 and bz2 here. The
    calculation is performed by applying the operator 'op' on elements of
    bz2 and looking up the indexes of resulting elements in bz1.  *)
    
  If[MPVCFAST == True,  
    mm = matrixrepresentationvcfast[op, bz1, bz2],
    mm = matrixrepresentationvc[op, bz1, bz2]
  ];

  (* Step 2: perform the unitary transformation so as to obtain the matrix
     representation of the operator 'op', <inv1|op|inv2>. *)
  b1 = bazavctransf[inv1];
  b2 = bazavctransf[inv2];

  (* For rotation symmetries such as C_3, the basis is complex!  We conjugate the bras, transpose the kets. *)
  Conjugate[b1].mm.Transpose[b2]
];

op2matrix[op_, inv_] := op2matrix[op, inv, inv];

(* Generate the matrix form of the Hamiltonian in a given invariant
subspace. Function ham[] is called from diagvc[]. *)

GENERATEHAM = If[option["READHAM"], False, True];

hamfn[inv_] := hamfilename <> "_" <> Invar2String[inv];

ham[inv_] := Module[{fn, rep},
  timestart["ham"];
  fn = hamfn[inv];

  If[GENERATEHAM == False,
    MyPrint["Reading matrix from " <> fn];
    rep = silentGet[fn, Path -> dumppath];
    MyVPrint[2, "rep=", rep];

    (* Fall back to matrix generation *)
    If[rep === $Failed,
      GENERATEHAM = True;
    ];
  ];
    
  If[GENERATEHAM && !option["PARAMPRE"],
    MyPrint["Generating matrix: " <> fn];
    rep = matrixrepresentationvc[H, bazavc[inv]];

    (* Simplification improves numerical precision! *)
    rep = Simplify[rep];
    MyVPrint[2, "rep=", rep];
  
    (* Save the *generic* Hamiltonian matrix to a file *)
    MyPut[rep, fn, option["GENERATE_TEMPLATE"]];
  ];
    
  (* If option PARAMPRE is specified, apply parameters now! *)
  If[GENERATEHAM && option["PARAMPRE"],
    MyPrint["Generatic numerical matrix ", inv];
    Hnum = H /. params;
    Hnum = Hnum /. 0. -> 0;
    MyVPrint[3, "Hnum=", Hnum];
    rep = matrixrepresentationvc[Hnum, bazavc[inv]];
    (* Don't save to file in this case!!! *)
  ];

  timeadd["ham"];
  rep
];

(* Returns True if 'a' is a numeric matrix *)
isnumericmat[a_] := And @@ Map[NumericQ, Flatten[a]];

(* diagvc[] diagonalises the Hamiltonian in selected subspace. It returns
two lists: eigenvalues (absolute energy units [D], no factors!), and
eigenvectors as numeric arrays (i.e. a matrix). The results are cached. *)

diagvc[inv_] := diagvc[inv] = Module[{hamil, dim, nr, val, vec},
  MyPrint["diagvc[", inv, "]"];

  hamil = ham[inv];
  MyVPrint[1, "hamil=", hamil];

  dim = Dimensions[hamil];
  MyVPrint[1, "dim=", dim];
  MyAssert[dim[[1]] == dim[[2]]];
  nr = dim[[1]];

  If[option["TEMPLATE"] || option["GENERATE_TEMPLATE"],
    (* Fake a result with correct structure *)
    Return[{Range[nr], IdentityMatrix[nr]//N}];
  ];
  
  (* Generation of numeric matrices *)
  hamil = hamil /. params; (* Apply parameters here! *)
  MyVPrint[2, "hamil=", hamil];
  hamil = N[hamil];

  If[!isnumericmat[hamil], (* Additional error trap *)
    MyError["Matrix is not numeric: ", hamil];
  ];

  diff = hamil - ComplexTranspose[hamil];
  If[Total[Abs[diff]] != 0,
    MyError["Non-hermitian Hamiltonian: diff=", diff];
  ];
  
  (* Mathematica's built-in Eigensystem[] has problems with numeric
  matrices having eigenvalues with high multiplicity (or even matrices
  with nearly degenerate eigenvalues, leading to lack of orthogonality
  between the eigenvectors), whereas SchurDecomposition[] works just
  fine in such instances! As of 28.7.2011, the Schur decomposition
  is the default behavior and it can be turned off using the NOSCHUR 
  option.*)
  
  If[option["NOSCHUR"],
    eigsys = Eigensystem[hamil];
    eigsys = N[eigsys]; (* Enforce numericity! *)
  
    (* Remove (infinitesimal) imaginary parts. *)
    {val, vec} = eigsys;
    val = Re[val],

  (* else *)

    (* Perform Schur decomposition (QR) *)   
    {Sq, St} = SchurDecomposition[hamil];
    val = Re @ Diagonal[St];
    vec = Transpose[Sq];
  ];

  eigsys = {val, vec};
     
  (* Sort in the ascending order of eigenvalues! *)
  {val, vec} = Transpose[ Sort[ Transpose[eigsys] ] ];

  If[DEBUG >= 3,
    MyPrint["val=", val];
  ];
  
  (* Sanity check 1 *)
  vecnormnonzero = Map[(Norm[#] != 0)&, vec];
  If[And @@ vecnormnonzero =!= True,
    MyError["Zero eigenvector detected."];
  ];
    
  (* Normalize eigenvectors *)
  vec = Map[#/Norm[#]&, vec];

  (* Sanity check 2 *)
  DETLIMIT = 10^-10;
  det = Det[vec];
  detdiff = 1.0-Abs[det];
  MyPrint["det[vec]=", det, " 1-abs=", detdiff];
  If[Abs[detdiff] > DETLIMIT, MyError["Det[] error"]];
    
  (* Sanity check 3 *)
  ORTHLIMIT = 10^-10 * nr; (* Don't be too stringent! Rescale by matrix size! *)
  res = Total[Abs[Dot[vec, 
    ConjugateTranspose[vec]]-IdentityMatrix[nr]], 2];
  MyPrint["orthogonality check=", res];
  If[Abs[res] > ORTHLIMIT, MyError["orhogonality error"]];

  {val, vec}
];


(* coupledQ? Are two subspaces coupled by creation operator? *)
(* spincoupledQ? Are two subspaces coupled by spin (triplet) operator? *)

coupledQ["QS" | "QSLR", {{q1_, ss1_, i1___}, {q2_, ss2_, i2___}}] :=
  If[q1 == q2+1 && Abs[ss1-ss2] == 1, True, False, MyError["oops"]];
spincoupledQ["QS" | "QSLR", {{q1_, ss1_, i1___}, {q2_, ss2_, i2___}}] := 
  If[q1 == q2 && (ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2), True, False,
    MyError["oops"]];

(* f[p] of three types. Need finer grained checks. *)
coupledQ["QSC3", {{q1_, ss1_, p1_, i1___}, {q2_, ss2_, p2_, i2___}}] :=
  If[q1 == q2+1 && Abs[ss1-ss2] == 1, True, False, MyError["oops"]];
spincoupledQ["QSC3", {{q1_, ss1_, p1_, i1___}, {q2_, ss2_, p2_, i2___}}] := 
  If[q1 == q2 && (ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2) && p1 == p2, True, False,
    MyError["oops"]];

(* Since f is an orbital triplet operator, the case t1=t2=0 counts as not coupled. *)
coupledQ["QST", {{q1_, ss1_, t1_, i1___}, {q2_, ss2_, t2_, i2___}}] :=
  If[q1 == q2+1 && Abs[ss1-ss2] == 1 && (Abs[t1-t2] <= 1 && !(t1==0 && t2==0)), True, False, MyError["oops"]];
spincoupledQ["QST", {{q1_, ss1_, t1_, i1___}, {q2_, ss2_, t2_, i2___}}] := 
  If[q1 == q2 && (ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2) && t1 == t2, True, False,
    MyError["oops"]];
orbcoupledQ["QST", {{q1_, ss1_, t1_, i1___}, {q2_, ss2_, t2_, i2___}}] := 
  If[q1 == q2 && (t1 == t2 || t1 == t2-1 || t1 == t2+1) && ss1 == ss2, True, False,
    MyError["oops"]];

coupledQ["QSTZ", {{q1_, ss1_, tz1_, i1___}, {q2_, ss2_, tz2_, i2___}}] :=
  If[q1 == q2+1 && Abs[ss1-ss2] == 1 && Abs[tz1-tz2] <= 1, True, False, MyError["oops"]];
spincoupledQ["QSTZ", {{q1_, ss1_, tz1_, i1___}, {q2_, ss2_, tz2_, i2___}}] := 
  If[q1 == q2 && (ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2) && tz1 == tz2, True, False,
    MyError["oops"]];

coupledQ["QSZTZ", {{q1_, ssz1_, tz1_, i1___}, {q2_, ssz2_, tz2_, i2___}}] :=
  If[q1 == q2+1 && Abs[ssz1-ssz2] == 1 && Abs[tz1-tz2] <= 1, True, False, MyError["oops"]];
spincoupledQ["QSTZ", {{q1_, ssz1_, tz1_, i1___}, {q2_, ssz2_, tz2_, i2___}}] := 
  If[q1 == q2 && (ssz1 == ssz2 || ssz1 == ssz2-2 || ssz1 == ssz2+2) && tz1 == tz2, True, False,
    MyError["oops"]];

(* Depends on j! *)
coupledQ["QJ", {{q1_, jj1_}, {q2_, jj2_}}, j_] :=
  If[q1 == q2+1 && Abs[jj1-jj2] <= 2j && Abs[JJ2J[jj1]-j] <= JJ2J[jj2] && Abs[JJ2J[jj2]-j] <= JJ2J[jj1], True, False, MyError["oops"]];
coupledQ["QJ", _] := False;

(* Since f is an orbital triplet operator, the case t1=t2=0 counts as not coupled. *)
coupledQ["SPSU2T", {{ss1_, t1_, i1___}, {ss2_, t2_, i2___}}] :=
  If[Abs[ss1-ss2] == 1 && (Abs[t1-t2] <= 1 && !(t1==0 && t2==0)), True, False, MyError["oops"]];
spincoupledQ["SPSU2T", {{ss1_, t1_, i1___}, {ss2_, t2_, i2___}}] := 
  If[(ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2) && t1 == t2, True, False,
    MyError["oops"]];

coupledQ["SPSU2" | "SPSU2LR" | "SPSU2C3", {{ss1_, i1___}, {ss2_, i2___}}] := 
  If[Abs[ss1-ss2] == 1, True, False, MyError["oops"]];
(* TODO: parity should not change! *)
spincoupledQ["SPSU2" | "SPSU2LR" | "SPSU2C3", {{ss1_, i1___}, {ss2_, i2___}}] := 
  If[ss1 == ss2 || ss1 == ss2-2 || ss1 == ss2+2, True, False,
    MyError["oops"]];

coupledQ["SPU1" | "SPU1LR", {{ssz1_, i1___}, {ssz2_, i2___}}] := 
  If[Abs[ssz1-ssz2] == 1, True, False, MyError["oops"]];
spincoupledQ["SPU1" | "SPU1LR", {{ssz1_, i1___}, {ssz2_, i2___}}] := 
  If[ssz1 == ssz2 || ssz1 == ssz2-2 || ssz1 == ssz2+2, True, False,
     MyError["oops"]];

coupledQ["ISO" | "ISOLR" | "ISO2" | "ISO2LR", 
  {{ii1_, ss1_, j1___}, {ii2_, ss2_, j2___}}] := 
  If[Abs[ii1-ii2] == 1 && Abs[ss1-ss2] == 1, True, False,
     MyError["oops"]];
    
coupledQ["ISOSZ" | "ISOSZLR", {{ii1_, ssz1_, j1___}, {ii2_, ssz2_, j2___}}] := 
  If[Abs[ii1-ii2] == 1 && Abs[ssz1-ssz2] == 1, True, False,
     MyError["oops"]];

coupledQ["SU2", {{ii1_, j1___}, {ii2_, j2___}}] := 
  If[Abs[ii1-ii2] == 1, True, False, MyError["oops"]];

coupledQ["DBLSU2", {{ii11_, ii21_, j1___}, {ii12_, ii22_, j2___}}] := 
  If[(Abs[ii11-ii12] == 1) ~ Xor ~ (Abs[ii21-ii22] == 1),
  True, False, MyError["oops coupledQ"]];

coupledQ["DBLISOSZ", {{ii11_, ii21_, ssz1_, j1___}, {ii12_, ii22_, ssz2_, j2___}}] := 
  If[((Abs[ii11-ii12] == 1 && Abs[ii21-ii22] == 0) ~ Xor ~ (Abs[ii21-ii22] == 1 && Abs[ii11-ii12] == 0)) && 
  Abs[ssz1-ssz2] == 1,
  True, False, MyError["oops coupledQ"]];

coupledQ["QSZ" | "QSZLR", {{q1_, ssz1_, i1___}, {q2_, ssz2_, i2___}}] := 
  If[q1 == q2+1 && Abs[ssz1-ssz2] == 1, True, False, MyError["oops"]];

coupledQ["P", {{p1_}, {p2_}}] := (p1 != p2); (* Opposite parity *)

coupledQ["PP", {{pa1_, pb1_}, {pa2_, pb2_}}] := ((pa1 != pa2) && (pb1 == pb2)) || ((pa1 == pa2) && (pb1 != pb2));

coupledQ["NONE", {_, _}] := True;

coupledQ["U1", {{q1_, i1___}, {q2_, i2___}}] := 
  If[q1 == q2+1, True, False, MyError["oops"]];

coupledQ["SL", {{q1_, i1___}, {q2_, i2___}}] := 
  If[q1 == q2+1, True, False, MyError["oops"]];

coupledQ["SL3", {{q11_, q21_, q31_, i1___}, {q12_, q22_, q32_, i2___}}] := 
  If[
   (q11 == q12+1 && q21 == q22 && q31 == q32) ~ Xor ~
   (q11 == q12 && q21 == q22+1 && q31 == q32) ~ Xor ~
   (q11 == q12 && q21 == q22 && q31 == q32+1),
  True, False, MyError["oops"]];
 
coupledQ["ANYJ", {{q1_, ssz1_, i1___}, {q2_, ssz2_, i2___}}] := 
 If[q1 == q2+1 && MemberQ[Table[spi, {spi, -BANDSPIN, BANDSPIN}], ssz1-ssz2], 
  True, False, MyError["oops"]];
 
(* Bug trap *)
coupledQ[s_String, a___] := MyError["coupledQ not defined for ", {s,a}];
(* spincoupledQ[String_, ___] := MyError["spincoupledQ not defined"]; *)

(* ======================= *)

Flatten1[list_] := Flatten[list, 1];

(* All pairs of subspaces *)
subspacepairs = Reverse @ Flatten1 @ Table[{subspaces[[i]], subspaces[[j]]}, 
  {i, nrsub}, {j, nrsub}];
nrp = Length[subspacepairs];

(* Subspaces coupled by creation operator *)
coupledpairs = Select[subspacepairs, coupledQ[SYMTYPE, #]&];
nrcp = Length[coupledpairs];

MyVPrint[3, "nrp=", nrp];
MyVPrint[3, "subspacepairs=", subspacepairs];

MyVPrint[3, "nrcp=", nrcp];
MyVPrint[3, "coupledpairs=", coupledpairs];

(* Subspaces coupled by spin operator *)
spincoupledpairs = Select[subspacepairs, spincoupledQ[SYMTYPE, #]&];
nrspincp = Length[spincoupledpairs];

(* Subspaces coupled by orbital operator *)
orbcoupledpairs = Select[subspacepairs, orbcoupledQ[SYMTYPE, #]&];
nrorbcp = Length[orbcoupledpairs];



(*********** Calculation of various matrix elements ************)

(* CONVENTIONS:

 [name] - routine for a single element calculation. Mostly for testing.
 [name]Matrix - create a matrix of operator in a given subspace or 
                combination of subspaces 
 [name]Table - call [name]Matrix for every subspace or combination of 
               subspaces, used to build an exportable table of all 
               required matrix elements
*)

(*                     --- DOUBLET OPERATORS ---               *)

(* 
 Call hierarchy: 
 maketable -> ireducTable -> ireducMatrixSpeedy -> optransform -> op2matrix
*)

(* Note: for QSZ, both spin projections are considered. They are output as a
single block: the spin projection is defined by the difference ssz1-ssz2 of 
the two invariant subspaces involved. The C++ part of the program should
disentangle the spectral densities for the separate spin projections. *)

(* inv1 and inv2 are the sets of quantum numbers corresponding
   to an invariant subspace. The actual eigenstates are then
   stored in diagvc[inv1] and diagvc[inv2]. *)


(* Compatibility wrapper. #1 will be replaced by CR/AN (isospin)
   index, while #2 will be replaced by UP/DO (spin) index. *)

ireducMatrixSpeedy[str_String, op_?operatorQ[j___], {inv1_, inv2_}, opt___] :=
  ireducMatrixSpeedy[str, op[#1, j, #2]&, {inv1, inv2}, opt];

(* optransform[] transforms an operator to a matrix representation
<inv1|op|inv2>/factor and rotates it by a suitable unitary transformation to
the eigenbases in inv1 and inv2. *)

optransform[op_, inv1_, inv2_, factor_:1] := Module[{mat, vecs1, vecs2, res},
  mat = op2matrix[op, inv1, inv2] / factor;
  vecs1 = Conjugate[ diagvc[inv1] [[2]] ]; (* bras are conjugated *)
  vecs2 = Transpose[ diagvc[inv2] [[2]] ]; (* kets are transposed *)

  vecs1 . mat . vecs2
];

udf[szop_] := If[szop == 1/2, UP, DO, MyError["oops, szop=", szop]];
duf[szop_] := If[szop == 1/2, DO, UP, MyError["oops, szop=", szop]];
signf[szop_] := If[szop == 1/2, +1, -1, MyError["oops, szop=", szop]];

getIsospinQN["ISO" | "ISOLR" | "ISO2" | "ISO2LR" | "ISOSZ" | "ISOSZLR" | 
             "SU2", inv_] := inv[[1]];

getSpinQN["QS" | "QSLR" | "QSC3" | "QSZ" | "QSZLR" | "ISO" | "ISOLR" | 
          "ISO2" | "ISO2LR" | "ISOSZ" | "ISOSZLR" | "QST" | "QSTZ" | "QSZTZ", inv_] := inv[[2]];
getSpinQN["SPSU2" | "SPSU2LR" | "SPSU2T" | "SPU1" | "SPU1LR" | "SPSU2C3", inv_] := inv[[1]];

getOrbitalQN["QST" | "QSTZ" | "QSZTZ", inv_] := inv[[3]];
getOrbitalQN["SPSU2T", inv_] := inv[[2]];

getJQN["QJ", inv_] := inv[[2]];

ireducMatrixSpeedy[symtype:("QSZ" | "QSZLR" | "SPU1" | "SPU1LR"),
                   op_, {inv1_, inv2_}, ___] := 
Module[
  {ssz1, ssz2, szop, ud, op1},
  ssz1 = getSpinQN[symtype, inv1];
  ssz2 = getSpinQN[symtype, inv2];
  szop = 1/2(ssz1-ssz2);
  ud = udf[szop];
  op1 = op[CR, ud]; (* IMPORTANT: the operator is a creation operator
    with the spin determined by the invariant subspaces linked by it. *)
  optransform[op1, inv1, inv2]
];

ireducMatrixSpeedy[symtype:("QS" | "QSLR" | "QSC3" | "SPSU2" | "SPSU2LR" | "SPSU2C3"),
                   op_, {inv1_, inv2_}, ___] :=
Module[{ss1, ss2, sz1, sz2, szop, op1, ud, cg},
  ss1 = getSpinQN[symtype, inv1];
  ss2 = getSpinQN[symtype, inv2];
  (* subspacesz[] returns the S_z for given S. The representative 
     states in the basis have this value of S_z. *)
  sz1 = subspacesz[ss1];
  sz2 = subspacesz[ss2];
  szop = sz1-sz2;
  ud = udf[szop];
  op1 = op[CR, ud]; (* IMPORTANT: CR by default! *)
  cg = ClebschGordan[{SS2S[ss2],sz2},{1/2,szop},{SS2S[ss1],sz1}];
  optransform[op1, inv1, inv2, cg]
];

(* Two types of f operators, j=3/2 and j=1/2 *)
ireducMatrixSpeedy[symtype:("QJ"),
                   op_, {inv1_, inv2_}, j_, ___] :=
Module[{jj1, jj2, jz1, jz2, jzop, op1, cg},
  jj1 = getJQN[symtype, inv1];
  jj2 = getJQN[symtype, inv2];
  jz1 = subspacejz[jj1];
  jz2 = subspacejz[jj2];
  jzop = jz1-jz2;
  op1 = op[CR, jzop];
  cg = ClebschGordan[{JJ2J[jj2],jz2},{j,jzop},{JJ2J[jj1],jz1}];
  optransform[op1, inv1, inv2, cg]
];

ireducMatrixSpeedy[symtype:("QST" | "SPSU2T"),
                   op_, {inv1_, inv2_}, ___] :=
Module[{ss1, ss2, sz1, sz2, szop, op1, ud, t1, t2, tzop, cg1, cg2},
  ss1 = getSpinQN[symtype, inv1];
  ss2 = getSpinQN[symtype, inv2];
  (* subspacesz[] returns the S_z for given S. The representative 
     states in the basis have this value of S_z. *)
  sz1 = subspacesz[ss1];
  sz2 = subspacesz[ss2];
  szop = sz1-sz2;
  ud = udf[szop];
  
  t1 = getOrbitalQN[symtype, inv1];
  t2 = getOrbitalQN[symtype, inv2];
  tzop = t1-t2;
  
  op1 = op[CR, tzop, ud]; (* IMPORTANT: CR by default! *)
  cg1 = ClebschGordan[{SS2S[ss2],sz2},{1/2,szop},{SS2S[ss1],sz1}];
  cg2 = ClebschGordan[{t2,t2},{1,tzop},{t1,t1}];
  optransform[op1, inv1, inv2, cg1 cg2]
];

ireducMatrixSpeedy[symtype:("QSTZ"),
                   op_, {inv1_, inv2_}, ___] :=
Module[{ss1, ss2, sz1, sz2, szop, op1, ud, tz1, tz2, tzop, cg1},
  ss1 = getSpinQN[symtype, inv1];
  ss2 = getSpinQN[symtype, inv2];
  (* subspacesz[] returns the S_z for given S. The representative 
     states in the basis have this value of S_z. *)
  sz1 = subspacesz[ss1];
  sz2 = subspacesz[ss2];
  szop = sz1-sz2;
  ud = udf[szop];
  
  tz1 = getOrbitalQN[symtype, inv1];
  tz2 = getOrbitalQN[symtype, inv2];
  tzop = tz1-tz2;
  
  op1 = op[CR, tzop, ud]; (* IMPORTANT: CR by default! *)
  cg1 = ClebschGordan[{SS2S[ss2],sz2},{1/2,szop},{SS2S[ss1],sz1}];
  optransform[op1, inv1, inv2, cg1]
];

ireducMatrixSpeedy[symtype:("QSZTZ"),
                   op_, {inv1_, inv2_}, ___] :=
Module[{ssz1, ssz2, szop, op1, ud, tz1, tz2, tzop, cg1},
  ssz1 = getSpinQN[symtype, inv1];
  ssz2 = getSpinQN[symtype, inv2];
  (* subspacesz[] returns the S_z for given S. The representative 
     states in the basis have this value of S_z. *)
  szop = 1/2(ssz1-ssz2);
  ud = udf[szop];
  
  tz1 = getOrbitalQN[symtype, inv1];
  tz2 = getOrbitalQN[symtype, inv2];
  tzop = tz1-tz2;
  
  op1 = op[CR, tzop, ud]; (* IMPORTANT: CR by default! *)
  optransform[op1, inv1, inv2]
];

(* f[] operators are singlets wrt the (trivial) symmetry group. *)
ireducMatrixSpeedy["NONE" | "P" | "PP", op_, {inv1_, inv2_}, fnr_, ___] := Module[{op1},
  op1 = Switch[fnr, 
    0, op[CR, DO],
    1, op[CR, UP],
    2, op[AN, DO],
    3, op[AN, UP],
    _, MyError["Critical error in ireducMatrixSpeedy[]"]
  ];
  optransform[op1, inv1, inv2]
];

(* IMPORTANT: two different matrices are required, one for each spin
   orientation! *)
ireducMatrixSpeedy["U1", op_, {inv1_, inv2_}, spin_, ___] := Module[
{op1, mat},
  op1 = op[CR, spin]; (* Recall: DO=0, UP=1 *)
  mat = optransform[op1, inv1, inv2]
];

ireducMatrixSpeedy["SL", op_, {inv1_, inv2_}, ___] := Module[{op1},
  op1 = op[CR, UP]; (* Spinless = only spin-up electrons retained. *)
  optransform[op1, inv1, inv2]
];

ireducMatrixSpeedy["SL3", op_, {inv1_, inv2_}, ___] := Module[{op1},
  op1 = op[CR, UP]; (* Spinless = only spin-up electrons retained. *)
  optransform[op1, inv1, inv2]
];

ireducMatrixSpeedy["ANYJ", op_, {inv1_, inv2_}, ___] := Module[
  {ssz1, ssz2, szop, ud, op1},
  ssz1 = inv1[[2]];
  ssz2 = inv2[[2]];
  szop = ssz1-ssz2;
  If[BANDSPIN == 1/2,
    ud = udf[szop];
    op1 = op[CR, ud],
  (* else *)
    op1 = op[CR, szop],
  (* *)
    MyError["oops"]
  ];
  optransform[op1, inv1, inv2]
];

(* Code audited: Rok, 22. 8. 2006 *)
ireducMatrixSpeedy[symtype:("ISO" | "ISOLR" | "ISO2" | "ISO2LR"),
                   op_, {inv1_, inv2_}, 
                   optnn___] := 
Module[{ii1, ii2, iz1, iz2, ss1, ss2, sz1, sz2, szop, 
        op1, ud, tip, cgs, cgi, faktor},
  ii1 = getIsospinQN[symtype, inv1];
  ii2 = getIsospinQN[symtype, inv2];
  iz1 = subspaceiz[ii1];
  iz2 = subspaceiz[ii2];
  ss1 = getSpinQN[symtype, inv1];
  ss2 = getSpinQN[symtype, inv2];
  sz1 = subspacesz[ss1];
  sz2 = subspacesz[ss2];

  (* If iz1 = iz2 + 1/2, a particle was created (CR) *)
  izop = iz1-iz2;
  tip = If[izop == 1/2, CR, AN, MyError["oops"]];

  (* If a particle with spin up is *created*, sz1 = sz2 + 1/2, if
  a particle wiht spin down is *annihilated*, sz1 = sz2 - 1/2 ! *)
  szop = If[izop == 1/2, sz1-sz2, sz2-sz1, MyError["oops"]];
  ud = udf[szop];

  op1 = op[tip, ud];

  (* Note: iz1-iz2=izop always! *)
  cgi = ClebschGordan[{II2I[ii2], iz2}, {1/2, iz1-iz2}, {II2I[ii1], iz1}];

  (* Note: sz1-sz2=szop if izop=1/2, otherwise sz1-sz2=-szop !! *)
  cgs = ClebschGordan[{SS2S[ss2], sz2}, {1/2, sz1-sz2}, {SS2S[ss1], sz1}];

  (* Alternating index for isospin basis *)
  If[{optnn} === {},
    nn = nnop[ ReleaseHold[
      op[HoldComplete[Sequence[]], HoldComplete[Sequence[]]] ]], (* HACK! *)
    nn = optnn;
  ];
    
  (* In the Nambu spinor, the isospin-up component (creation operator) has
  always the same phase, while the isospin-down component (annihilation
  operator) has an alternating phase, dependant on the position of the site
  in the (bipartite) lattice. *)
  (* NOTE: mind the (-2szop) faktor in the expression for the irreducible
  matrix elements!! *)
  faktor = If[izop == 1/2, 1, (-1)^nn (-2szop), MyError["oops"]];

  optransform[op1, inv1, inv2, faktor * cgs * cgi]
];


ireducMatrixSpeedy[symtype:"ISOSZ" | "ISOSZLR", op_, {inv1_, inv2_}, optnn___] := 
Module[{ii1, ii2, iz1, iz2, ssz1, ssz2, szop, 
        op1, ud, tip, cgs, cgi, faktor},
  ii1 = getIsospinQN[symtype, inv1];
  ii2 = getIsospinQN[symtype, inv2];
  iz1 = subspaceiz[ii1];
  iz2 = subspaceiz[ii2];
  ssz1 = getSpinQN[symtype, inv1];  
  ssz2 = getSpinQN[symtype, inv2];
  szop = 1/2(ssz1-ssz2);

  (* If iz1 = iz2 + 1/2, a particle was created (CR) *)
  izop = iz1-iz2;
  tip = If[izop == 1/2, CR, AN, MyError["oops1"]];

  (* If a particle with spin up is *created*, sz1 = sz2 + 1/2, if
  a particle wiht spin down is *annihilated*, sz1 = sz2 - 1/2 ! *)
  szop = If[izop == 1/2, szop, -szop, MyError["oops2"]];
  ud = udf[szop];

  op1 = op[tip, ud];

  (* Note: iz1-iz2=izop always! *)
  cgi = ClebschGordan[{II2I[ii2], iz2}, {1/2, iz1-iz2}, {II2I[ii1], iz1}];

  (* Alternating index for isospin basis *)
  If[{optnn} === {},
    nn = nnop[ ReleaseHold[
      op[HoldComplete[Sequence[]], HoldComplete[Sequence[]]] ]], (* HACK! *)
    nn = optnn;
  ];
    
  (* In the Nambu spinor, the isospin-up component (creation operator) has
  always the same phase, while the isospin-down component (annihilation
  operator) has an alternating phase, dependant on the position of the site
  in the (bipartite) lattice. *)
  (* NOTE: mind the (-2szop) faktor in the expression for the irreducible
  matrix elements!! *)
  faktor = If[izop == 1/2, 1, (-1)^nn (-2szop), MyError["oops"]];

  optransform[op1, inv1, inv2, faktor * cgi]
];

(* For symtype=SU2, there are two different f operator doublets,
distinguished by the value of 'type', taking two values (UP,DO). One doublet
consists of the [f^dag_UP, f_DO] (szop=1/2) pair, the other of the
[f^dag_DO, -f_UP] (szop=-1/2) pair of operators. Note the sign! *)

ireducMatrixSpeedy[symtype:"SU2", 
                   op_, 
                   {inv1_, inv2_}, 
                   type_, 
                   optnn___] := 
 Module[{ii1, ii2, iz1, iz2, izop, tip, spin, szop, op1, ud, cgs, cgi, faktor},
  ii1 = getIsospinQN[symtype, inv1];
  ii2 = getIsospinQN[symtype, inv2];
  iz1 = subspaceiz[ii1];
  iz2 = subspaceiz[ii2];

  (* If iz1 = iz2 + 1/2, i.e. izop=1/2, a particle was created (CR) *)
  izop = iz1-iz2;
  tip = If[izop == 1/2, CR, AN, MyError["oops1"]];

  If[type == UP, (* cf. matrixnihop1 in su2.m *)
    spin = If[tip == CR, UP, DO]
  ];
  If[type == DO, (* cf. matrixnihop2 in su2.m *)
    spin = If[tip == CR, DO, UP]
  ];

  (* Spin quantum number of the operator 'op1' *)
  szop = If[spin == UP, 1/2, -1/2] * If[izop == 1/2, +1, -1];

  op1 = op[tip, spin];

  (* Note: iz1-iz2=izop always! *)
  cgi = ClebschGordan[{II2I[ii2], iz2}, {1/2, iz1-iz2}, {II2I[ii1], iz1}];

  (* Alternating index for isospin basis *)
  If[{optnn} === {},
    nn = nnop[ ReleaseHold[
      op[HoldComplete[Sequence[]], HoldComplete[Sequence[]]] ]], (* HACK! *)
    nn = optnn;
  ];
    
  (* Compare with the ISOSZ symmetry type. The difference is that szop is
  not determined by the quantum numbers of the invariant subspaces involved,
  but by the parameter type=1,2. *)

  faktor = If[izop == 1/2, 1, (-1)^nn (2szop), MyError["oops"]];

  optransform[op1, inv1, inv2, faktor * cgi]
];


ireducMatrixSpeedy[symtype:"DBLSU2", 
                   op_, 
                   {inv1_, inv2_}, 
                   type_, 
                   optnn___] := 
 Module[{ii11, ii12, iz1, iz2, ii21, ii22, ii1, ii2,
         izop, tip, spin, szop, op1, ud, cgs, cgi, faktor},
  ii11 = inv1[[1]]; (* first index: channel; second index: 1st or 2nd arg *)
  ii21 = inv1[[2]];
  ii12 = inv2[[1]];
  ii22 = inv2[[2]];
  If[Abs[ii11-ii12] == 1, (* f in ch 1 *)
    iz1 = subspaceiz[ii11];
    iz2 = subspaceiz[ii12];
    ii1 = ii11;
    ii2 = ii12;
  ];
  If[Abs[ii21-ii22] == 1, (* f in ch 2 *)
    iz1 = subspaceiz[ii21];
    iz2 = subspaceiz[ii22];
    ii1 = ii21;
    ii2 = ii22;
  ];

  (* If iz1 = iz2 + 1/2, i.e. izop=1/2, a particle was created (CR) *)
  izop = iz1-iz2;
  tip = If[izop == 1/2, CR, AN, MyError["oops1"]];

  If[type == UP, (* cf. matrixnihop1 in su2.m *)
    spin = If[tip == CR, UP, DO]
  ];
  If[type == DO, (* cf. matrixnihop2 in su2.m *)
    spin = If[tip == CR, DO, UP]
  ];

  (* Spin quantum number of the operator 'op1' *)
  szop = If[spin == UP, 1/2, -1/2] * If[izop == 1/2, +1, -1];

  op1 = op[tip, spin];

  (* Note: iz1-iz2=izop always! *)
  cgi = ClebschGordan[{II2I[ii2], iz2}, {1/2, iz1-iz2}, {II2I[ii1], iz1}];

  (* Alternating index for isospin basis *)
  If[{optnn} === {},
    nn = nnop[ ReleaseHold[
      op[HoldComplete[Sequence[]], HoldComplete[Sequence[]]] ]], (* HACK! *)
    nn = optnn;
  ];
    
  (* Compare with the ISOSZ symmetry type. The difference is that szop is
  not determined by the quantum numbers of the invariant subspaces involved,
  but by the parameter type=1,2. *)

  faktor = If[izop == 1/2, 1, (-1)^nn (2szop), MyError["oops"]];

  optransform[op1, inv1, inv2, faktor * cgi]
];

ireducMatrixSpeedy[symtype:"DBLISOSZ", 
                   op_, 
                   {inv1_, inv2_}, 
                   optnn___] := 
 Module[{ii11, ii12, iz1, iz2, ii21, ii22, ii1, ii2, ssz1, ssz2,
         izop, tip, spin, szop, op1, ud, cgs, cgi, faktor},
  ii11 = inv1[[1]]; (* first index: channel; second index: 1st or 2nd arg *)
  ii21 = inv1[[2]];
  ssz1 = inv1[[3]];
  ii12 = inv2[[1]];
  ii22 = inv2[[2]];
  ssz2 = inv2[[3]];
  
  ii1 = Indeterminate;
  If[Abs[ii11-ii12] == 1 && Abs[ii21-ii22] == 0, (* f in ch 1 *)
    iz1 = subspaceiz[ii11];
    iz2 = subspaceiz[ii12];
    ii1 = ii11;
    ii2 = ii12;
  ];
  If[Abs[ii21-ii22] == 1 && Abs[ii11-ii12] == 0, (* f in ch 2 *)
    iz1 = subspaceiz[ii21];
    iz2 = subspaceiz[ii22];
    ii1 = ii21;
    ii2 = ii22;
    ];
  If[ii1 === Indeterminate,
    MyPrint["inv1=", inv1, " inv2=", inv2];  
    MyError["ireduc[] error."];
  ];
    
  szop = 1/2(ssz1-ssz2);

  (* If iz1 = iz2 + 1/2, i.e. izop=1/2, a particle was created (CR) *)
  izop = iz1-iz2;
  tip = If[izop == 1/2, CR, AN, MyError["oops1. izop=", izop]];

  (* If a particle with spin up is *created*, sz1 = sz2 + 1/2, if
  a particle wiht spin down is *annihilated*, sz1 = sz2 - 1/2 ! *)
  szop = If[izop == 1/2, szop, -szop, MyError["oops2"]];
  ud = udf[szop];

  op1 = op[tip, ud];

  (* Note: iz1-iz2=izop always! *)
  cgi = ClebschGordan[{II2I[ii2], iz2}, {1/2, iz1-iz2}, {II2I[ii1], iz1}];
  If[cgi === Indeterminate || cgi === 0,
    MyPrint["inv1=", inv1, " inv2=", inv2];
    MyPrint["I2=", II2I[ii2]];
    MyPrint["I1=", II2I[ii1]];
    MyError["Aborting. cgi=", cgi];
  ];

  (* Alternating index for isospin basis *)
  If[{optnn} === {},
    nn = nnop[ ReleaseHold[
      op[HoldComplete[Sequence[]], HoldComplete[Sequence[]]] ]], (* HACK! *)
    nn = optnn;
  ];
    
  faktor = If[izop == 1/2, 1, (-1)^nn (-2szop), MyError["oops"]];

  optransform[op1, inv1, inv2, faktor * cgi]
];

(* Create a table of irreducible matrix elements of a doublet tenzor
   operator. The Table returned is suitable for inclusion in the 'data' file
   for NRG. *)

(* This routine is used both for f[0], f[1],... for the Wilson chain site,
which is used in the NRG iteration to determine the eigenstates, and also
for spectral densities (i.e. d[]). For debugging purposes it should thus be
noted that if the eigenvalue flows are correct, but the spectra seem wrong,
then the corresponding ireducMatrixSpeedy routine works properly and the bug
should be sought after elsewhere (such as i. coefficients, ii. recalculation
code, iii. specdens_factor() routine). *)
   
ireducTable[op_, 
            optional___] :=  (* optional is passed to ireducMatrixSpeedy[] *)
Module[{t, cp, i},
  MyPrint["ireducTable: ", op, {optional}];
  t = {{nrcp}};
  For[i = 1, i <= nrcp, i++,
    (* coupledpairs is a list of subspace pairs that are coupled
       by doublet operators [that increase charge, when charge conservation 
       is explicitly taken into account], i.e. creation operators! *)
    cp = coupledpairs[[i]];
    AppendTo[t, Flatten[cp]];
    t = Join[t, ireducMatrixSpeedy[SYMTYPE, op, cp, optional] ];
  ];
  t (* Return *)
];


(* Second version: ops is now a list of pairs {weight, op} and we calculate
the doublet irreducible matrix elements for an orbital described by a linear
combinations of "atomic" orbitals. The appropriate normalization is
automatically computed. *)

ireducTable[ops_List] := Module[{norm, t, cp, i},
  MyPrint[ops];
  norm = Norm @ ops[[All,1]];
  t = {};
  AppendTo[t, {nrcp}];
  For[i = 1, i <= nrcp, i++,
    cp = coupledpairs[[i]];
    AppendTo[t, Flatten[cp]];
    t = Join[t, Plus @@ Map[#[[1]]/norm *
      ireducMatrixSpeedy[SYMTYPE, #[[2]], cp] &, ops] ];
  ];
  t (* Return *)
];


(*                     --- TRIPLET OPERATORS ---               *)

(* Irreducible matrix elements of triplet tenzor operator op *)
ireducsigma[SYMTYPE:("SPSU2" | "QS" | "QSLR" | "QSC3" | "QST" | "QSTZ" | "SPSU2T" | "SPSU2C3"),
 op_, {inv1_, inv2_}] := 
  Module[{ss1, ss2, sz1, sz2, szop, op1, xa, xb},
   ss1 = getSpinQN[SYMTYPE, inv1];
   ss2 = getSpinQN[SYMTYPE, inv2];
   sz1 = subspacesz[ss1];
   sz2 = subspacesz[ss2];
   szop = sz1-sz2;

   (* Spherical operators! *)
   If[szop == 0, op1 = spinz[op] ];
   If[szop == 1, op1 = -1/Sqrt[2] spinplus[op] ];
   If[szop == -1, op1 = 1/Sqrt[2] spinminus[op] ];

   xa = optransform[op1, inv1, inv2];
   xb = If[ss1 == 1 && ss2 == 1, 1, 
     ClebschGordan[{SS2S[ss2],sz2},{1,szop},{SS2S[ss1],sz1}]];
   xa/xb
];

(* Bug trap *)
ireducsigma[SYMTYPE_, args___] := 
  MyError["ireducsigma not implemented for " <> SYMTYPE <> 
  " args=" <> ToString[{args}]];

(* Table form for all irreducible matrix elements of a triplet operator sigma *)
ireducsigmaTable[op_] := Module[{t, i, cp},
  t = {};
  AppendTo[t, {nrspincp}];
  For[i = 1, i <= nrspincp, i++,
    cp = spincoupledpairs[[i]];
    AppendTo[t, Flatten[cp]];
    t = Join[t,  ireducsigma[SYMTYPE, op, cp]];
  ];
  t (* return *)
];


(*                     --- ORBITAL TRIPLET OPERATORS ---               *)

ireducorbsigma[SYMTYPE:("QST"), op_, {inv1_, inv2_}] := 
  Module[{t1, t2, tz1, tz2, tzop, op1, xa, xb},
   t1 = getOrbitalQN[SYMTYPE, inv1];
   t2 = getOrbitalQN[SYMTYPE, inv2];
   tz1 = subspacetz[t1];
   tz2 = subspacetz[t2];
   tzop = tz1-tz2;

   (* Spherical operators! *)
   If[tzop == 0, op1 = orbmomentz[op] ];
   If[tzop == 1, op1 = -1/Sqrt[2] orbmomentplus[op] ];
   If[tzop == -1, op1 = 1/Sqrt[2] orbmomentminus[op] ];

   xa = optransform[op1, inv1, inv2];
   xb = If[t1 == 0 && t2 == 0, 1, ClebschGordan[{t2,tz2},{1,tzop},{t1,tz1}]];
   xa/xb
];

(* Bug trap *)
ireducorbsigma[SYMTYPE_, args___] := 
  MyError["ireducorbsigma not implemented for " <> SYMTYPE <> 
  " args=" <> ToString[{args}]];

(* Table form for all irreducible matrix elements of a triplet operator sigma *)
ireducorbsigmaTable[op_] := Module[{t, i, cp},
  t = {};
  AppendTo[t, {nrorbcp}];
  For[i = 1, i <= nrorbcp, i++,
    cp = orbcoupledpairs[[i]];
    AppendTo[t, Flatten[cp]];
    t = Join[t,  ireducorbsigma[SYMTYPE, op, cp]];
  ];
  t (* return *)
];


(*                       --- SINGLET OPERATORS ---                   *)

(* 
 Call hierarchy: 
 maketable -> mtSingletOp -> mtOp -> singletopTable -> singletopMatrixSpeedy -> op2matrix
 NOTE: optransform is not used here!
*)

singletopMatrixSpeedy[op_, inv_] := Module[{mat, vecs, res},
  mat = op2matrix[op, inv];
  AppendTo[opdata, {inv, mat}]; (* opdata is global! *)
  vecs = diagvc[inv] [[2]];
  (* Unitarity transformation to the eigenbasis. This is essentially the same code as in optransform[]. *)
  res = Conjugate[vecs] . mat . Transpose[vecs]; (* conjugate bras, transpose kets *)
  res
];

generalopMatrixSpeedy[op_, inv1_, inv2_] := Module[{mat0, mat, vecs1, vecs2},
  mat0 = op2matrix[op, inv1, inv2];
  AppendTo[opdata, {inv1, inv2, mat0}]; (* opdata is global! *)
  vecs1 = diagvc[inv1] [[2]];
  vecs2 = diagvc[inv2] [[2]];
  mat = Conjugate[vecs1] . mat0 . Transpose[vecs2];
  If[Total[mat^2, 2] == 0, Return[0]];  (* Null matrix? *)
  mat
];

(* Makes a table with singlet operator irreducible matrix elements
   for every (inv) subspace. *)
singletopTable[op_] := Module[{t, i, inv, tmp},
  t = {{nrsub}};
  For[i = 1, i <= nrsub, i++,
    inv = subspaces[[i]];
    AppendTo[t, Flatten[{inv, inv}]];
    tmp = singletopMatrixSpeedy[op, inv];
    t = Join[t, N[tmp]]; (* NUMERICAL *)
  ];
  t
];

(* Strategy: test all combinations of spaces if they give non-zero matrix elements. *)
generalopTable[op_] := Module[{t, cnt, i, cp, mat},
  t = {};
  cnt = 0;
  For[i = 1, i <= nrp, i++, (* nrp = length of subspacepairs *)
    cp = subspacepairs[[i]];
    mat = generalopMatrixSpeedy[op, First[cp], Last[cp]];
    If[mat =!= 0,
      AppendTo[t, Flatten[cp]];
      t = Join[t, mat];
      cnt++;
    ];
  ];
  PrependTo[t, {cnt}];
  t (* Return *)
];

(* low-level functions called from mtSingletOp, mtGlobalOp & mtGeneralOp with different OPTABLEFNC *)
(* It first checks if the calculation of operator 'opname' was requested! *)
mtOp[opname_String, opinput_, prefix_, OPTABLEFNC_] :=  Module[{t, op},
  If[calcopq[opname],
    op = Expand[opinput];
    MyPrint[prefix, ": ", opname, " ", op];
    t = {};
    opfn = opfilename <> "." <> opname;
    opdata = {}; (* Global variable ! *)
    AppendTo[t, {prefix <> opname}];
    t = Join[t, OPTABLEFNC[ op ] ];
    MyPut[opdata, opfn, option["GENERATE_TEMPLATE"]];
    t,  
  (* else *) {}
  ]
];

mtSingletOp[opname_String, opinput_] := mtOp[opname, opinput, "s", singletopTable];
mtGeneralOp[opname_String, opinput_] := mtOp[opname, opinput, "p", generalopTable];
mtGlobalOp[opname_String, opinput_] := mtOp[opname, opinput, "g", singletopTable];
mtDoubletOp[opname_String, opinput_] := mtOp[opname, opinput, "d", ireducTable];
mtTripletOp[opname_String, opinput_] := mtOp[opname, opinput, "t", ireducsigmaTable];
mtOrbTripletOp[opname_String, opinput_] := mtOp[opname, opinput, "ot", ireducorbsigmaTable[op]];

(************* DIAGONALIZATION *************)

(* Find the ground state energy of the system. As a byproduct, all
diagonalisations will be performed and all eigenvalue/eigenvectors pairs
will be cached in memory for later use. *)
   
calcgsenergy[] := Module[{all, i, val, vec},
  MyPrint["calcgsenergy[]"];
  all = {}; (* List of all energy levels *)
  For[i = 1, i <= nrsub, i++,
    {val, vec} = diagvc @ subspaces[[i]];
    all = Join[all, val];
  ];
  all = Sort[all];
  first20 = Take[all, Min[Length[all], 20]];
  MyPrint["Lowest energies (absolute):", first20];
  GSenergy = First[all]; (* WARNING (side effect): GSenergy is a global variable! *)
  MyPrint["Lowest energies (GS shifted):", first20-GSenergy];
  MyPrint["Scale factor SCALE(Ninit):", SCALE[Ninit]];
  MyPrint["Lowest energies (shifted and scaled):", (first20-GSenergy)/SCALE[Ninit]];
];

(***********************************************************************)

ClearAll[thetaCh]; (* Bug honey-pot *)

(* ---- Tridiagonalisation approach: 
old - direct use of the recursion relations
sc - as above, but extended for superconducting hosts with arbitrary DOS 
     and even-frequency pairing function
sc2 - as above, but for fully general pairing function
cpp - tridiagonalisation performed in the C++ part of the code 
orth - recursion relations + orthogonality requirements
none - don't output the coeffcient table
manual - load the discretization tables from a file
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
If[TRI == "manual" || TRI == "manual_nambu",
  defaultprec = 30; (* Should be enough *)
  dothelanczos = loaddiscretizationtables;
];
If[option["GENERATE_TEMPLATE"],  
  dothelanczos = None;
  Nmax = 0;
];

(* Use arbitrary precision arithmetics *)
PREC = paramdefaultnum["prec", defaultprec];
MyPrint["PREC=", PREC];
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
MyPrint["mMAX=", mMAX];
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
  MyAssert[COEFCHANNELS == 1];
  imp1=Flatten[ImportTable["scdelta.dat"]];
  Print["scdelta=", Table[ scdelta[_][n] = imp1[[n+1]], {n,0,DISCNMAX} ]];
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
MyPrint["Diagonalisation."];

hookfile["hook_pre_lanczosinit"];
If[TRI != "manual" && TRI != "manual_nambu" && !option["GENERATE_TEMPLATE"],
  lanczosinit[];
];
hookfile["hook_post_lanczosinit"];

dothelanczos[];

(* Handle exceptions *)
If[TRI == "manual_nambu",
  If[SYMTYPE == "SPU1",
    sckappa[a_][n_] := 0;
    scdelta[a_][n_] := 1(zetaR[a][n]); (* on-site 'rung' term is the same as pairing *)
  ];
  If[SYMTYPE == "SPSU2",
    sckappa[a_][n_] := 0;
    loadtablescdelta[];
  ];
  If[SYMTYPE == "P" || SYMTYPE == "PP",
    sckappa[a_][n_] := 0;
    loadtablescdelta[];
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

(* IMPORTANT: As of 12.5.2010: use bcsgap=xx for spu1lr, and
bcsgap1=bcsgap2=xx for spu1/spsu2, etc. For a regression test, see
spu1lr4/. 

For f_0 site, included in the initial Hamiltonian, the
sign is positive by definition.

For symtype=SPSU2, the sign for isospin is inverted for f_1, etc.!

For single-channel problems, the sign issue is completely
irrelevant.

TO DO: This is all quite messy. One should clean this up, even if it
breaks backwards compatibility with param files. 
*)

gapdefined = False;

hookfile["hook_bcs"];

If[!gapdefined && paramexists["bcsgap"],
  gapdefined = True;
  bcssignch2 = If[paramdefaultbool["bcsinvertphase", False], -1, 1];
  scdelta[1][_] = paramnum["bcsgap"];
  
  scdelta[2][_] = bcssignch2 * paramnum["bcsgap"];
  
  scdelta[3][_] = paramnum["bcsgap"];
  sckappa[_][_] = 0;
];

(* Added 4.5.2016 - new possibility for SYMTYPE=SPSU2.
Equivalent to bcsgap1=bcsgap2, just saves some typing. *)
If[!gapdefined && paramexists["bcsgapspsu2"],
  gapdefined = True;
  scdelta[1][_] = paramnum["bcsgapspsu2"];
  
  scdelta[2][0] = paramnum["bcsgapspsu2"];
  scdelta[2][_] = -paramnum["bcsgapspsu2"];
  
  scdelta[3][_] = paramnum["bcsgapspsu2"];
  sckappa[_][_] = 0;
];

If[!gapdefined && paramexists["bcsgap1"] && paramexists["bcsgap2"] && !paramexists["bcsgap3"],
  gapdefined = True;
  scdelta[1][_] = paramnum["bcsgap1"];

  scdelta[2][0] = paramnum["bcsgap2"];
  scdelta[2][_] = -paramnum["bcsgap2"];

  sckappa[_][_] = 0;
];

(* See spsu2-3ch.m for sign conventions. *)
If[!gapdefined && paramexists["bcsgap1"] && paramexists["bcsgap2"] && paramexists["bcsgap3"],
  gapdefined = True;
  scdelta[1][_] = paramnum["bcsgap1"];

  scdelta[2][0] = paramnum["bcsgap2"];
  scdelta[2][_] = -paramnum["bcsgap2"];

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
  
Module[{a, precxi, preczeta},
  For[a = 1, a <= COEFCHANNELS, a++,
      MyPrintForm["Discretization (channel ``)", a];

      OUTPREC = 20;
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

(***** Add an arbitrary perturbation term to the Hamiltonian *****)

perturbhamiltonian[] := Module[{},
  If[PERTURB != "",
    If[option["READHAM"],
      MyError["If PERTURB is set, then the saved Hamiltonian matrices SHOULD NOT BE USED"];
    ];

    MyPrint["PERTURBATION: ", PERTURB];
    Hprime = Expand @ ToExpression[PERTURB];
    MyPrint["Hprime: ", Hprime];
    H = H + Hprime;
    MyPrint["New total Hamiltonian: ", H];
  ];
];

(***********************************************************************)

makeheader[] := {
  {"# Input file for NRG Ljubljana"},
  {"# symtype ", SYMTYPE},
  {"# Using sneg version ", $SnegVersion},
  {"#!9"}, (* Data file version: increased when incompatible changes are made *)
  {"# Number of channels, chain sites, subspaces: "},
  { CHANNELS, Nmax, nrsub }
};

(* Eigenenergies in each invariant subspace. IMPORTANT: We subtract ground
state energy and multiply by the scaling factor 1/SCALE[0], see PRB 41 9403,
Eq. (39) and Campo, Olivira PRB 72 104432. This factor is required to cancel
the N-dependence of xi[N] for large N. When more than one site of the Wilson
chain is included in the initial cluster (i.e. Ninit!=0), SCALE[0] is to be
replaced with SCALE[Ninit]. *)

makeenergies[] := Module[{i, inv, val, line},
  Flatten1 @ Table[
      inv = subspaces[[i]];
      val = First[diagvc[inv]];
      If[option["TEMPLATE"] || option["GENERATE_TEMPLATE"],
        line = {"DIAG ", hamfn[inv]},
      (* else *)
        line = (val-GSenergy)/SCALE[Ninit]; (* FC *)
      ];
      { inv, {dim[inv]}, line},
    {i, nrsub}
  ]
];

Flatten2[l_List] := Flatten[l, 2];
Flatten3[l_List] := Flatten[l, 3];

(* The last operator f on the Wilson chain. ch=0,1,2 *)
lastf[ch_] := If[Ninit == 0, f[ch], f[ch, Ninit]];


(* For U1, we need to distinguish between 
   spin-up and spin-down matrices, since we can't simply take 
   the difference of S_z of the subspaces in bra and ket. *)

   (* NOTE: in sneg, DO=0, UP=1. In NRG Ljubljana, and in indexing
   the xi coefficient tables, the order is reversed, first UP,
   then DO. For this reason, the j sum is reversed and we transform
   j -> 1-j in labeling the result. *)

makeireducf["U1"] := Module[{},
  MyPrint["makeireducf U1"];
  Flatten3 @ Table[{{{"f " <> ToString[i] <> " " <> ToString[1-j]}},
                   ireducTable[ lastf[i], j ]}, (* Recall: DO=0, UP=1 *)
                   {i, 0, CHANNELS-1}, {j, UP, DO, -1}]
];

(* We need to distinguish between two different 
operators wrt isospin symmetry. Recall: DO=0, UP=1. *)
makeireducf["SU2" | "DBLSU2"] := Module[{},
  MyPrint["makeireducf SI2 & DBLSU2"];
  Flatten3 @ Table[{{{"f " <> ToString[i] <> " " <> ToString[1-j]}},
                   ireducTable[ lastf[i], j ]},
                   {i, 0, CHANNELS-1}, {j, UP, DO, -1}]
];

(* For NONE & P, there are 4 types of operators *)
makeireducf["NONE" | "P" | "PP"] := Module[{},
  MyPrint["makeireducf NONE/P/PP"];
  Flatten3 @ Table[{{{"f " <> ToString[i] <> " " <> ToString[j]}},
                   ireducTable[ lastf[i], j ]},
                   {i, 0, CHANNELS-1}, {j, 0, 3}]
];

(* For QST & SPSU2T, there is a single coefficient table,
   since the three channels correspond to the orbital quantum
   number. *)
makeireducf["QST" | "SPSU2T"] := Module[{},
  MyPrint["makeireducf QST or SPSU2T"];
  Flatten1 @ {{ "f 0 0" }, ireducTable[ f[#1,#2,#3]& ], {"f 1 0"}, {"0"}, {"f 2 0"}, {"0"}}
];

(* Channel information is here carried by the Delta T_z quantum number !! *)

makeireducf["QSTZ" | "QSZTZ"] := Module[{},
  MyPrint["makeireducf QSTZ/QSZTZ"];
  Flatten1 @ {{ "f 0 0" }, ireducTable[ f[#1,#2,#3]& ], {"f 1 0"}, {"0"}, {"f 2 0"}, {"0"}}
];

(* p quantum number handling in QSC3 and SPSU2C3 *)
makeireducf["xxxxxQSC3"] := Module[{x, u, o},
  MyPrint["makeireducf QSC3 and SPSU2C3"];

  u = Exp[I 2 Pi/3];

  (* see c3_representations_NEW.nb for details *)
  x[s_] = {f[CR, 0, s], f[CR, 1, s], f[CR, 2, s]};
  o[0,s_] = x[s] . {1,1,1}/Sqrt[3];
  o[1,s_] = x[s] . {u^2,u,1}/Sqrt[3];
  o[2,s_] = x[s] . {u,u^2,1}/Sqrt[3];
  
  Flatten2 @ Table[{{{ "f " <> ToString[p] <> " 0" }}, ireducTable[ o[p,#]& ]}, {p, 0, 2}]
];

If[is["QJ"],
 fj[0, 1/2, -1/2] := -(Sqrt[2/3]*f[0, -1, 1]) + f[0, 0, 0]/Sqrt[3];
 fj[0, 1/2, 1/2] := -(f[0, 0, 1]/Sqrt[3]) + Sqrt[2/3]*f[0, 1, 0];
 fj[0, 3/2, -3/2] := f[0, -1, 0];
 fj[0, 3/2, -1/2] := f[0, -1, 1]/Sqrt[3] + Sqrt[2/3]*f[0, 0, 0];
 fj[0, 3/2, 1/2] := Sqrt[2/3]*f[0, 0, 1] + f[0, 1, 0]/Sqrt[3];
 fj[0, 3/2, 3/2] := f[0, 1, 1];
 fj[1, 1/2, -1/2] := -(Sqrt[2/3]*f[1, -1, 1]) + f[1, 0, 0]/Sqrt[3];
 fj[1, 1/2, 1/2] := -(f[1, 0, 1]/Sqrt[3]) + Sqrt[2/3]*f[1, 1, 0];
 fj[1, 3/2, -3/2] := f[1, -1, 0];
 fj[1, 3/2, -1/2] := f[1, -1, 1]/Sqrt[3] + Sqrt[2/3]*f[1, 0, 0];
 fj[1, 3/2, 1/2] := Sqrt[2/3]*f[1, 0, 1] + f[1, 1, 0]/Sqrt[3];
 fj[1, 3/2, 3/2] := f[1, 1, 1];

 coupledpairsQJ[1/2] = Select[subspacepairs, coupledQ[SYMTYPE, #, 1/2]&];
 coupledpairsQJ[3/2] = Select[subspacepairs, coupledQ[SYMTYPE, #, 3/2]&];

 ireducTableQJ[op_, j_] := Module[{t, nrcp, coupledpairs, cp, i},
  MyPrint["ireducTableQJ: ", op, j];
  coupledpairs = coupledpairsQJ[j];
  nrcp = Length[coupledpairs];
  t = {{nrcp}};
  For[i = 1, i <= nrcp, i++,
    cp = coupledpairs[[i]];
    AppendTo[t, Flatten[cp]];
    t = Join[t, ireducMatrixSpeedy[SYMTYPE, op, cp, j] ];
  ];
  t (* Return *)
 ];
];

(* Channel QN masquerading as j. ch=0 => j=1/2, ch=1 => j=3/2 *)
makeireducf["QJ"] := Module[{},
 MyPrint["makeireducf QJ"];
 Flatten1 @ {{ "f 0 0" }, ireducTableQJ[ fj[#1, 1/2, #2]&, 1/2 ],
             { "f 1 0" }, ireducTableQJ[ fj[#1, 3/2, #2]&, 3/2 ],
             { "f 2 0" }, {"0"}}
];

(* Generic: ireducible matrix elements <||f||> *)
makeireducf[_] := Module[{},
  MyPrint["makeireducf GENERAL"];
  Flatten2 @ Table[{{{ "f " <> ToString[i] <> " 0" }}, 
                  ireducTable[ lastf[i] ]}, {i, 0, CHANNELS-1}]
];

(* Discretization parameter tables *)
makedisctables[] := Module[{t, a},
   t = {{"z"}};
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, xitable[a]];
   ];
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, zetatable[a]];
   ];
   t
];

makeRdisctables[] := Module[{t, a},
   t = {{"X"}};
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, xiRtable[a]];
   ];
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, zetaRtable[a]];
   ];
   t
];

makescdisctables[] := Module[{t, a},
   t = {{"Z"}};
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, scdeltatable[a]];
   ];
   For[a = 1, a <= COEFCHANNELS, a++,
     AppendTo[t, {DISCNMAX}];
     t = Join[t, sckappatable[a]];
   ];
   t
];

maketritables[] := Module[{t, a},   
  t = {{"T"}};
  (*  Coefficients Epsilon_j^z *)
  For[a = 1, a <= COEFCHANNELS, a++,
    AppendTo[t, {mMAX}];
    t = Join[t, eptable[a]];
  ];
  For[a = 1, a <= COEFCHANNELS, a++,
    AppendTo[t, {mMAX}];
    t = Join[t, emtable[a]];
  ];
  (* Integrals of hybridisation *)
  For[a = 1, a <= COEFCHANNELS, a++,
    AppendTo[t, {mMAX}];
    t = Join[t, u0ptable[a]];
  ];
  For[a = 1, a <= COEFCHANNELS, a++,
    AppendTo[t, {mMAX}];
    t = Join[t, u0mtable[a]];
  ];
  t
];   
   
(* maketable[] builds a table which is to be saved as a file named 'data' 
and used as input to the NRG iteration routines in the C++ part of the
software package. The beginning block of this function is the right place 
to perform various extra tweaks, log parameters and other debugging info, 
check if parameters make any sense, etc. *) 

maketable[]:=Module[{t},
  MyPrint["\n\nmaketable[]\n", makeheader[]];
  timestart["maketable"];                   

  addexnames[]; (* To be on the safe side... *)
  perturbhamiltonian[];
  inittheta0ch[];
  If[!option["GENERATE_TEMPLATE"], checkdefinitions[]];

  (* Perform all diagonalisations *)  
  calcgsenergy[]; 

  t = Join[makeheader[],
           {{"# SCALE ", SCALE[Ninit]}},
           {{"# Energies (GS energy subtracted, multiplied by 1/SCALE):"}},
           makeenergies[],
           {{"# Irreducible matrix elements for Wilson chains:"}},
           makeireducf[SYMTYPE],
           {{"# GS energy in absolute units:"}},
           {{"e"}, {GSenergy}},
           {{"# Irreducible matrix elements for other operators:"}}
          ];

  (* Operator definitions *)
  tops = loadmodule["operators.m"];
  t = Join[t, tops];

  tops = loadmodule["customoperators.m", False];
  If[ tops =!= $Failed,
    t = Join[t, tops];
  ];

  (* Model specific operators *)  
  tops = loadmodule["modeloperators.m", False];
  If[ tops =!= $Failed,
    t = Join[t, tops];
  ];
  
  If[TRI != "none" && !option["GENERATE_TEMPLATE"],
    t = Join[t,
             {{"# Discretization tables:"}},
             makedisctables[]
    ];
    If[isSC[],
      t = Join[t, makescdisctables[]];
    ];
    If[RUNGS,  
      t = Join[t, makeRdisctables[]];
    ];    
  ];
  
  If[TRI == "cpp" && !option["GENERATE_TEMPLATE"],
    t = Join[t, maketritables[]];
  ];

  (* Enforce linefeed on the last line. *)  
  AppendTo[t, {}]; 

  (* Some C++ compilers/libraries generate code which has trouble parsing
  very small (non-representable) floating point numbers. In this case one
  can use option EPSCLIP to truncate all small floating point numbers to
  zero. *)
  If[option["EPSCLIP"],
    EPSREPRESENTABLE = 10^-300; (* Actually smallest double is approx. 10^-324. *)
    t = t /. { x_Real /; Abs[x] < EPSREPRESENTABLE -> 0};
    t = t /. { z_Complex /; Abs[z] < EPSREPRESENTABLE -> 0};  
    t = t /. { 0. -> 0, Complex[0.,0.] -> 0 };
  ];

  If[option["CHOP"],
    EPSCHOP = 10^-12;
    t = t /. { x_Real /; Abs[x] < EPSCHOP -> 0};
    t = t /. { z_Complex /; Abs[z] < EPSCHOP -> 0};
    t = t /. { z_Complex /; Abs[Im[z]] < EPSCHOP :> Re[z] };
    t = t /. { 0. -> 0, Complex[0.,0.] -> 0 };
  ];

  MyPrint["-- maketable[] done --"];

  timeadd["maketable"];
  timereport[];

  t   (* RETURN *)
];

(* Actually generate and write datafile *)
makedata[filename_]:=Module[{suffix, fn, tabelca},
  suffix = If[option["TEMPLATE"] || option["GENERATE_TEMPLATE"], ".in", ""];
  fn = filename <> suffix;
  tabelca = maketable[];
  iscomplex = !FreeQ[tabelca, Complex[_,_]] || option["COMPLEX"];
  If[iscomplex,
    tabelca = tabelca /. Complex[x_,y_] :> 
      "(" <> ToString[CForm[x]] <> "," <> ToString[CForm[y]] <> ")";
    tabelca[[3]] = "# COMPLEX";
  ];
  MyPrint[Export[fn, tabelca, "Table"]];

  (* Various troubleshooting aids *)
  MyPrint["gammaPol=", gammaPol /. params];
];

MyPrint["--EOF--"];

"initial.m loaded"
