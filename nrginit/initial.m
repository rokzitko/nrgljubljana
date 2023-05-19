(*
   NRG Ljubljana -- Numerical renormalization group code

   Copyright (C) 2005-2023 Rok Zitko

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

VERSION = "2023.05";

(* Logging of Mathematica output: this is useful for bug hunting *)
If[!ValueQ[mmalog],
  mmalog = OpenWrite["mmalog"];
  AppendTo[$Output, mmalog];
  AppendTo[$Messages, mmalog]; (* Important! *)
];

SetOptions[$Output, PageWidth -> 240];
SetOptions[$Messages, PageWidth -> 240];

DEBUG=1;
Get["misc.m", Path->PACKAGEPATH];

Print2["NRG Ljubljana ", VERSION, " (c) Rok Zitko, rok.zitko@ijs.si, 2005-2023"];
Print2["Mathematica version: ", $Version];
Print2["sneg version: ", $SnegVersion];

Off[InterpolatingFunction::dmvali];
Off[InterpolatingFunction::dmval];

(* No recursion limit *)
$RecursionLimit = Infinity;

(* List of directories in which we look for modules. *)
modulespath = {".", ".."};
If[ValueQ[PACKAGEPATH], modulespath = Join[modulespath, PACKAGEPATH]];

(* Parse parameter file, if not already performed by nrginit *)
loadmodule["initialparse.m"];
If[PARSED =!= True,
  MyPrint["Parsing parameters"];
  parse["param"];
];

loadmodule["sym.m"];

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

(* Define extraPARAM=VALUE for all PARAM=VALUE lines in [extra]
   block of the input file (for backward compatibility). *)
If[listdata["extra"] =!= {},
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
during the NRG iteration. It is sorted alphabetically to ensure consistency
of the generated data file irrespective of the ordering in the input file. *)
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
BANDSPIN = 1/2;

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
Conjugate[coefkappa[i__]] ^= coefkappa[i];

(* ! *)
(* Conjugate[coefdelta[i__]] ^= coefdelta[i]; *)

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

Get["hamiltonian.m", Path->PACKAGEPATH];

(************************ Parameter handling *****************************)

(* All PARAM=VALUE lines in [extra] block of the input file get transformed
   into PARAM->VALUE rules in 'params'. Since rules are applied in the order
   of their appearance in the list, this implies that the preexisting rules
   take precendence over the automatically appended ones. *)
If[ValueQ[listdata["extra"]],
  params2 = Map[ToExpression[First[#]] -> parsevalue @ Last[#] &,
                listdata["extra"]];
  params = Join[params, params2];
];
MyVPrint[2,"params=", params];

(* Important: all left sides of the table 'params' should be defined as sneg
parameters. Currently they are all real constants! *)
paramnames = params[[All, 1]];
paramnames = Select[paramnames, AtomQ]; (* Drop parametrized rules *)
paramnames = Select[paramnames, (# =!= "Gamma")&]; (* Gamma is Protected *)
snegrealconstants @@ paramnames;

If[NRDOTS == -1, MyError["NRDOTS = -1"]]; (* Bug trap *)
If[CHANNELS == -1, MyError["CHANNELS = -1"]]; (* Bug trap *)

MyPrint["NRDOTS:", NRDOTS, " CHANNELS:", CHANNELS, " COEFCHANNELS:", COEFCHANNELS];

basopsch = Table[f[n], {n, 0, CHANNELS-1}];
Do[basopsch = Join[basopsch, Table[f[n, i], {n, 0, CHANNELS-1}]], {i, 1, Ninit}];
If[Length[basopsdot]==0 && Length[impuritybasis]==0,
  basopsdot = Take[{d[], a[], b[], e[], g[]}, NRDOTS];
];
If[Length[impuritybasis]>0, basopsdot = impuritybasis];
basisops = basops = Sort[Join[basopsch, basopsdot]];

MyPrint["basis:", basisops];
If[lrchain =!= {}, MyPrint["lrchain:", lrchain]];
If[lrextrarule =!= {}, MyPrint["lrextrarule:", lrextrarule]];

(* Check for consistency between the number of operators and the numbers of
   dots and channels *)
NROPS = Length[basisops];
If[NROPS != NRDOTS + CHANNELS (1+Ninit),
  MyError["Number of operators does not match what we were expecting."];
];

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
If[paramdefaultbool["checkHc", True] && Not[option["GENERATE_TEMPLATE"]],
  Module[{Hcheck, Hdiff},
    Hcheck = conj[H];
    Hdiff = Simplify[Expand[H-Hcheck]];
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

Get["basis.m", Path->PACKAGEPATH];

(********************** EXACT DIAGONALIZATION ************************)

(* At this point the basis has been generated (or it has been read
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

(* List of directories in which we look for dump files of basis vectors and
Hamiltonian matrices. In addition to the current directory, we also look
in the parent; this is useful for parameter sweeps. *)
dumppath = {".", ".."};

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
(* WARNING: <i| op |j>; q1,s1 correspond to j, q2,s2 correspond to i *)

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

coupledQ["DBLQSZ", {{q11_, q21_, ssz1_, j1___}, {q12_, q22_, ssz2_, j2___}}] :=
  If[    (* ((q11 == q12+1 && q21 == q22) ~ Xor ~ (q11 == q12 && q21 == q22+1))  && *)
  Abs[ssz1-ssz2] == 1,
  True, False, MyError["oops coupledQ"]];

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

(* Bug trap *)
coupledQ[s_String, a___] := MyError["coupledQ not defined for ", {s,a}];

(* ======================= *)

MyVPrint[3, "subspaces=", subspaces];

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

getSpinQN["DBLQSZ", inv_] := inv[[3]];
getSpinQN["QS" | "QSLR" | "QSC3" | "QSZ" | "QSZLR" | "ISO" | "ISOLR" |
          "ISO2" | "ISO2LR" | "ISOSZ" | "ISOSZLR" | "QST" | "QSTZ" | "QSZTZ", inv_] := inv[[2]];
getSpinQN["SPSU2" | "SPSU2LR" | "SPSU2T" | "SPU1" | "SPU1LR" | "SPSU2C3", inv_] := inv[[1]];

getOrbitalQN["QST" | "QSTZ" | "QSZTZ", inv_] := inv[[3]];
getOrbitalQN["SPSU2T", inv_] := inv[[2]];

getJQN["QJ", inv_] := inv[[2]];

ireducMatrixSpeedy[symtype:("QSZ" | "QSZLR" | "SPU1" | "SPU1LR" | "DBLQSZ"),
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
Module[{t, cp, i, mat, opfnsub},
  t = {{nrcp}};
  For[i = 1, i <= nrcp, i++,
    (* coupledpairs is a list of subspace pairs that are coupled
       by doublet operators [that increase charge, when charge conservation
       is explicitly taken into account], i.e. creation operators! *)
    cp = coupledpairs[[i]];
    mat = Expand @ ireducMatrixSpeedy[SYMTYPE, op, cp, optional];
    AppendTo[t, Flatten[cp]];
    AppendTo[opdata, {cp, mat}];
    If[!option["GENERATE_TEMPLATE"] || opfn === "",
      t = Join[t, mat],
    (* else *)
      opfnsub = opfn <> "_" <> Invar2String[cp[[1]]] <> "_" <> Invar2String[cp[[2]]];
      t = Join[t, {opfnsub}];
      Put[mat, opfnsub];
    ];
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

ireducSPIN[SYMTYPE:("SPSU2" | "QS" | "QSLR" | "QSC3" | "QST" | "QSTZ" | "SPSU2T" | "SPSU2C3"),
 SPIN_, {inv1_, inv2_}] :=
  Module[{ss1, ss2, sz1, sz2, szop, op1, xa, xb},
   ss1 = getSpinQN[SYMTYPE, inv1];
   ss2 = getSpinQN[SYMTYPE, inv2];
   sz1 = subspacesz[ss1];
   sz2 = subspacesz[ss2];
   szop = sz1-sz2;

   (* Spherical operators! *)
   If[szop == 0, op1 = spinketbraZ[SPIN] ];
   If[szop == 1, op1 = -1/Sqrt[2] spinketbraP[SPIN] ];
   If[szop == -1, op1 = 1/Sqrt[2] spinketbraM[SPIN] ];

   xa = optransform[op1, inv1, inv2];
   xb = If[ss1 == 1 && ss2 == 1, 1,
     ClebschGordan[{SS2S[ss2],sz2},{1,szop},{SS2S[ss1],sz1}]];
   xa/xb
];

ireducSPINList[SYMTYPE:("SPSU2" | "QS" | "QSLR" | "QSC3" | "QST" | "QSTZ" | "SPSU2T" | "SPSU2C3"),
 {ndx_}, {inv1_, inv2_}] :=
  Module[{ss1, ss2, sz1, sz2, szop, op1, xa, xb, SPIN, id, skl},
   ss1 = getSpinQN[SYMTYPE, inv1];
   ss2 = getSpinQN[SYMTYPE, inv2];
   sz1 = subspacesz[ss1];
   sz2 = subspacesz[ss2];
   szop = sz1-sz2;

   Which[
    ndx == 1, SPIN = SPIN1; skl = {1,0}; id = spinketbraI[SPIN2, {0,1}],
    ndx == 2, SPIN = SPIN2; skl = {0,1}; id = spinketbraI[SPIN1, {1,0}]
   ];

   (* Spherical operators! *)
   If[szop == 0, op1 = spinketbraZ[SPIN, skl] ];
   If[szop == 1, op1 = -1/Sqrt[2] spinketbraP[SPIN, skl] ];
   If[szop == -1, op1 = 1/Sqrt[2] spinketbraM[SPIN, skl] ];

   op1 = op1 ~ nc ~ id;

   xa = optransform[op1, inv1, inv2];
   xb = If[ss1 == 1 && ss2 == 1, 1,
     ClebschGordan[{SS2S[ss2],sz2},{1,szop},{SS2S[ss1],sz1}]];
   xa/xb
];

(* Bug trap *)
ireducsigma[SYMTYPE_, args___] :=
  MyError["ireducsigma not implemented for " <> SYMTYPE <>
  " args=" <> ToString[{args}]];

ireducSPIN[SYMTYPE_, args___] :=
  MyError["ireducSPIN not implemented for " <> SYMTYPE <>
  " args=" <> ToString[{args}]];

ireducSPINList[SYMTYPE_, args___] :=
  MyError["ireducSPINList not implemented for " <> SYMTYPE <>
  " args=" <> ToString[{args}]];

(* Table form for all irreducible matrix elements of a triplet operator sigma *)
ireducsigmaTable[op_] := Module[{t, i, cp, mat},
  t = {};
  AppendTo[t, {nrspincp}];
  For[i = 1, i <= nrspincp, i++,
    cp = spincoupledpairs[[i]];
    Which[
      NumberQ[op], mat = ireducSPIN[SYMTYPE, op, cp],
      Head[op] === List, mat = ireducSPINList[SYMTYPE, op, cp],
      True, mat = ireducsigma[SYMTYPE, op, cp]
    ];
    AppendTo[t, Flatten[cp]];
    t = Join[t, mat];
    AppendTo[opdata, {cp, mat}];
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
ireducorbsigmaTable[op_] := Module[{t, i, cp, mat},
  t = {};
  AppendTo[t, {nrorbcp}];
  For[i = 1, i <= nrorbcp, i++,
    cp = orbcoupledpairs[[i]];
    mat = ireducorbsigma[SYMTYPE, op, cp];
    AppendTo[t, Flatten[cp]];
    t = Join[t, mat];
    AppendTo[opdata, {cp, mat}];
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
singletopTable[op_] := Module[{t, i, inv, mat},
  t = {{nrsub}};
  For[i = 1, i <= nrsub, i++,
    inv = subspaces[[i]];
    AppendTo[t, Flatten[{inv, inv}]];
    mat = Simplify @ singletopMatrixSpeedy[op, inv]; (* simplify! *)
    If[!option["GENERATE_TEMPLATE_ALL"] || opfn === "",
      t = Join[t, N[mat]],  (* NUMERICAL *)
    (* else *)
      opfnsub = opfn <> "_" <> Invar2String[inv] <> "_" <> Invar2String[inv];
      t = Join[t, {opfnsub}];
      Put[mat, opfnsub];
    ];
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
mtDoubletOp[opname_String, opinput_, opt_] := mtOp[opname, opinput, "d", ireducTable[#1,opt]& ];
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

Get["wilson.m", Path->PACKAGEPATH];

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

If[paramdefaultbool["data_has_rescaled_energies", True], energiesscale = SCALE[Ninit], energiesscale = 1];

makeenergies[] := Module[{i, inv, val, line},
  Flatten1 @ Table[
      inv = subspaces[[i]];
      val = First[diagvc[inv]];
      If[option["TEMPLATE"] || option["GENERATE_TEMPLATE"],
        line = {"DIAG ", hamfn[inv]},
      (* else *)
        line = (val-GSenergy)/energiesscale
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

(* New interface *)
makewilsonchain["QS"] := Module[{mat},
  t = {{"W"}};
  AppendTo[t, {DISCNMAX, COEFCHANNELS, COEFCHANNELS}];
  (* epsilon *)
  For[n = 0, n <= DISCNMAX, n++,
    mateps = DiagonalMatrix[Table[bandrescale zeta[a][n], {a, 1, COEFCHANNELS}]];
    t = Join[t, N[mateps, OUTPREC]];
  ];
  (* t *)
  For[n = 0, n <= DISCNMAX, n++,
    matt = DiagonalMatrix[Table[bandrescale xi[a][n], {a, 1, COEFCHANNELS}]];
    t = Join[t, N[matt, OUTPREC]];
  ];
  t
];

(* Old interface *)
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
  timestart["maketable"];

  perturbhamiltonian[];
  inittheta0ch[];
  If[!option["GENERATE_TEMPLATE"], checkdefinitions[]];

  (* Perform all diagonalisations *)
  calcgsenergy[];

  opfn=""; opdata={}; (* Prior to makeireducf[] call to silence errors *)

  t = Join[makeheader[],
           {{"# SCALE ", energiesscale}},
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
    If[WILSONCHAIN == "legacy",
      t = Join[t,
               {{"# Discretization tables (legacy form):"}},
               makedisctables[]
      ];
      If[isSC[],
        t = Join[t, makescdisctables[]];
      ];
      If[RUNGS,
        t = Join[t, makeRdisctables[]];
      ];
    ];
    If[WILSONCHAIN == "matrix",
      t = Join[t,
               {{"# Discretization tables (matrix form):"}},
               makewilsonchain[SYMTYPE]
      ];
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
];
