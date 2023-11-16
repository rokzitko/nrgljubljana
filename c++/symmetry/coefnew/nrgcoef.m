(*
Derivation of NRG coefficient tables
Rok Zitko, rok.zitko@ijs.si, 2005-2023

CHANGE LOG
15.1.2008 - generalisation to all symmetry types
18.1.2008 - matrixniNEW[] function
21.1.2008 - QN.dat generation
22.1.2008 - tos2[], tos3[], some code cleanup
28.4.2008 - fix for "U1" symtype
13.11.2008 - recalcopONELOOP[]
14.11.2008 - debugging instrumentalization
12.1.2009 - dorecalcfLOOP["NONE"]
18.8.2009 - new defstr[]
11.8.2015 - more diagnostics
12.8.2015 - mynewbk[a,b,c] instead of vev[nc[a,b,c]]
          - HoldPattern added
30.10.2015 - SUMCH[] functions
9.11.2015 - optimizations
19.1.2016 - ONELOOP[] removed
2.2.2016 - formatting
26.4.2016 - predznak[] bug fix
14.5.2023 - cleanup
*)

SetOptions[$Output, PageWidth -> 1000];
SetOptions[$Messages, PageWidth -> 1000];
SetOptions[$Output, FormatType -> InputForm];
Setoptions[$Messages, FormatType -> InputForm];

Off[General::"spell", ClebschGordan::"tri", ClebschGordan::"phy", Simplify::"infd"];

(* fnprep[] must be called after PREFIX had been defined. *)
fnprep[] := Module[{},
  FNSINGLET = PREFIX <> "-singlet";
  FNDOUBLET = PREFIX <> "-doublet";
  FNTRIPLET = PREFIX <> "-triplet";

  FNRECALCF = PREFIX;

  FNSPINUP =   FNRECALCF <> "-spinup";
  FNSPINDOWN = FNRECALCF <> "-spindown";

  FNSPINUPISOUP =     FNRECALCF <> "-spinup-isoup";
  FNSPINDOWNISOUP =   FNRECALCF <> "-spindown-isoup";
  FNSPINUPISODOWN =   FNRECALCF <> "-spinup-isodown";
  FNSPINDOWNISODOWN = FNRECALCF <> "-spindown-isodown";
];

(* Ensure that the last line is termined with a linefeed. *)
MyExport[fn_, l_List, st_:"Lines"] := Export[fn, Append[l,""], st];

(* HINT: insert MyPrint lines at suitable places for debugging assistance. *)
If[!ValueQ[debug], debug = 0];
MyPrint[msg__] := If[debug >= 1, Print[msg]];
MyPrint[i_Integer, msg__] := If[debug >= i, Print[msg]];

(* Generic support for up to five channels, a-e! 
f denotes the states on the previous site of the Wilson chain. *)

snegfermionoperators[a, b, c, d, e, f];
allops = {a[], b[], c[], d[], e[]};

expr2list[HoldPattern[x : Plus[l__]]] := {l};
expr2list[x_] = {x};

ClearAll[newbk1, newbk2];
newbk1[a___, z_?NumberQ b_, c___] := z newbk1[a, b, c];
newbk1[a_, b_, c_] := newbk2[nc[a, b, c]];
newbk2[z_?NumberQ] := z;
newbk2[z_. a_nc] := z vevwicknew[a];
newbk2[z_. x:o_?operatorQ[___]] := 0; (* Can happen... *)

(* <b|op|k> calculation using a new implementation of VEV with Wick's
   theorem *)
newbk[b_, op_, k_] := Module[{bl, opl, kl, r},
 (* Approach: expand the sums and calculate VEV term by term *)
 bl = expr2list @ Expand[b];
 opl = expr2list @ Expand[op];
 kl = expr2list @ Expand[k];
 r = Outer[newbk1, bl, opl, kl];
 Total[Flatten[r]]
];

(* Optimized <b|op|k> *)
mynewbk[x_, op_, y_] := Module[{},
  ordering[a] = NONE;
  ordering[b] = NONE;
  ordering[c] = NONE;
  ordering[d] = NONE;
  ordering[e] = NONE;
  ordering[f] = NONE;
  res = newbk[x, op, y];
  ordering[a] = EMPTY;
  ordering[b] = EMPTY;
  ordering[c] = EMPTY;
  ordering[d] = EMPTY;
  ordering[e] = EMPTY;
  ordering[f] = EMPTY;
  res
];

(* Call basisprep[] after the basis had been generated using one of
the ???basis[] function from package "sneg". mult[] takes as argument
the quantum numbers which designate the invariant subspace and
returns the degeneracy of each state in that subspace. *)

basisprep[] := Module[{},
  MyPrint["basisprep[]"];
  MyPrint[basis // MatrixForm];
  MyPrint[nrstates = Plus @@ ((mult[ #[[1]] ] Length[#[[2]]])& /@ basis)];
];

(* Generally usefull utility functions *)

toc[koef_] := ToString[CForm[koef]];
tos[koef_] := ToString[koef];

toi[koef_] := ToString[koef, InputForm]; (* Useful for debugging output! *)

tos2[x_, y_] := "(" <> ToString[x] <> ", " <> ToString[y] <> ")";
tos3[x_, y_, z_] := "(" <> ToString[x] <> ", " <> ToString[y] <> 
                       ", " <> ToString[z] <> ")";
tos4[x_, y_, w_, z_] := "(" <> ToString[x] <> ", " <> ToString[y] <>
                       ", " <> ToString[w] <> ", " <> ToString[z] <> ")";

toc2[x_, y_] := "(" <> toc[x] <> ", " <> toc[y] <> ")";
toc3[x_, y_, z_] := "(" <> toc[x] <> ", " <> toc[y] <> ", " <> toc[z] <> ")";
toc4[x_, y_, w_, z_] := "(" <> toc[x] <> ", " <> toc[y] <> ", " <> toc[w] <> ", " <> toc[z] <> ")";

(* Simplification of rational expressions. *)
ratsimpl[expr_, rules_] :=
  FullSimplify[expr /. (a_. Sqrt[b_. x_] + c_. Sqrt[d_. x_]) / Sqrt[e_. x_] /;
    Simplify[x > 0, rules] :> (a Sqrt[b]+c Sqrt[d])/Sqrt[e], rules];

spinofop[izn_, szn_] := 
  If[izn == CR && szn == UP, 1/2, 0] + 
  If[izn == CR && szn == DO, -1/2, 0] +
  If[izn == AN && szn == UP, -1/2, 0] +
  If[izn == AN && szn == DO, 1/2, 0];

isoof[CR] = 1/2;
isoof[AN] = -1/2;

CG[x__] := ClebschGordan[x];

skpdtsimpl[x_] := x;

skpdt[x_, y_] := Module[{sx=skpdtsimpl[x], sy=skpdtsimpl[y], res},
  MyPrint[3, x, " -> ", sx];
  MyPrint[3, y, " -> ", sy];
  res=If[sx === sy, 1, 0];
  MyPrint[3, "res=", res];
  res
];

(* Convert a coefficient into a form suitable for inclusion in C++ code. *)
koefstr[koef_] := Module[{str = toc[koef]},
  str = StringReplace[str, "pow("->"Power("];
  (* NOTE: avoid upper case S, which is used to denote spin. *)
  str = StringReplace[str, "Sign(" -> "sign("];
  str = StringReplace[str, "Sqrt("->"sqrt("]
];

(* String for a C++ #define line. *)
defstr[key_, val_] := "const unsigned int " <> tos[key] <>
  " =  " <> tos[val] <> ";";

(* Display the table, store it in a file.
MyExport takes care of the LF in the last line. *)

store[l_, filename_] := Module[{},
  Print[MyExport[filename, l, "Lines"]];
];

(* fixspin[] takes an operator for maximal sz for a given spin s and 
returns an operator with lower sz. To be used in fixstate[] and similar
functions. *)

fixspin[s_,sz_,op_] := Nest[spindown, op, s-sz];

fixisospin[i_, iz_, op_] := Nest[isospindown, op, i-iz];

(* newstates[]: mapped on the basis. qn: quantum numbers, states: the list of states *)

(* Loop over all states in an invariant subspace. Called from donew[]
generic routine. newstates1[] must be defined in the symmetry-specific
script file. newstates1[], in turn, would typically call newstates2[]. *)

newstates[{qn_, states_}] := Scan[newstates1[qn, #]&, states];

(* Low level: form the suitable combination of states. 
Called from newstates1[] symmetry-specific routine. 
'qn' are the preserved quantum numbers, 'd' are the additional
quantum numbers which drop out of the problem due to the existance
of the symmetry, 'state' is the state of the N-th site electrons
that is being considered.
getexpr1[] is a symmetry-specific routine which returns the
combination of f[N] states and states from the previous iteration
with well defined quantum numbers. *)

newstates2[qn_, d_, state_] := Module[{expr},
  expr = getexpr1[qn, d, state];
  If[SIMPLIFYNEW, expr = simplifynew @ expr];
  new[nrcomb] = expr;

  string = "In[" <> tos[nrcomb] <> "] = " <> Invar[d];
  Print[string];
  AppendTo[code, string];

  string = "QN[" <> tos[nrcomb] <> "] = " <> InvarQN[qn];
  Print[string];
  AppendTo[codeQN, string];

  string = "[" <> tos[nrcomb] <> "] qn=" <> InvarQN[qn] <> " d=" <> Invar[d];
  string = string <> " state=" <> tos[state];
  AppendTo[details, string];
  string = tos[expr];
  AppendTo[details, string];

  (* TRICK: we store qn, d and state for later use. *)
  newqn[nrcomb] = qn;
  newd[nrcomb] = d;
  newstate[nrcomb] = state;

  nrcomb++;

  Print["expr=", expr];
];

(* Do the calculation: determine all P::combs possible combinations
of states in the new NRG iteration. *)
donew[] := Module[{},
  MyPrint["donew[]"];
  nrcomb = 1;
  code = {};
  codeQN = {};
  details = {};
  Map[newstates, basis];
  nrcomb--;
  If[nrcomb != nrstates,
    Print["Error in donew[]!"];
    Exit[1];
  ];
  store[code, PREFIX <> "-In2.dat"];
  store[codeQN, PREFIX <> "-QN.dat"];
  store[details, PREFIX <> "-details.dat"];

  (* Store for reload using loadnew[]. NEW, 4.3.2016 *)
  Definition[new] >> "new.dat";
  Definition[newqn] >> "newqn.dat";
  Definition[newd] >> "newd.dat";
  Definition[newstate] >> "newstate.dat";
  Definition[nrcomb] >> "nrcomb.dat";
];

(* NEW, 4.3.2016 *)
loadnew[] := Module[{},
 MyPrint["Loading precomputed coefficients"];
 Get["new.dat"];
 Get["newqn.dat"];
 Get["newd.dat"];
 Get["newstate.dat"];
 Get["nrcomb.dat"];
 MyPrint["Loaded."];
];

(* Factor out numerical coefficients. *)
factorcoef[expr_] := expr //. HoldPattern[ nc[a___, coef[x_], b___] ] :> simpl2[x] nc[a, b];

(* pbconj[] allow for overloading the conjugation operation: this is
convenient, for example, for isospin symmetry due to the (-1)^N factors,
which need to be simplified after the conjugation. *) 

pbconj[expr_] := conj[expr];

(* Prepare a bra state: perform a Hermitian conjugation! *)
(* 9.11.2015 - caching *)
(* 9.11.2015 - simpl2 *)

prepbra[expr_] := prepbra[expr] = simpl2[ expr //. dirpdt[{a___}, state_] :> nc[ old[a], pbconj[state] ] ];

prepket[expr_] := prepket[expr] = simpl2[ expr //. dirpdt[{a___}, state_] :> nc[state,  old[a] ] ];

prepall[i1_, i2_, ops_] := Module[{st1, st2},
  st1 = prepbra @ new[i1];
  st2 = prepket @ new[i2];
  factorcoef @ nc[st1, ops, st2]
];

(* Returns -1 for an odd number of fermions operators in 'state'
and 1 for an even number of fermion operators. Called from prepbrasign[]. *)

(* TODO: this only works for cases where the charge is conserved!!! FIX ME!!! *)
countXX[a_ /; FreeQ[a, XX]] = 0;
countXX[a_?isnumericQ b_] := countXX[b];
countXX[XX] := 1;
countXX[HoldPattern[l : nc[XX ..]]] := Length[l];
predznakOLD[state_] := (-1)^countXX[state /. (op_[CR, ___] /; operatorQ[op] -> XX)];

ClearAll[opparity];
opparity[z_?NumericQ] := 0;
opparity[z_?NumericQ a_] := opparity[a];
opparity[o_?operatorQ[a___]] := 1;
opparity[HoldPattern[l_nc]] := Mod[Length[l], 2];
opparity[(-1)^NN] := 0;
opparity[(-1)^(NN+1)] := 0;
opparity[(-1)^NN a__] := opparity[a];
AllSame[{1 ..}] := 1;
AllSame[{0 ..}] := 0;
opparity[HoldPattern[Plus[a__]]] := AllSame[Map[opparity, {a}]];
opsign[state_] := (-1)^opparity[state];
predznak[state_] := opsign[state];

predznak[xyz[a__, state_]] := predznak[state];
conj[xyz[a__]] := xyz[a];
nc[xyz[a__], xyz[b__]] := If[{a} === {b}, 1, 0];
mynewbkX[a_?NumberQ, b_?NumberQ, c_?NumberQ] := a b c;

(* Prepare a bra state, taking into account necessary permutations to 
anticommute out the operator being considered. Used in recalcop[] for
non-integer spin of the operator. *)
(* 9.11.2015 - caching *)

prepbrasign[expr_] := prepbrasign[expr] = 
  expr /. dirpdt[{a___}, state_] :> predznak[state] nc[ old[a], conj[state] ];

(* vev3[] anticommutates f to the last position in the operator string and
evaluates vev[] of the remainder of the expression. f can be either
annihilation or creation operator; only head f matters. *)
(* Used in matrixel[] and matrixelV2[] *)

vev3[expr_] := Module[{expr2},
  (* Note: there will always be only one f operator*)
  expr2 = (expr /. nc -> nc2) //. 
   { nc2[a___, x1:f[__], x2:op_[__], b___] :> 
       -nc2[a, x2, x1, b] /; fermionQ[op] };
  expr2 = expr2 /. nc2[a__, x:f[__]] :> nc2[a] x;
  expr2 /. nc2[a__] :> vevwicknew[nc[a]]   (* WAS: vev *)
];

(* matrixel[] is called from matrixni[], matrixniNEW[], matrixanomalous[] to
   calculate the offdiagonal matrix element coefficients.

    Note, however, that dodiag[], matrixopisospin[], matrixopspinz[] etc. call
    matrixop1[] instead. CGmatrixel is a symmetry-specific function which
    computes the required Clebsch-Gordan coefficients. *)

(* WARNING, WARNING: In matrixel[] we assume that f and a operators have the same spin.
   Note that the results of matrixel[] are instructions to obtain
   the out-of-diagonal block matrix elements of the Hamiltonian matrix. *)

matrixel[i1_, i2_, ops_] := Module[{expr},
  (* Prepare states (incl. conjugation to get bra) and factor out
     numerical constants. *)
  expr = prepall[i1, i2, ops];
  MyPrint[5,">>>> ",expr];

  (* Anticommute f out of the expression and evaluate vev of the   
     remaining expression. Recall that f is the operator on the PREVIOUS
     site, while a,b,c are operators on the NEW site. *)
  expr = expr //. HoldPattern[ nc[old[x__], c___, old[y__]] ] :>
  nc[OLD[x], vev3[nc[c]], OLD[y]];

  If[expr =!= 0,
    MyPrint[4,"____  ",expr];
  ];

  (* Transform matrix element to IRREDUCIBLE matrix element. Note
     that entire operator f is passed as argument to CGmatrixel[]. *)
  expr = expr /. 
    HoldPattern[ nc[OLD[qn1_,diff1_,r1_], opf:f[__], OLD[qn2_,diff2_,r2_]] ] :>
      CGmatrixel[qn1,diff1,qn2,diff2,opf];

  If[expr =!= 0,
    MyPrint[3,"+___  ",expr];
  ];

  (* Simplify. *)
  expr = simpl5[expr];

  If[expr =!= 0,
    MyPrint[2,"==== ", expr];
  ];

  matrixelsimpl[expr]
];

matrixelsimpl[x_] := x; (* Do nothing by default *)

(* Here we do NOT assume that the two operators carry the same spin
   index! *)

matrixelV2[i1_, i2_, ops_] := Module[{expr},
  MyPrint[2,"matrixelV2[] i1=", i1, " i2=", i2, " ops=", ops];   

  (* Prepare states (incl. conjugation to get bra) and factor out
     numerical constants. *)
  expr = prepall[i1, i2, ops];
  MyPrint[3, "expr=", expr];

  (* Anticommute f out of the expression and evaluate vev of the   
     remaining expression. Recall that f is the operator on the PREVIOUS
     site, while a,b,c are operators on the NEW site. *)
  expr = expr //. HoldPattern[ nc[old[x__], c___, old[y__]] ] :>
  nc[OLD[x], vev3[nc[c]], OLD[y]];

  If[expr =!= 0,
    MyPrint[3, "____  ",expr];
  ];

  (* Transform matrix element to IRREDUCIBLE matrix element. Note
     that entire operator f is passed as argument to CGmatrixel[]. *)
  expr = expr /.
    HoldPattern[ nc[OLD[qn1_,diff1_,r1_], opf:f[_,_,_], OLD[qn2_,diff2_,r2_]] ] :>
      CGmatrixelV2[qn1,diff1,qn2,diff2,opf];

  If[expr =!= 0,
    MyPrint[3,"+___  ",expr];
  ];

  (* Simplify. *)
  expr = simpl5[expr]
];

(* Called from dooffdiag[]. Take into account the hopping matrix
elements for all channels. *)
(* NOTE: simpl4[] is applied on the final result in all routines.
Thus simpl4 must simplify expressions as far as possible, using
FullSimplify[] for example. *)

matrixni[i1_, i2_, ch_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], MyPrint[3,"skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op[AN,UP]] ];
  expr += matrixel[i1, i2, nc[f[CR,ch,DO], op[AN,DO]] ];
  simpl4 @ expr
];

matrixniup[i1_, i2_, ch_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], MyPrint[3,"skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,UP], op[AN,UP]] ];
  simpl4 @ expr
];

matrixnidown[i1_, i2_, ch_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], MyPrint[3,"skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,DO], op[AN,DO]] ];
  simpl4 @ expr
];

(* Different approach, to be used for in2 >= in1 only
in order to avoid double counting! *)

matrixniNEWCR[i1_, i2_, ch_, spin_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], MyPrint[3,"skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[f[CR,ch,spin], op[AN,spin]] ];
  simpl4 @ expr
];

matrixniNEWAN[i1_, i2_, ch_, spin_] := Module[{expr=0, op},
  If[MATRIXNICHECKSYM[i1, i2, ch], MyPrint[3,"skip"]; Return[0]];
  op = ch2op[ch];
  expr += matrixel[i1, i2, nc[op[CR,spin], f[AN,ch,spin]] ];
  simpl4 @ expr
];

matrixniNEW[i1_, i2_, ch_, spin_] := Module[{},
  matrixniNEWCR[i1,i2,ch,spin] + matrixniNEWAN[i1,i2,ch,spin]
];

matrixniNEW[i1_, i2_, ch_] := Module[{},
  simpl4[matrixniNEW[i1,i2,ch,UP] + matrixniNEW[i1,i2,ch,DO]]
];

(* Generic calculation code. Function 'function' is called for all (in1,in2)
combinations and additional outer loop over channel index 'ch' exists. If
upper == True, only upper diagonal blocks will be considered, i.e. in2 >=
in1. Called from dooffdiag[] and other generator functions. *)

(* If INCLUDEZEROONDIAG=True, zero coeffcients are also output for in1=in2. *)

dogeneralloop[function_, prefix_, filename_, upper_:False, INCLUDEZEROONDIAG_:False] := Module[
  {code, ch, in1, startin2, in2, mx, koef, string},
  MyPrint["dogeneralloop[", function, " ", prefix, " ", filename, "]"];
  code = {};
  For[ch = 0, ch < channels, ch++,
    MyPrint["ch=", ch];
    For[in1 = 1, in1 <= nrstates, in1++,
      startin2 = If[upper, in1, 1];
      For[in2 = startin2, in2 <= nrstates, in2++,
        MyPrint[2, in1, " ", in2, " (ch=", ch, ")"];
        mx = function[in1, in2, ch];
        If[mx =!= 0 || (INCLUDEZEROONDIAG && in1 == in2),
          koef = koefstr1[mx];
          string = prefix <> "(" <> 
            tos[in1] <> ", " <> tos[in2] <> ", " <>
            tos[ch] <> ", " <> koef <> ");";
          Print[string];
          AppendTo[code,string];
        ];
      ];
    ];
  ];
  store[code, filename];
];

(* No sum over channels. Performed in function[]. *)
dogeneralloopSUMCH[function_, prefix_, filename_, upper_:False, INCLUDEZEROONDIAG_:False] := Module[
  {code, in1, startin2, in2, mx, koef, string},
  MyPrint["dogeneralloopSUMCH[", function, " ", prefix, " ", filename, "]"];
  code = {};
  For[in1 = 1, in1 <= nrstates, in1++, (* in1 <= nrstates *)
      startin2 = If[upper, in1, 1];
      For[in2 = startin2, in2 <= nrstates, in2++, (* in2 <= nrstates *)
	  MyPrint["##### dogeneralloopSUMCH ##### ", in1, " ", in2];
	  mx = function[in1, in2];
	  If[mx =!= 0 || (INCLUDEZEROONDIAG && in1 == in2),
	     koef = koefstr1[mx];
	     string = prefix <> "(" <>
	       tos[in1] <> ", " <> tos[in2] <> ", " <> koef <> ");";
	     Print[string];
	     AppendTo[code,string];
	     ];
	  ];
      ];
  store[code, filename];
];

(*
Compute the coefficients for building the off-diagonal blocks
of the Hamiltonian matrix. NOTE: the procedure iterates 
through ALL (in1, in2) combinations. This is OK, since
recent implementations of offdiag_function() handle both
in1>in2 and in1<in2 cases by doing necessary transpositions.
For this reason, only one part of the hopping operator is needed,
but not its Hermitian conjugate! 
Alternatively, we use operators which are Hermitian, but only
loop over in2 >= in1, which is what dooffdiagNEW[] does!
*)

dooffdiag[] := Module[{},
  MyPrint["dooffdiag[]"];
  dogeneralloop[matrixni, 
                "OFFDIAG", 
                PREFIX <> "-offdiag.dat", False]; (* default=False *)
];

dooffdiagNEW[] := Module[{},
  MyPrint["dooffdiagNEW[]"];
  dogeneralloop[matrixniNEW, (* different low-level routine! *)
                "OFFDIAG", 
                PREFIX <> "-offdiag.dat", True];
  (* True -> in2 >= in1 *)
];

(* Called from matrixopdiag[], matrixopspinz[], etc.. Computes the matrix
element of an operator op between bra determined by i1 and ket determined by
i2. *)

matrixop1[i1_, i2_, op_] := Module[{expr},
  MyPrint[2, "matrixop1 i1=", i1, " i2=", i2, " op=", op];
  expr = prepall[i1, i2, op];
  MyPrint[4, "1>>> ", expr];
  expr = expr //. 
    HoldPattern[ nc[b___, old[x__], c___, old[y__], d___] ] :> 
    Module[{sk},
      sk = skpdt[{x},{y}]; 
      If[sk == 1, mynewbk[nc[b],nc[c],nc[d]], 0]
    ];
  MyPrint[3, "2>>> ", expr];  
  expr = simpl2 @ expr;
  MyPrint[2, "3>>> ", expr];
  expr
];

(* NOTE: make sure symbols a,b,c,d,e are not used for other purposes (use
in patterns or as local variables is OK, though). *)

ch2op[ch_] := Module[{op},
  op = Null;
  If[ch == 0, op = a];
  If[ch == 1, op = b];
  If[ch == 2, op = c];
  If[ch == 3, op = d];
  If[ch == 4, op = e];
  op
];

(* Called from dodiag[]. *)
matrixopdiag[i_, j_, ch_] := matrixop1[i, j, number[ch2op[ch]] ];

(* Calculate the diagonal matrix elements, i.e. the number of
electrons. This is required in building the Hamiltonian matrix
in the case of p-h asymmetric bands. *)

dodiag[] := Module[{},
  code = {};
  For[ch = 0, ch <channels, ch++,
    For[in1 = 1, in1 <= nrstates, in1++,
      mx = matrixopdiag[in1, in1, ch];
      string = "DIAG(" <> tos[in1] <> ", " <> tos[ch] <> 
        ", " <> toc[mx] <> ");";
      Print[string];
      AppendTo[code,string];
    ];
  ];
  store[code,PREFIX<>"-diag.dat"];
];

(* Irreducible matrix element <||a||> and <||b||>,...
Called from makerecalcf[]. *)

ireducel[i1_, i2_, op_, rule_, cg_] := Module[{st1, st2, expr},
  (* Note: S is Sp in our nomenclature *)
  MyPrint[2, "ireducel i1=", i1, " i2=", i2];
  st1 = prepbra @ (new[i1] /. rule);
  st2 = prepket @ new[i2];
  expr = factorcoef @ nc[st1, op, st2];
  MyPrint[3, "1>>>", expr];
  expr = expr //. HoldPattern[ nc[b___, old[x__], c___, old[y__], d___] ] :>
    Module[{sk},
      sk = skpdt[{x}, {y}]; 
      If[sk == 1, mynewbk[nc[b],nc[c],nc[d]], 0]
    ];
  MyPrint[3, "2>>>", expr];
  expr = expr/cg;
  expr = simpl5 @ expr;
  MyPrint[3, "3>>>", expr];
  expr
];

(* Determine the coefficients for recalculation of the
irreducible matrix elements.  *)
(* 'qn' are quantum numbers of the operator given as 
argument 'op'. 'filename' will contain the resulting 
transformation rules. Be careful not to overwrite things! *)

makerecalcf[op_, qn_, filename_] := Module[
 {rule, cg, code, in1, in2, mx, koef, string},
  MyPrint["makerecalcf op=", op, " qn=", qn, " filename=", filename];
  If[debug >=1,  
    filename2 = filename <> ".mma";
    Put[{}, filename2];
  ];
  rule = irrule @@ qn;
  cg = ircg @@ qn;
  code = {};
  For[in1 = 1, in1 <= nrstates, in1++,
    For[in2 = 1, in2 <= nrstates, in2++,
      mx = ireducel[in1, in2, op, rule, cg];
      MyPrint[{in1,in2} -> mx];
      If[mx =!= 0,
        string = "{ " <> tos[in1] <> ", " <> tos[in2] <> ", " <>
          koefstr2[mx] <> "},";
        AppendTo[code, string];
        If[debug >= 1, PutAppend[{in1, in2, mx}, filename2]];
      ];
    ];
  ];
  If[filename != "",
    store[code, filename <> ".dat"]
  ];
  If[ValueQ[allcode],
    allcode = Join[allcode, code];
  ];
];

(* Distinguish between spin-up and spin-down operators! *)
dorecalcfLOOP["U1"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {}, FNRECALCF <> "-" <> tos[op] <> "-UP"];
    makerecalcf[op[CR, DO], {}, FNRECALCF <> "-" <> tos[op] <> "-DO"];
  ];
];

dorecalcfLOOP["NONE"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    allcode = {};
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {}, FNRECALCF <> "-" <> tos[op] <> "-CR-UP"];
    makerecalcf[op[CR, DO], {}, FNRECALCF <> "-" <> tos[op] <> "-CR-DO"];
    makerecalcf[op[AN, UP], {}, FNRECALCF <> "-" <> tos[op] <> "-AN-UP"];
    makerecalcf[op[AN, DO], {}, FNRECALCF <> "-" <> tos[op] <> "-AN-DO"];
  ];
];

(* recalcf[] in the case of SU(2)_spin symmetry. *)

dorecalcfLOOP["QS"] := Module[{ch, op},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2}, 
                FNSPINUP <> tos[op]];
    makerecalcf[op[CR, DO], {-1/2},
                FNSPINDOWN <> tos[op] ];
  ];
];

dorecalcfLOOP["QSLR", p_] := Module[{str, ch, op},
  str = If[p == 1, "", "diff"];
  For[ch = 0, ch <= 1, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2, p}, 
                FNSPINUP <> str <> tos[op]];
    makerecalcf[op[CR, DO], {-1/2, p},
                FNSPINDOWN <> str <> tos[op] ];
  ];
];

dorecalcfLOOP["QSLR"] := Module[{},
  dorecalcfLOOP["QSLR", 1];
  dorecalcfLOOP["QSLR", -1];
];

dorecalcfLOOP["ISO"] := Module[{ch, op, nn},
  For[ch = 0, ch < channels, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2, 1/2},
                FNSPINUPISOUP <> tos[op]];
    makerecalcf[op[CR, DO], {1/2, -1/2},
                FNSPINDOWNISOUP <> tos[op] ];
    makerecalcf[op[AN, DO], {-1/2, 1/2}, 
                FNSPINUPISODOWN <> tos[op]];
    makerecalcf[op[AN, UP], {-1/2, -1/2},
                FNSPINDOWNISODOWN <> tos[op] ];
  ];
];

dorecalcfLOOP["ISOLR", p_] := Module[{str, ch, op},
  str = If[p == 1, "", "diff"];
  For[ch = 0, ch <= 1, ch++,
    op = ch2op[ch];
    makerecalcf[op[CR, UP], {1/2, 1/2, p},
                FNSPINUPISOUP <> str <> tos[op]];
    makerecalcf[op[CR, DO], {1/2, -1/2, p},
                FNSPINDOWNISOUP <> str <> tos[op] ];
    makerecalcf[op[AN, DO], {-1/2, 1/2, p}, 
                FNSPINUPISODOWN <> str <> tos[op]];
    makerecalcf[op[AN, UP], {-1/2, -1/2, p},
                FNSPINDOWNISODOWN <> str <> tos[op] ];
  ];
];

dorecalcfLOOP["ISOLR"] := Module[{},
  dorecalcfLOOP["ISOLR", 1];
  dorecalcfLOOP["ISOLR", -1];
];

dorecalcf[sym_:"QS"] := Module[{ch, op},
  dorecalcfLOOP[sym];
];

dorecalcf[sym_:"QS", nr_] := Module[{ch, op},
  dorecalcfLOOP[sym, nr];
];

checkdiff[d_, diffs_] := (MyVPrint[3, "cd", d, diffs]; If[MemberQ[diffs, d], 1, 0]);

(*
'm' characterizes the symmetry properties of the operator considered 
(such as spin of the operator). NOTE: the operators must
be tensor operators with respect to all symmetries taken into account!

'diffs' is the list of differences in quantum numbers and associated
filenames for the results. 

QNrecalcop[] determines the relation between quantum numbers of
bra and ket.

RULE1recalcop[] and RULE2recalcop[] are related to QNrecalcop and
performs the necessary mapping.

CG1recalcop[] and CG2recalcop[] are the Clebsch-Gordan coefficients
that appear in the denominator and numberator of the recursion
formula.
*)

(* TO DO: RULE1 argument m. m in others, too! *)

(* A single value of mz in {m,mz} is sufficient due to the Wigner-Eckart
theorem! *)

checksignfn[{m_}] := !(m \[Element] Integers); (* True for fermionic operators *)

recalcopONE[m_, diffs_, {diffone_, filename_}] := Module[
  {seznam, checksign, i, j, expr, string},

  MyPrint[3,"recalcopONE"];
  MyPrint[3,"m=",m];
  MyPrint[3,"diffs=",diffs];
  MyPrint[3,"diffone=", diffone];
  Myprint[3,"filename=", filename];

  seznam = {}; (* In this array we accumulate the results *)
  seznamdebug = {}; (* Details of the calculation.. *)

  checksign = checksignfn[m]; (* True for fermionic operators *)
  string = "checksign=" <> toi[checksign];
  AppendTo[seznamdebug, string];

  QNrecalcop[ diffone ]; (* Relation between QNs in bra and QNs in ket. *)

  For[i = 1, i <= nrstates, i++, (* All combinations are considered! *)
    MyPrint[2, "############## i=", i];
    st1 = new[i] /. RULE1recalcop[m];
    If[checksign,
      st1 = prepbrasign @ st1,
      st1 = prepbra @ st1
    ];
    st1 = st1 /. nc[coef[0], ___] -> 0;
    MyPrint[3, "st1=", st1];

    AppendTo[seznamdebug, ""];
    string = "*** <i=" <> toi[i] <> "|=" <> toi[st1];
    AppendTo[seznamdebug, string];

    For[j = 1, j <= nrstates, j++,
      MyPrint[2, "------------ j=", j];
      st2 = prepket[ new[j] /. RULE2recalcop[m]];
      st2 = st2 /. nc[coef[0], ___] -> 0;
      MyPrint[3, "st2=", st2];

      AppendTo[seznamdebug, ""];
      string = "|j=" <> toi[j] <> ">=" <> toi[st2];
      AppendTo[seznamdebug, string];

      expr = factorcoef @ nc[st1, st2];
      expr = expr //. HoldPattern[ nc[b___, old[x__], c___, old[y__], d___] ] :>
        pdt[{x}, {y}] mynewbk[nc[b], nc[c], nc[d]];

      MyPrint[2, i, " ", j, " expr=", expr];
      string = "expr=" <> toi[expr];  
      AppendTo[seznamdebug, string];

      cg1 = CG1recalcop[m];
      string = "cg1=" <> toi[cg1];
      AppendTo[seznamdebug, string];

      expr = expr /. pdt[{qn1_, diff1_, r1_}, {qn2_, diff2_, r2_}] :>
        checkdiff[diffqn[qn1, qn2], diffs] ireduc[i, j, qn1, qn2] *
        CG2recalcop[qn1, diff1, qn2, diff2, m];

      MyPrint[3, "expr(a)=", expr];  
      expr = expr/cg1;
      MyPrint[3, "expr(b)=", expr];
      expr = simpl4[expr];
      MyPrint[3, "expr(c)=", expr];

      string = "final expr=" <> toi[expr];
      AppendTo[seznamdebug, string];

      If[expr =!= 0,
        string = expr /. faktor_. ireduc[i1_, ip_, IN1_, INp_] :>
          "{ " <> tos[i1] <> ", " <> tos[ip] <> ", " <>
          outmake[IN1] <> ", " <> outmake[INp] <>
          ", " <> koefstr3[ simpl4[faktor] ] <> "},";
        Print["#### string=", string];

        AppendTo[seznamdebug, string];
        AppendTo[seznam, string];
      ];
    ];
  ];
  If[debug >= 1, store[seznamdebug, filename <> "-debug"]];
  store[seznam, filename];
];

(* 
'm' is a set of quantum numbers that characterize the
symmetry properties of the operator; one quantity per 
symmetry group! It is used as the argument to various
auxiliary functions (checksignfn[], CG[12]recalcopLOOP[],
RULE[12]recalcopLOOP[]).

'list' is a list of {qn, filename} pairs for all possible
differences in the values of quantum numbers of bras and kets.
*)

recalcop[m_, list_] := Module[{diffs},
  (* Differences of conserved quantum numbers (i.e. those that characterize 
     the invariant subspaces. *)
  diffs = list[[All, 1]];
  MyPrint["diffs=", diffs];
  Scan[recalcopONE[m, diffs, #]&, list];
];

(* matrixopdiagSUMCH[] must be defined! *)
dodiagSUMCH[] := Module[{},
   code = {};
   For[in1 = 1, in1 <= nrstates, in1++,
       mx = matrixopdiagSUMCH[in1, in1];
       string = "DIAG(" <> tos[in1] <> ", " <> toc[mx] <> ");";
       Print[string];
       AppendTo[code,string];
   ];
   store[code,PREFIX<>"-diag.dat"];
];

(* Sum over channels performed in the inner-most loop, i.e.,
 in matrixniSUMCH[], rather than in dogeneralloop[]. *)
dooffdiagSUMCH[] := Module[{},
    MyPrint["dooffdiagSUMCH[]"];
    dogeneralloopSUMCH[matrixniSUMCH, (* different low-level routine! *)
		       "OFFDIAG",
		       PREFIX <> "-offdiag.dat", True];
			   (* True -> in2 >= in1 *)
];

(* New implementation, Oct 2016 *)

recalcopX[m_, list_] := Module[{diffs},
  MyPrint["recalcopX[]"];
  (* Differences of conserved quantum numbers (i.e. those that characterize
     the invariant subspaces. *)
  diffs = list[[All, 1]];
  MyPrint["diffs=", diffs];
  Scan[recalcopONEX[m, diffs, #]&, list];
];

(* JUST ONE ! *)
recalcopX[m_, list_, nr_] := Module[{diffs},
  MyPrint["recalcopX[]"];
  (* Differences of conserved quantum numbers (i.e. those that characterize
     the invariant subspaces. *)
  diffs = list[[All, 1]];
  MyPrint["diffs=", diffs];
  recalcopONEX[m, diffs, list[[nr]] ];
];

(* check[] must be defined. Example:
check[{a_, b_, {d_, e_, f_}, {g_, h_, i_}, j_}] := True;
check[arg___] := (Print["Check failed:", args[arg]]; Exit[1]);
*)

recalcopONEX[m_, diffs_, {diffone_, filename_}] := Module[
  {seznam, checksign, i, j, expr, string},

  fn2 = filename <> ".mma";
  fn3 = filename <> ".debug";

  seznam = {}; (* In this array we accumulate the results *)

  checksign = checksignfn[m]; (* True for fermionic operators *)
  MyPrint["recalcopONEX[] ", m, diffs, {diffone, filename}, {fn2, fn3}, "checksign=", checksign];

  (* Relation between QNs in bra and QNs in ket. This defines Qp,Sp,Tp in terms of Q,S,T plus
     the differences diffone={deltaq,deltas,delta). *)
  QNrecalcop[ diffone ];

  cnt = 0;

  For[i = 1, i <= nrstates, i++,  (* All combinations are considered! *)
    st1 = new[i] /. RULE1recalcop[m];
    If[checksign,
      st1 = prepbrasign @ st1,
      st1 = prepbra @ st1
    ];
    st1 = st1 /. nc[coef[0], ___] -> 0;
    st1 = simpl4[st1]; (* 21.4.2016 *)

    For[j = 1, j <= nrstates, j++,
      st2 = prepket[ new[j] /. RULE2recalcop[m]];
      st2 = st2 /. nc[coef[0], ___] -> 0;

      expr = factorcoef @ nc[st1, st2];
      expr = expr //. HoldPattern[ nc[b___, old[x__], c___, old[y__], d___] ] :>
        pdt[{x}, {y}] mynewbk[nc[b], nc[c], nc[d]];

      MyPrint[{i, j}];

      (* If evaluates only the argument determined by the value of the condition. *)
      expr = expr /. pdt[{qn1_, diff1_, r1_}, {qn2_, diff2_, r2_}] :>
        If[checkdiff[diffqn[qn1, qn2], diffs] == 1, ireduc[i, j, qn1, qn2] * CG2recalcop[qn1, diff1, qn2, diff2, m], 0];

     If[expr =!= 0, (* Don't waste time *)

      cg1 = CG1recalcop[m];
      expr = expr/cg1;

      expr = simpl4[expr];

      If[expr =!= 0,
        cnt = cnt + 1;
        string = {};

        p = Position[expr, ireduc[_,_,_,_]];
        el = Extract[expr, p];
        el = Union[el]; (* Must be unique! *)
        If[Length[el] != 1, Print["error el=", el]; Exit[1]];
        expr = expr //. ireduc[_,_,_,_] -> 1;
        res = Append[List @@ First[el], expr];

        PutAppend[res, fn2]; (* simpl form *)
        MyPrint["res=", res];

        PutAppend[expr, fn3]; (* debug: unprocessed form *)

        check[res];

        AppendTo[seznam, string];
        ]; (* If =!= 0 *)
      ]; (* If =!= 0 *)
    ];  (* j *)
  ]; (* i *)
  store[seznam, filename];
];
