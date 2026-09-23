(* SIAM with U(1) charge symmetry only.

   Spin is not conserved: the impurity may carry a magnetic field of
   arbitrary direction in the x-z plane, and (with pol2x2=true) the Wilson
   chain has a full 2x2 matrix structure in the spin space, i.e. the chains
   for spin up and spin down are coupled by cross hopping coefficients. *)

def1ch[nrimp=1];

(* coefV[i,j] is read by "matrix -s" from V{i}{j}{ch}.dat and is real; it is
   not the hybV of nrginit, and carries no factor of 1/Sqrt[Pi], so the files
   hold the physical amplitude directly. Declaring it real makes Conjugate
   collapse -- the parser has no Conjugate. *)
snegrealconstants[coefV[1, 1], coefV[1, 2], coefV[2, 1], coefV[2, 2]];

Heps = eps1 number[d[]] + Bz1 spinz[d[]] + Bx1 spinx[d[]];
Hint = U1 hubbard[d[]];
Himp = Heps + Hint;

(* H_hyb = Sum_{i,j} V[i,j] f^dag_{0,i} d_j + h.c., first index the bath
   channel and second the impurity spin, so that Delta = V^dag V g.

   Hc from HC[] is not used: it is built from gammaPolCh = Sqrt[theta/Pi],
   which is non-negative by construction (tools/matrix/parser.cc asserts
   theta >= 0), so it cannot represent a V with a negative or vanishing
   off-diagonal entry -- and the discretization produces exactly that, since
   V is fixed only up to a unitary on the bath index.

   Each term carries its own conjugate. The 4-argument
   hop[a,b,s1,s2] = a^dag_s1 b_s2 + b^dag_s2 a_s1 pairs V[1,2] with itself
   rather than with its conjugate, and pairing V[1,2] with V[2,1] is Hermitian
   only for a Hermitian V. genhop does the right thing but takes a single spin,
   so it covers only the diagonal. *)
Hhyb = genhop[coefV[1, 1], f[0], d[], UP] +
       genhop[coefV[2, 2], f[0], d[], DO] +
       coefV[1, 2] nc[f[CR, 0, UP], d[AN, DO]] +
       Conjugate[coefV[1, 2]] nc[d[CR, DO], f[AN, 0, UP]] +
       coefV[2, 1] nc[f[CR, 0, DO], d[AN, UP]] +
       Conjugate[coefV[2, 1]] nc[d[CR, UP], f[AN, 0, DO]];

H = H0 + Himp + Hhyb;

Print["Hhyb=", Hhyb];

selfopd = ( Chop @ Expand @ komutator[ Hint, d[#1, #2] ] )&;

(* High-frequency (Hartree) limit of the auxiliary correlator F used in the
   self-energy trick. With a spin-mixing bath the self-energy is a 2x2 matrix
   in the spin space, hence all four components are required. *)
(* Chop is required: the anticommutator leaves a residual `0.` additive term,
   which op2matrix cannot digest (First/Last of a bare scalar), and nrginit then
   writes an unevaluated expression into `data`. *)
SigmaHartree = ( Chop @ Expand @ antikomutator[ selfopd[CR, #1], d[AN, #2] ] )&;
SigmaHartreeAvg := Chop @ Expand @ ((SigmaHartree[UP, UP] + SigmaHartree[DO, DO]) / 2);

Print["SigmaHartree[UP,UP]=", SigmaHartree[UP, UP] ];
Print["SigmaHartree[DO,DO]=", SigmaHartree[DO, DO] ];
Print["SigmaHartree[UP,DO]=", SigmaHartree[UP, DO] ];
Print["SigmaHartree[DO,UP]=", SigmaHartree[DO, UP] ];
Print["SigmaHartree=", SigmaHartreeAvg ];
