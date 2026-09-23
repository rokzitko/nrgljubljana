def1ch[1];

(* Impurity Hamiltonian *)
snegrealconstants[eps1, U1, coefV[1, 1], coefV[1, 2]];

(* Currently the parser does not support conjugated coefficients *)
Conjugate[coefdelta[i__]] ^= coefdelta[i];
(* Define the model parameters *)

(* Define the parameters for the first site of the Wilson chain *)

(* Define the Hamiltonian *)

(* Define the Hamiltonian for the first site of the Wilson chain *)

Heps = eps1 number[d[]];
Hint = U1 hubbard[d[]];

Himp = Heps + Hint;

(* Hybridization Hamiltonian *)

Hhyb = genhop[coefV[1, 1], f[0], d[]] + genanhop[coefV[1, 2], f[0], d[]];

(* Total Hamiltonian *)
(* H0 is the first site of Wilson chain *)
H = H0 + Himp + Hhyb;

(* Define the operators *)

(* Definition of auxilary operators *)
(* Only the interacting part! Hartree term needs to be calculated separately *)
Hselfd = Hint;

selfopd =  ( Chop @ Expand @ komutator[Hselfd, d[#1, #2]] )&;
selfopcd = ( Chop @ Expand @ komutator[Hselfd, ((-1)^#2 d[1-#1, 1-#2])] )&;

(* Constant Hartree term *)
SigmaHd = Expand @ antikomutator[ selfopd[#1, #2], d[#3, #4] ] &;

SigmaHdAvg11 := Expand @ (SigmaHd[CR, UP, AN, UP]+SigmaHd[CR, DO, AN, DO])/2;
SigmaHdAvg12 := Expand @ (SigmaHd[AN, DO, AN, UP]+SigmaHd[CR, DO, CR, UP])/2;
SigmaHdAvg21 := Expand @ (SigmaHd[CR, UP, CR, DO]+SigmaHd[AN, UP, AN, DO])/2;
SigmaHdAvg22 := Expand @ (SigmaHd[AN, DO, CR, DO]+SigmaHd[AN, UP, CR, UP])/2;

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];

Print["selfopcd[CR,UP]=", selfopcd[CR, UP]];
Print["selfopcd[CR,DO]=", selfopcd[CR, DO]];
Print["selfopcd[AN,UP]=", selfopcd[AN, UP]];
Print["selfopcd[AN,DO]=", selfopcd[AN, DO]];

Print["SigmaH11[CR, UP, AN, UP]=", SigmaHd[CR, UP, AN, UP] ];
Print["SigmaH12[AN, DO, AN, UP]=", SigmaHd[AN, DO, AN, UP] ];
Print["SigmaH21[CR, UP, CR, DO]=", SigmaHd[CR, UP, CR, DO] ];
Print["SigmaH22[AN, DO, CR, DO]=", SigmaHd[AN, DO, CR, DO] ];

Print["SigmaHdAvg11=", SigmaHdAvg11];
Print["SigmaHdAvg12=", SigmaHdAvg12];
Print["SigmaHdAvg21=", SigmaHdAvg21];
Print["SigmaHdAvg22=", SigmaHdAvg22];