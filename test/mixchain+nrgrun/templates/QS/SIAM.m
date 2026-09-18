def1ch[nrimp=1];
Heps = eps1 number[d[]];
Hint = U1 hubbard[d[]];
Himp = Heps + Hint;
Hhyb = gammaPolCh[1] hop[f[0], d[]];
H = H0 + Himp + Hhyb;
selfopd = ( Chop @ Expand @ komutator[ Hint, d[#1, #2] ] )&;

SigmaHartree = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #1] ] &;
SigmaHartreeAvg := Expand @ (SigmaHartree[UP] + SigmaHartree[DO]) / 2;

Print["SigmaHartree[UP]=", SigmaHartree[UP] ];
Print["SigmaHartree[DO]=", SigmaHartree[DO] ];
Print["SigmaHartree=", SigmaHartreeAvg ];