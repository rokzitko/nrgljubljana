def1ch[1];

snegrealconstants[eps, U];

Hpot = U hubbard[d[]];
H1 = eps number[d[]] + Hpot;
Himp = H1;
Hselfd = Hpot;

H = H0 + H1 + Hc;

selfopd = ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;

(* Evaluate *)
Print["selfopd[CR,UP]=", selfopd[CR, UP]];
Print["selfopd[CR,DO]=", selfopd[CR, DO]];
Print["selfopd[AN,UP]=", selfopd[AN, UP]];
Print["selfopd[AN,DO]=", selfopd[AN, DO]];

SigmaHd = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #1] ] /. params &;
SigmaHdAvg := Expand @ (SigmaHd[UP]+SigmaHd[DO])/2;

Print["SigmaH[UP]=", SigmaHd[UP] ];
Print["SigmaH[DO]=", SigmaHd[DO] ];
Print["SigmaH=", SigmaHdAvg ];
