def2ch[2];
nnop[f[1]] = 0;
nnop[a[]] = -1;
nnop[f[0]] = 0;
nnop[d[]] = -1;
lrchain = {f[0], d[], a[], f[1]};

H1 = U1 hubbard[d[]] + epsilon1 number[d[]] + B spinz[d[]];
H2 = U2 hubbard[a[]] + epsilon2 number[a[]] + B spinz[a[]];
Hc = gammaPol1 hop[f[0], d[]] + gammaPol2 hop[f[1], a[]];

sz = spinz[d[]] + spinz[a[]];
H12 = -J spinspin[a[], d[]] + anD pow[sz, 2] + U12 nc[number[d[]], number[a[]]];

H = H0 + H1 + H2 + H12 + Hc;
Hhyb = Hc;

Hselfd = H1 + H12;
Hselfa = H2 + H12;
selfopd = ( Chop @ Expand @ komutator[Hselfd /. params, d[#1, #2]] )&;
selfopa = ( Chop @ Expand @ komutator[Hselfa /. params, a[#1, #2]] )&;

Print["self_1[CR,UP]=", selfopd[CR,UP] ];
Print["self_2[CR,UP]=", selfopa[CR,UP] ];

SigmaHd = Expand @ antikomutator[ selfopd[CR, #1], d[AN, #1] ] &;
SigmaHa = Expand @ antikomutator[ selfopa[CR, #1], a[AN, #1] ] &;

Print["SigmaH_1[UP]=", SigmaHd[UP] ];
Print["SigmaH_1[DO]=", SigmaHd[DO] ];
Print["SigmaH_2[UP]=", SigmaHa[UP] ];
Print["SigmaH_2[DO]=", SigmaHa[DO] ];

params = Join[params, {
gammaPol1 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma1) gammaA],
gammaPol2 -> Sqrt[(1/Pi) thetaCh[2] (extraGamma2) gammaA]
}];
