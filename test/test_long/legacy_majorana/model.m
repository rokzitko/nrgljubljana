def1ch[1];

op1 = nc[ (f[CR, 0, DO]+f[AN,0, DO])/Sqrt[2], (d[CR, DO] + d[AN, DO])/Sqrt[2] ];

Hc = gammaPolCh[1] (hop[f[0], d[], UP] + op1 + conj[op1]);

H = H0 + H1 + Hc;
