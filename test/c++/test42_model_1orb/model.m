def1ch[1];

H1 = eps1 number[d[]] + U1 hubbard[d[]];
Hc1 = gammaPol1 hop[f[0], d[]];

H = H0 + H1 + Hc1;

params = Join[params, {
  gammaPol1 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma1) gammaA]
}];
