def2ch[2];

H1 = eps1 number[d[]] + U1 hubbard[d[]];
H2 = eps2 number[a[]] + U2 hubbard[a[]];
Hc1 = gammaPol1 hop[f[0], d[]];
Hc2 = gammaPol2 hop[f[1], a[]];

H = H0 + H1 + Hc1 + H2 + Hc2;

params = Join[params, {
  gammaPol1 -> Sqrt[(1/Pi) thetaCh[1] (extraGamma1) gammaA],
  gammaPol2 -> Sqrt[(1/Pi) thetaCh[2] (extraGamma2) gammaA]
}];
