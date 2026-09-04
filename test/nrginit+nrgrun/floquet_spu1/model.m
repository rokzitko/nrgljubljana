def1ch[1];

MAKEFLOQUET = 1;
nph = 0;
MAKEPHONONTENSOR = 1;
snegrealconstants[A, Omega];

Hfl = A floquetlift[number[d[]], 1, ncut] +
      A floquetlift[number[d[]], -1, ncut] +
      Omega floquetnumber[ncut];
opm = floquetnumber[ncut];
opm2 = pow[opm, 2];

H1 = delta number[d[]] + U/2 pow[number[d[]]-1, 2] + B spinz[d[]];
H = H0 + Hc + H1 + Hfl;
