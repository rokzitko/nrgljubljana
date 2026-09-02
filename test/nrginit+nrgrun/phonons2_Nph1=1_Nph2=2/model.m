def1ch[1];

nph1 = ToExpression @ optionvalue["Nph1"];
nph2 = ToExpression @ optionvalue["Nph2"];
cutoffs = {nph1, nph2};

Himp = delta number[d[]] + U/2 pow[number[d[]]-1, 2];
Himp = Himp + omega1 phononnumber[1, cutoffs] + omega2 phononnumber[2, cutoffs];
Himp = Himp + g1 nc[number[d[]]-1, phononx[1, cutoffs]];
Himp = Himp + g2 nc[number[d[]]-1, phononx[2, cutoffs]];
MAKEPHONON = 2; (* Two phonon modes *)

H = Himp + Hc + H0;
