def1ch[1];

nph = ToExpression @ optionvalue["Nph"];

Himp = delta number[d[]] + U/2 pow[number[d[]]-1, 2];
Himp = Himp + omega phononnumber[nph] + g nc[number[d[]]-1, phononplus[nph] + phononminus[nph]];
MAKEPHONON = 1; (* One phonon mode *)

PR[i_] := PR[i] = nc[ket[i], bra[i]];
PRlast = PR[nph];
 
H = Himp + Hc + H0;
