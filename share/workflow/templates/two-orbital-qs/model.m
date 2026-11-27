(* SPDX-License-Identifier: GPL-3.0-or-later *)
def2ch[nrimp = 2];
impuritybasis = {d[1], d[2]};

SetAttributes[bathCoupling, NumericFunction];
bathCoupling /: Conjugate[bathCoupling[i_, j_]] := bathCoupling[i, j];

pairTransfer[i_, j_] := nc[d[CR, i, UP], d[CR, i, DO], d[AN, j, DO], d[AN, j, UP]];
Hint = U1 hubbard[d[1]] + U2 hubbard[d[2]] +
       (U12 - J12/2) nc[number[d[1]], number[d[2]]] -
       2 J12 spinspin[d[1], d[2]] + J12 (pairTransfer[1, 2] + pairTransfer[2, 1]);
Heps = eps1 number[d[1]] + eps2 number[d[2]] + V12 hop[d[1], d[2]];
Hhyb = Sum[bathCoupling[i, j] hop[d[i], f[j - 1]], {i, 1, 2}, {j, 1, 2}];
H = H0 + Heps + Hint + Hhyb;

selfopd1 = (Chop @ Expand @ komutator[Hint, d[#1, 1, #2]]) &;
selfopd2 = (Chop @ Expand @ komutator[Hint, d[#1, 2, #2]]) &;
hartree[i_, j_, q_] := Expand @ Sum[antikomutator[d[AN, i, spin], q[CR, spin]], {spin, {UP, DO}}]/2;
