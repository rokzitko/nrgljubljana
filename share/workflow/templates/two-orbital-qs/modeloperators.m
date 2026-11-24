(* SPDX-License-Identifier: GPL-3.0-or-later *)
Join[
  mtSingletOp["n_d1", number[d[1]]],
  mtSingletOp["n_d2", number[d[2]]],
  mtSingletOp["D_d1", hubbard[d[1]]],
  mtSingletOp["D_d2", hubbard[d[2]]],
  mtSingletOp["n_d1n_d2", nc[number[d[1]], number[d[2]]]],
  mtSingletOp["S_d1S_d2", spinspin[d[1], d[2]]],
  mtSingletOp["SigmaHartree11", hartree[1, 1, selfopd1]],
  mtSingletOp["SigmaHartree12", hartree[1, 2, selfopd2]],
  mtSingletOp["SigmaHartree21", hartree[2, 1, selfopd1]],
  mtSingletOp["SigmaHartree22", hartree[2, 2, selfopd2]],
  mtDoubletOp["A_d1", d[1]],
  mtDoubletOp["A_d2", d[2]],
  mtDoubletOp["self_d1", selfopd1],
  mtDoubletOp["self_d2", selfopd2],
  mtTripletOp["sigma_d1", d[1]],
  mtTripletOp["sigma_d2", d[2]]
]
