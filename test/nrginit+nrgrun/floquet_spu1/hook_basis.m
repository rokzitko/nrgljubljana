Module[{hookbz},
  hookbz = transformtoPHtensor[bzvc2bzop[bvc], {0}];
  bvc = bzop2bzvc[hookbz, vak];
]
