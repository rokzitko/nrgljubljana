Modules[{ t={} },
  (* 'Occupancy' of electron gas on the dot *)
  t = Join[t, mtSingletOp["n_F", Sum[number[f[tz]], {tz, -1, 1}] ]];
  
  texportable = t;
];

MyPrint["modelopetators.m DONE"];

texportable
