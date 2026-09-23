tt={};

tt = Join[tt, mtSingletOp["SigmaHartree",    SigmaHartreeAvg ]];
tt = Join[tt, mtSingletOp["SigmaHartree-uu", SigmaHartree[UP, UP] ]];
tt = Join[tt, mtSingletOp["SigmaHartree-dd", SigmaHartree[DO, DO] ]];
tt = Join[tt, mtSingletOp["SigmaHartree-ud", SigmaHartree[UP, DO] ]];
tt = Join[tt, mtSingletOp["SigmaHartree-du", SigmaHartree[DO, UP] ]];

tt
