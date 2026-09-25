(* Test-only finite-star injection, before the real normalization/reconstruction. *)
Clear[de, deminus, df, dfminus, thetaCh, theta0Ch];
starRead[file_] := setpr[SetPrecision[Flatten[Import[file, "Table"]], Infinity]];
starEp = starRead["de_pos.dat"];
starEm = starRead["de_neg.dat"];
starUp = starRead["du_pos.dat"];
starUm = starRead["du_neg.dat"];
starTheta = First[starRead["theta.dat"]];
If[Length /@ {starEp, starEm, starUp, starUm} =!= {2, 2, 2, 2}, Exit[1]];
mMAX = 1;
thetaCh[a_] = starTheta;
de[a_, m_] := starEp[[m+1]];
deminus[a_, m_] := starEm[[m+1]];
df[a_, m_] := thetaCh[a] starUp[[m+1]]^2;
dfminus[a_, m_] := thetaCh[a] starUm[[m+1]]^2;
