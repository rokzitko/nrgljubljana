#!/bin/sh
# Compare a scalar mixchain chain with the one nrgchain computed from the same input.
#
#   compare_nrgchain.sh TOLERANCE
#
# run in a directory that holds nrgchain's theta.dat, xi.dat and zeta.dat and mixchain's V11.dat, T11.dat and E11.dat
# (discretization_files=true). For one channel the polar gauge makes V and every T_n positive, so the numbers are
# directly comparable, and both are in the units of the input: nrgchain multiplies xi and zeta by bandrescale, and
# mixchain writes E and T the same way; theta and V^2 are both the integral of the input over the range the mesh
# covers.
#
# The comparison is relative, since mycomp.pl's absolute tolerance cannot see an error in xi_n at the late sites:
#   theta:  |V^2 - theta| <= TOLERANCE * theta
#   xi_n:   |T_n - xi_n|  <= TOLERANCE * xi_n
#   zeta_n: |E_n - zeta_n| <= TOLERANCE * max(|zeta_n|, xi_n), on the scale of the site, since zeta_n vanishes for a
#           particle-hole symmetric input and is then rounding.
set -eu
tolerance=$1

for f in theta.dat xi.dat zeta.dat V11.dat T11.dat E11.dat; do
  [ -f "$f" ] || { echo "compare_nrgchain: $f is missing"; exit 1; }
done

awk -v tol="$tolerance" '
  function abs(x) { return x < 0 ? -x : x }
  function max(a, b) { return a > b ? a : b }
  FILENAME == "theta.dat" { theta = $1; next }
  FILENAME == "V11.dat"   { v = $1; next }
  FILENAME == "xi.dat"    { xi[nxi++] = $1; next }
  FILENAME == "zeta.dat"  { zeta[nzeta++] = $1; next }
  FILENAME == "T11.dat"   { t[nt++] = $1; next }
  FILENAME == "E11.dat"   { e[ne++] = $1; next }
  END {
    failed = 0
    if (nxi != nt || nzeta != ne) {
      printf "row counts differ: xi %d, T11 %d, zeta %d, E11 %d\n", nxi, nt, nzeta, ne
      exit 1
    }
    d = abs(v * v - theta) / theta
    printf "theta: nrgchain %.17g, mixchain V11^2 %.17g, relative difference %.2e\n", theta, v * v, d
    if (d > tol) failed = 1

    worst_xi = 0; worst_zeta = 0
    for (n = 0; n < nxi; n++) {
      dx = abs(t[n] - xi[n]) / xi[n]
      if (dx > worst_xi) { worst_xi = dx; at_xi = n }
      dz = abs(e[n] - zeta[n]) / max(abs(zeta[n]), xi[n])
      if (dz > worst_zeta) { worst_zeta = dz; at_zeta = n }
    }
    printf "xi:    %d sites, largest relative difference %.2e at site %d\n", nxi, worst_xi, at_xi
    printf "zeta:  %d sites, largest difference relative to the site %.2e at site %d\n", nzeta, worst_zeta, at_zeta
    if (worst_xi > tol || worst_zeta > tol) failed = 1
    if (failed) { printf "FAILED: tolerance %g\n", tol; exit 1 }
    printf "passed: tolerance %g\n", tol
  }' theta.dat V11.dat xi.dat zeta.dat T11.dat E11.dat
