#!/bin/bash
# The hybridisation of an s-wave BCS bath in Nambu space, psi = (c_up, c_down^dagger), as in 01_star_test/016 but
# with the anomalous part in the other basis phase:
#
#   Gamma(omega) = rho0/sqrt(omega^2 - delta^2) [[ |omega|,  -delta sgn(omega) ],
#                                                [ -delta sgn(omega),  |omega| ]]     for |omega| > delta,
#   Gamma(omega) = 0                                                                  inside the gap.
#
# The sign of the anomalous part is a choice of Nambu basis phase: psi = (c_up, c_down^dagger) against
# (c_up, -c_down^dagger) turns it over, and with it E(1,2) of the chain and the sign of the induced <pair_d>. The
# physics -- the gap, the spectrum, |<pair_d>| -- does not depend on it.
#
# This case measured which of the two phases matches nrg's scdelta, whose own sign comes from ISOSPINX and conj_me:
# it is this one, negative above the gap. Against the established flat-band-plus-constant-gap route, it reproduces
# <pair_d> in sign and to 0.3% in size, while the other phase gives the same magnitude with the opposite sign. The
# staging step passes E(1,2) through unchanged; the convention is stated here, where it belongs, since a user
# supplying their own Gamma chooses it themselves.
#
# Sampled from the gap outward on a mesh geometric in u = |omega| - delta, so the square-root edge is resolved. The
# first node is the gap edge itself with Gamma = 0; there are no nodes inside the gap, and below the edge the loader
# continues the input at that zero, so the gap stays empty.
#
# The normal-state weight is theta = 2 rho0 sqrt(1 - delta^2), which the reference route needs in order to put the
# same hybridisation strength on its flat band.
set -eu
RHO0=${RHO0:-0.1}
DELTA=${DELTA:-0.01}
D=${D:-1.0}
UMIN=${UMIN:-1e-15}
RATIO=${RATIO:-1.01}

awk -v rho0=$RHO0 -v delta=$DELTA -v D=$D -v umin=$UMIN -v ratio=$RATIO 'BEGIN {
  n = 0;
  grid[n] = delta; diag[n] = 0; off[n] = 0; n++;
  for (u = umin; u < D - delta; u *= ratio) {
    w = delta + u;
    grid[n] = w; diag[n] = rho0 * w / sqrt(w * w - delta * delta); off[n] = rho0 * delta / sqrt(w * w - delta * delta);
    n++;
  }
  grid[n] = D; diag[n] = rho0 * D / sqrt(D * D - delta * delta); off[n] = rho0 * delta / sqrt(D * D - delta * delta);
  n++;

  for (i = n - 1; i >= 0; i--) printf "%.17g %.17g\n", -grid[i], diag[i] > "Gamma_11-re.dat";
  for (i = 0; i < n; i++)      printf "%.17g %.17g\n",  grid[i], diag[i] > "Gamma_11-re.dat";
  for (i = n - 1; i >= 0; i--) printf "%.17g %.17g\n", -grid[i], diag[i] > "Gamma_22-re.dat";
  for (i = 0; i < n; i++)      printf "%.17g %.17g\n",  grid[i], diag[i] > "Gamma_22-re.dat";
  # the anomalous part is odd in omega, and negative above the gap: the phase that matches scdelta
  for (i = n - 1; i >= 0; i--) printf "%.17g %.17g\n", -grid[i],  off[i] > "Gamma_12-re.dat";
  for (i = 0; i < n; i++)      printf "%.17g %.17g\n",  grid[i], -off[i] > "Gamma_12-re.dat";
  for (i = n - 1; i >= 0; i--) printf "%.17g %.17g\n", -grid[i],  off[i] > "Gamma_21-re.dat";
  for (i = 0; i < n; i++)      printf "%.17g %.17g\n",  grid[i], -off[i] > "Gamma_21-re.dat";

  printf "make_input: rho0=%.17g delta=%.17g  normal-state theta = %.17g\n", rho0, delta, 2 * rho0 * sqrt(1 - delta * delta);
}'
