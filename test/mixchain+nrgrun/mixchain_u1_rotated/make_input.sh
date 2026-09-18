#!/bin/bash
# The hybridisation of this case, rotated in spin space by an angle:
#
#   Gamma(phi) = U Gamma(0) U^dagger,   U = [[c, -s], [s, c]],  c = cos(phi/2), s = sin(phi/2)
#
# with the two spin bands of *different shape*:
#
#   g_up(omega) = 1/2,   g_dn(omega) = 1/4 (1 - omega^2/2)      for |omega| <= 1,
#
# which is
#
#   Gamma_11 = c^2 g_up + s^2 g_dn,  Gamma_22 = s^2 g_up + c^2 g_dn,  Gamma_12 = Gamma_21 = c s (g_up - g_dn).
#
# U is the spin-1/2 rotation by phi about the y axis, so a vector quantity rotates by phi while the matrix carries
# the half angle. phi=0 gives back a diagonal Gamma, with vanishing off-diagonal files.
#
# The shapes have to differ for this case to test what it is here for. Two flat bands, however unequal their heights,
# have the same chain coefficients xi_n -- the height goes into theta alone -- so T_n is proportional to the identity
# and a rotation leaves it diagonal, putting the whole mixing into V. With an omega-dependent g_dn the two branches
# have different xi_n, the rotated T_n has off-diagonal elements, and the mixing coefficient sets xi3 and xi4 carry
# something. Both bands stay even in omega, so the on-site energies still vanish.
#
#   make_input.sh [PHI]     (radians, default pi/3; writes the four Gamma files into the working directory)
set -eu
phi=${1:-1.0471975511965976}   # pi/3, the angle the param of this case uses for the field

awk -v phi="$phi" 'BEGIN {
  c = cos(phi / 2); s = sin(phi / 2);
  for (k = -100; k <= 100; k++) {
    w = k / 100;
    gup = 0.5;
    gdn = 0.25 * (1 - w * w / 2);
    g11 = c * c * gup + s * s * gdn;
    g22 = s * s * gup + c * c * gdn;
    g12 = c * s * (gup - gdn);
    printf "%.17g %.17g\n", w, g11 > "Gamma_11-re.dat";
    printf "%.17g %.17g\n", w, g22 > "Gamma_22-re.dat";
    printf "%.17g %.17g\n", w, g12 > "Gamma_12-re.dat";
    printf "%.17g %.17g\n", w, g12 > "Gamma_21-re.dat";
  }
  printf "make_input: phi=%.17g  at omega=0: Gamma = [[%.17g, %.17g], [%.17g, %.17g]]\n",
    phi, c * c * 0.5 + s * s * 0.25, c * s * 0.25, c * s * 0.25, s * s * 0.5 + c * c * 0.25;
}'
