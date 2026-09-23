#!/bin/bash
# A spin-dependent hybridisation, diagonal in spin space:
#
#   Gamma(omega) = diag(1/2, 1/4) for |omega| <= 1, channel 1 = up, channel 2 = down.
#
# Both spins see a flat band, so the chain coefficients xi and zeta are the same for the two and are those of
# nrgchain's band=flat; only the weights differ, theta_up = 1 and theta_down = 1/2. That is what makes the reference
# route below exact without a second discretization.
#
# All four -re files are required, including the vanishing off-diagonal ones. No -im file is written, so Gamma is
# real and mixchain runs in real arithmetic.
set -eu
awk 'BEGIN { for (k = -100; k <= 100; k++) printf "%.17g 0.5\n", k / 100 }'  > Gamma_11-re.dat
awk 'BEGIN { for (k = -100; k <= 100; k++) printf "%.17g 0.25\n", k / 100 }' > Gamma_22-re.dat
awk 'BEGIN { for (k = -100; k <= 100; k++) printf "%.17g 0\n", k / 100 }'    > Gamma_12-re.dat
cp Gamma_12-re.dat Gamma_21-re.dat
