#!/bin/bash
# The hybridisation of the test case: a flat band, Gamma_11(omega) = 1/2 for |omega| <= 1.
#
# 1/2 on each frequency branch is what nrgchain's band=flat uses, so the nrgchain route below sees exactly the same
# band without a tabulated input, and theta = int Gamma domega = 1. The tabulated grid is the same 201 points as in
# test/tools/mixchain/mixchain1_flat; for a flat band the linear interpolant is exact, so nothing is lost by it.
set -eu
awk 'BEGIN { for (k = -100; k <= 100; k++) printf "%.17g 0.5\n", k / 100 }' > Gamma_11-re.dat
