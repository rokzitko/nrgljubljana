# QSTZ Dynamical Spin Regression

QSTZ dynamical spin susceptibility is implemented. The failure addressed here
was a coefficient-table bug, not absent triplet support. Generic large-spin
expressions contained removable radical singularities at allowed S=0 and
S=1/2 transitions. Eight QSTZ entries and one three-channel QS entry have been
replaced by finite equivalent expressions; no assertion or matrix element is
discarded to hide the problem.

In terms of target spin S, the three corrected forms are:

```text
-sqrt((2*S+3)/(2*S+1))/(S+1)        -> -2*sqrt(2)/3 at S=1/2
-sqrt(3*(S+2)/(2*S+3))/(S+1)        -> -sqrt(5/6) at S=1/2
sqrt(3)/((S+1)*sqrt((2*S+1)*(2*S+3))) -> 1 at S=0
```

The `recalc` unit test checks every structurally allowed entry in all three
QSTZ and three-channel QS triplet tables for spin multiplicities 1 through 21
against an independent GSL Clebsch-Gordan contraction. Populated runtime
endpoint tests cover both real and complex matrix elements. The Mathematica
test `nrginit_spin_recoupling` verifies that the shared coefficient simplifier
cancels these radical factors before emitting C++ expressions. Its use is
restricted to QSTZ and three-channel QS spin triplets. The generic
large-spin simplification is not by itself proof at low spin; the coefficient
and dynamical-response tests enforce those physical boundary values.

This small prepared-input regression is independent of saved solver outputs and
the SIAM/star oracles. It requires only the scientific suite's Python/NumPy
dependencies at test time, plus a solver built with `SYM_ALL=ON`. No Mathematica
license, external basis module, calculation directory, binary reference, or
HDF5 reference is needed to run it.

## Run

With `TEST_SCIENTIFIC=ON` and `SYM_ALL=ON`, configure and build as usual, then:

```sh
ctest --test-dir build -R '^scientific_qstz_spin_' \
  --output-on-failure --no-tests=error
python3 -B -m unittest discover -s test/scientific -p test_qstz_spin.py
```

The two CTest entries are `scientific_qstz_spin_T02` and
`scientific_qstz_spin_T07`. They use the build-tree executable/library,
single-thread numerical libraries, one CTest processor and a 180-second timeout
(120 seconds for each solver subprocess). Each invocation creates a fresh
`qstz-spin-*` directory under its build-tree work root. Inputs, `log`, and
`validation.json` are retained even on solver failure; a fresh `DONE` is required.
The Python self-checks are also discovered by `scientific_unit`.

Direct invocation, from the repository root on Linux:

```sh
LD_LIBRARY_PATH="$PWD/build/c++:${LD_LIBRARY_PATH:-}" \
  python3 -B test/scientific/validate_qstz_spin.py \
  --nrg build/c++/nrg --temperature 0.2 --work-root /tmp/opencode
```

## Model And Oracle

The impurity has `Tz=0`, `epsilon_d=0.15`, `V=0.3`, and no interaction or field.
It couples only to the central (`Tz=0`) bath. The `Tz=+1,-1` channels are free
spectators, but are retained in the solver Hilbert space. `Ninit=0` includes
the impurity and three bath sites (256 physical states); `Nmax=1` adds three
sites (16,384 states, 6,435 SU(2) multiplets). The seed has 51 sectors and 126
multiplets. The nontrivial local rank-one spin operator connects different
ancestor spins, exercising QSTZ triplet recalculation at the first extension.

The fixture uses `band=flat`, `Lambda=2`, `z=1`, `bandrescale=1`, and the `Z`
discretization. The actual first bath hopping is `0.54528747084262258`, not an
assumed `0.4`. The validator reads all three channels' physical hoppings/onsites
from the prepared `data` trailer and checks them against the analytic flat-band
formula, including the unused second hopping. `V` is explicit in `model.m`,
not inferred from `Gamma`. `Gamma=0.1` is only an initializer bath-normalization
parameter. The independent seed occupation spectrum checks the physical seed
energies and their ground-state offset. No NRG eigenvectors, reduced operator
matrices or coefficient/recalculation tables enter the oracle.

The oracle diagonalizes a 7x7 one-particle Hamiltonian with `numpy.linalg.eigh`.
For `p_i=|u_di|^2`, `f_i=1/(1+exp(e_i/T))`, and `omega_n=2*pi*n*T`, it evaluates

```text
Czz(i*omega_n) = (1/2) sum_ij p_i p_j (f_i-f_j)/(i*omega_n+e_i-e_j).
```

At zero frequency, degenerate pairs contribute
`-p_i*p_j*f_i*(1-f_i)/(2*T)`; at nonzero frequency they contribute zero.
The spin mean is zero, so no disconnected term remains. The sign convention is
the solver's bosonic Green function, `-<T_tau Sz(tau) Sz(0)>`, whose static value
is negative. Positive static susceptibility is `-Czz(0)`.

`sigma_d` is the standard `mtTripletOp["sigma_d", d[]]` in
`nrginit/operators.m`, not a Pauli matrix or a fermionic doublet callback.
`ireducsigma` in `nrginit/initial.m` uses `Sz`, `-Splus/sqrt(2)` and
`Sminus/sqrt(2)`. `dynamicsusceptibility_factor` in
`c++/symmetry/sym-QSTZ-impl.hpp` supplies the final spin multiplicity divided
by three. Thus the result is one Cartesian component, not `<S.S>`; there is
no fitted factor of three, four or sign. Analytic atomic and symmetric-dimer
self-checks lock down this normalization and the zero-frequency thermal term.

At `T=0.2` and `0.7`, all four rows of
`spin_FDMmats_dens_sigma_d-sigma_d.dat`, including `n=0`, are checked. Real and
imaginary components use the scientific suite's symmetric
`atol=1e-11, rtol=1e-10`; frequency tolerance is `1e-14` for both. Nonfinite,
missing, empty, wrong-shaped and numerically incorrect results fail. No
broadening, binning-to-Matsubara transform, or saved spectral reference is used.

The validator also requests `specchit=sigma_d-sigma_d` at runtime without
changing the fixture inputs. Both rows of its CHIT curve (seed and final
shell) must equal `-T*Czz(0)` evaluated independently at the reported shell
temperature with the corresponding 4x4 or 7x7 Hamiltonian. This distinguishes
the local Kubo response from an equal-time spin moment and from the fixed-T
FDM result. The `specchit` C++ unit test additionally checks free spin-1/2
and spin-1 Curie constants, a nonconserved Heisenberg-dimer local spin,
scalar-pairing preservation and real/complex spin-SU(2) dispatch.
Nonzero orbital multiplets in QST/SPSU2T test the spectator factor `2*T+1`,
including spin-changing transitions in the local-dimer response.

`keep=100000`, `keepenergy=-1`, full `dsyevd` diagonalization and explicit
`lastall=true` prevent truncation. The validator independently counts multiplets
from occupation states and requires the log to report all 256 seed states,
`nrkept=nrall=6435`, and `nrkeptmult=nrallmult=16384`, with exactly one NRG and
one DM extension. This is not a truncation/continuum/production benchmark or a
test of early-shell FDM discarded-kept terms. NumPy and NRG may share LAPACK.

## Regeneration

Only regeneration needs a licensed Wolfram kernel. It uses the bundled
`nrginit/QSTZbasis.m` by default; an explicit `--basis` argument can override
it. The model declares `qstzImpurityCharges={d[]->0}`, while the adapter uses
bath charges `{f[1]->1,f[0]->0,f[-1]->-1}` and quantum numbers `{Q,2S+1,Tz}`.
The prepared data already contains the reduced basis information and needs
no Mathematica or basis module at runtime.

```sh
python3 -B test/scientific/prepare_qstz_spin.py \
  --kernel /path/to/WolframKernel \
  --work-root /tmp/opencode
```

This uses repository `sneg.m`/`initial.m` directly in an isolated temporary
parent and child, preserving license discovery but excluding home-directory
initializer overrides. The candidate contains the exact `model.m`, `param`,
entry script, initializer log, `data`, and `provenance.json`; an explicit basis
override is also copied into the candidate. The default instead exercises
package lookup, with no calculation-local basis file.
Use `validate_qstz_spin.py --fixture <candidate>` at both temperatures to
qualify it. `--write` on the preparation command publishes only the generated
text `data` and provenance; it does not run NRG or imply qualification. Numerical
qualification, not eigenvector-phase/byte equality, judges regeneration.

The initial fixture was generated on 2026-11-18 with Mathematica 13.3, precision
100 and `mMAX=80`, source revision
`696d061438b04155b3a6279e8de0bdaa4dd092d2`. The development-only basis came from
`/project/zitko/calc/spinaron/concept/nrg_qstz_r0/QSTZbasis.m`, SHA-256
`b7f6e01b063f27511f04eb707561eac50e04055bb04a7cc80fac8eba8bc500b4`.
The manifest records initializer-source, preparation-script, basis and input
hashes. No production calculation was run.

The original unpatched build deterministically exits with status 1 in the DM
pass at `SymmetryQSTZ<double>::recalc_triplet`, asserting `isfinite(factor)` in
`recalc.hpp:186`, after reporting full retention. This is a normal failing
regression, not `WILL_FAIL` or an expected-failure exemption. Rerun the two CTest
entries after rebuilding the corrected solver. Both temperatures pass with
the corrected tables, including the static response and three nonzero
Matsubara frequencies.
