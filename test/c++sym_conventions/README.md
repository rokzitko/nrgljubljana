# Cross-symmetry convention tests

Each group of cases solves one physical problem with several symmetry types. All cases of a group share the same
references for the symmetry-independent outputs (`custom`, `customfdm`, `tdfdm`); `td` contains symmetry-specific
columns and has a per-case reference. The groups fix the conventions described in
`docs/docs/symmetry-framework.md` ("Conventions Shared By All Symmetry Types").

| Group | Symmetry types | Problem | Reference values from |
| --- | --- | --- | --- |
| `conv_pairing_1ch_*` | SPSU2, SPU1, P, NONE | Clean SC chain, `bcsgap=0.1+0.1i` | SPSU2 |
| `conv_pairing_2ch_*` | SPSU2, SPU1, P, NONE, PP | Two clean SC chains, `bcsgap1=0.1+0.1i`, `bcsgap2=0.2-0.05i` | SPSU2 |
| `conv_pairing_2ch_equal_*` | SPSU2, SPSU2LR, SPU1LR | Two clean SC chains, `bcsgap=0.1+0.1i` | SPSU2 |
| `conv_pairing_kappa_*` | SPSU2, SPU1 | Complex hopping, complex pairing and anomalous hopping (hand-edited tables in `coefficients`) | `../c++sym_spsu2/spsu2_kappa+eps+t-complex1` |
| `conv_polarized_1ch_*` | QSZ, U1, SPU1, P, NONE | SIAM, `polarized=true`, spin-dependent particle-hole asymmetric bath | QSZ |
| `conv_polarized_2ch_*` | QSZ, U1, P, NONE, PP | Two clean chains, `polarized=true`, spin-dependent particle-hole asymmetric bath | QSZ |
| `conv_unpolarized_1ch_*`, `conv_polarized_symmetric_1ch_*` | QS, P (unpolarized); QSZ, P, NONE (`polarized=true`) | SIAM, spin-independent particle-hole asymmetric bath | QS |
| `conv_3ch_*` | QS, QSZ, U1, SPU1, P | Three clean chains, particle-hole asymmetric bath | QS |
| `conv_pairing_3ch_*` | SPSU2, SPU1 | Three clean SC chains, `bcsgap=0.1+0.1i` | SPSU2 |
| `conv_doublet_1ch_*` | NONE, P | SIAM, spectral function of the doublet operator `A_d_1` | NONE |
| `conv_substeps_2ch_*` | QS, QSZ, U1, SPU1, P, NONE | Two-channel SIAM, interleaved diagonalization (`substeps=true`), particle-hole asymmetric bath | QS |
| `conv_globalop_1ch_*` | QSZ, U1 | SIAM with a local field (`siamB.m`), global operator `SZtot` (`modeloperators.m`) | QSZ |

What the groups check:

- **Channel number.** The coefficient tables of NONE, P, U1, QSZ, and SPU1 are generated at run time for any number
  of channels (`c++/abelian_shell.hpp`); for three channels they are compared with QS and SPSU2, whose tables are
  still the Mathematica-generated ones. The NONE cases with three channels are omitted for their run time.
- **Interleaved diagonalization.** With `substeps=true` one channel is added per step; the coefficient tables are
  those of a one-channel shell, applied to the active channel.
- **Fermion sign in operator recalculation.** A fermionic operator of the previous sites acquires the sign
  (-1)^n of the occupation n of the added site. The P doublet tables generated with Mathematica lacked this sign.
- **Spin labels of the U1 hopping.** U1 stores f^dag_up in block 0 and f^dag_down in block 1 of `<||f||>`. Without
  spin polarization the previous code paired these blocks with the opposite spin of the added site, which is a
  spin relabeling of the bath sites: impurity quantities were unaffected, but the bath part of global operators such
  as `SZtot` was wrong.

- **Pairing phase.** The pairing term of the sites added in the C++ iteration must have the same phase convention as
  the pairing term of the initial sites set up by `nrginit`. With a uniform complex gap, the expectation values of the
  isospin components of the first chain site obey $\langle I_y\rangle/\langle I_x\rangle = -\mathrm{Im}\,\Delta /
  \mathrm{Re}\,\Delta$ in each channel. The SPSU2 phase convention is validated by exact diagonalization in
  `../c++sym_spsu2/*/ref.nb`. An opposite convention for the added sites changes the relative phase of neighbouring
  pairing terms and hence all observables, including the energies in `tdfdm`.
- **On-site energy.** The polarized on-site term $\sum_\sigma \zeta_\sigma (n_\sigma - 1/2)$ must agree with the
  unpolarized $\zeta (n - 1)$. A different average-occupancy subtraction shifts all energies of a shell by a constant:
  expectation values, `C_fdm`, and `S_fdm` are unchanged, but `E_fdm` and `F_fdm` (absolute energies) are not. The
  particle-hole asymmetric bath (`band=dmft` with an energy-dependent `gamma`) makes $\zeta_n \neq 0$, so the offset
  is visible; for a particle-hole symmetric bath, or with `globalh` alone, $\zeta_{n\uparrow}+\zeta_{n\downarrow}=0$
  and the test would not discriminate.

The cases need `SYM_ALL`. Run them with

```sh
ctest --test-dir build -L '^conventions$' --output-on-failure --no-tests=error
```

## Regenerating the input data

The `data` files were generated with `nrginit` from each case's `param` and Mathematica files (`customoperators.m` for
the two-channel pairing cases, which defines the isospin operators `Ix_f1`, `Iy_f1` of the second channel; the model
`siamB.m` and `modeloperators.m` for the global-operator cases). For the
`conv_pairing_kappa_*` cases, the Wilson-chain blocks are then replaced by `coefficients`. To regenerate (requires
Mathematica):

```sh
test/c++sym_conventions/generate            # all cases
test/c++sym_conventions/generate CASE...    # selected cases
```

The references are not regenerated by this script.

## Adding a symmetry type

To check a new symmetry type against an existing group, add a case directory with the group's `param` (changing only
`symtype`, and adding `bcsgap=0` where `nrginit` requires a pairing definition), list it in `simple_tests`, run
`generate CASE`, copy the shared references of the group, and add the `td` reference from the new run.
