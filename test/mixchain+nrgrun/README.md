# End-to-end tests: a `mixchain` chain through a SIAM template

Each test discretizes a hybridisation with `mixchain`, renames the resulting chain into the coefficient sets a SIAM
template reads, lets the template's own Perl `instantiate` build `data`, and runs `nrg` on it. Mathematica is not
involved: the templates are checked in under `templates/` and read the chain from files.

```
Gamma_ij-re.dat -> mixchain -> V/E/T per element -> stage_<SYM>.pl -> the coefficient files of the template
                -> instantiate (matrix, diag, unitary) -> data -> nrg -> compared with ref/
```

Two kinds of test live here. Most run once and compare with a stored `ref/`, sharing the `runtest` script. The rest
build the same physics a second way and compare the two runs against each other, storing nothing; those carry their
own `run` script. The second kind states a relation that stays true when the numbers legitimately move, which a
refreshed reference would otherwise absorb unnoticed.

| Test | Kind | Gamma | Runs when |
| --- | --- | --- | --- |
| `mixchain_qs` | ref | a flat scalar band, QS | always |
| `mixchain_qs_vs_nrgchain_legacy`, `mixchain_qs_vs_nrgchain_rkpw` | pair | the same band through mixchain and through `nrgchain band=flat` | selected chain backends |
| `mixchain_u1` | ref | `diag(1/2, 1/4)`, field along z, U1 with `pol2x2` | `SYM_ALL` |
| `mixchain_u1_rotated` | ref | two bands of different shape, turned by pi/3, field tilted with them | `SYM_ALL` |
| `mixchain_u1_rotation` | pair | those two runs against each other, through the rotation identity | `SYM_ALL` |
| `mixchain_spsu2` | ref | an s-wave BCS bath in Nambu form, `chain_gauge=nambu`, SPSU2 | `SYM_MORE` |
| `mixchain_spsu2_nrginit_legacy`, `mixchain_spsu2_nrginit_rkpw` | pair | the same physics through `nrginit` with `bcsgap` | selected chain backends, `SYM_MORE` and Mathematica |

`TEST_CHAIN_LEGACY` and `TEST_CHAIN_RKPW` select the suffixed registrations (both default ON). Each has its own
work directory and reruns both physical routes, without using another test's outputs. Only the scalar nrgchain or
nrginit route selects the backend; mixchain's independent method is unchanged. The SC initializer route uses full
`tri=old` or `tri=rkpw` with constant `bcsgap`, not a runtime star handoff. The registrations carry `chain-legacy`
or `chain-rkpw` labels.

`mixchain_u1` leaves the mixing sets `xi3`, `xi4` vanishing, since two flat bands of unequal height share their chain
coefficients and a rotation cannot mix what is proportional to the identity. The rotated cases therefore give the two
bands different shapes, which makes `xi3` a few percent of `xi1`. A chain written as its own transpose swaps the two
mixing sets, leaves the thermodynamics untouched and reverses `SXd`, which `mixchain_u1_rotation` is what catches.

In `mixchain_spsu2` only `pair_d` carries the sign of the anomalous coefficients, which is why it is among the
operators; `mixchain_spsu2_nrginit` checks that sign and size against the established route.

## What is compared

`ref/` holds the staged coefficient files, so that a change in the chain shows up directly, and `td` and `custom`
from the run. `runtest` compares with `test/compare.pl --strict`, whose numerical tolerance is the usual 1e-5,
absolute.

## Refreshing the references

```sh
test/mixchain+nrgrun/refresh "$PWD" "$PWD/build" [TEST...]
```

runs the same pipeline with `NRG_REFRESH=1` and rewrites `ref/`. Tests of the second kind are skipped, having nothing
stored. It writes nothing else — the test lists are left
alone — and works in a scratch directory under the build tree. Read the diff before committing it.

## The templates

`templates/QS`, `templates/U1` and `templates/SPSU2` are copies of the SIAM templates from the separate
`nrgljubljana_templates` repository, generated there by `nrginit` with `GENERATE_TEMPLATE`. They are checked in for
the same reason `test/nrgspawn` checks in its own: the tests must run without Mathematica and without a second
repository. Their `mmalog` files are dropped. Regenerating them belongs upstream, not here.

## Conventions the tests depend on

- The staging of each symmetry type follows what `nrg` means by its coefficient tables: `c++/coef.hpp` for the
  order of the sets and the `sym-*-impl.hpp` macros for the terms they multiply.
- The SPSU2 case assumes the Nambu basis phase in which the anomalous part of Gamma is negative above the gap, which
  is what `make_input.sh` writes; the other phase reverses `pair_d`. `stage_SPSU2.pl` passes `E(1,2)` through as it
  stands and checks that the chain really is in the Nambu gauge.
- `channels` means the dimension of Gamma to `mixchain` and the number of physical channels to `nrg`, so each case
  has a `param` and a `param-mixchain`.
