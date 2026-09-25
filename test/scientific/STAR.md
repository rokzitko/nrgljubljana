# Asymmetric Finite-Star Qualification

This independent suite lives entirely under `test/scientific` and is enabled
by `TEST_SCIENTIFIC`. It needs only the existing optional Python/NumPy
dependencies for license-free runs. It does not use flat-SIAM reference outputs,
another backend's results, or any `mix*` inputs or templates.

## Coverage And Commands

For each enabled member of `CHAIN_TEST_BACKENDS`, the compact suite registers:

| Test suffix before `_legacy` or `_rkpw` | Candidate path |
| --- | --- |
| `scientific_star_nrgchain` | `nrgchain l`, then consume its coefficients with that backend's prepared seed |
| `scientific_star_runtime` | Prepared seed plus `T` block, actual solver `tri=cpp` reconstruction |
| `scientific_star_instantiate` | Full `instantiate`, including chain generation, seed diagonalization, operator transforms and serialization |
| `scientific_star_nrginit` | Fresh source-only Mathematica initialization with the test-only finite-star hook, then solver |

The first three require no Mathematica installation/license. The last is
registered only with `Mathematica_FOUND AND SYM_ALL`. All tests have labels
`scientific`, `chain-qualification`, and `chain-legacy` or `chain-rkpw`.
Fresh initialization also has `chain-generation`. There are no cross-backend
dependencies. Parent options `TEST_CHAIN_LEGACY` and `TEST_CHAIN_RKPW` determine
`CHAIN_TEST_BACKENDS`; the suite does not change their defaults.

`SCIENTIFIC_TEST_ENV` passes that list to Python as comma-separated
`NRG_TEST_CHAIN_BACKENDS`. `scientific_unit` reads/copies only enabled backend
seed directories; fixture contracts for disabled backends are skipped, while
mathematical oracle tests remain backend-neutral. An empty list is valid and
skips all backend fixture contracts. Standalone discovery defaults to both
backends unless the environment variable is set explicitly, for example
`NRG_TEST_CHAIN_BACKENDS=rkpw python3 -B -m unittest discover -s test/scientific`.
Discovery contracts run the complete Python suite in isolated copies with
disabled backend asset directories absent, including the empty-list case.

With both backends enabled, the compact suite adds six license-free tests and,
when available, two licensed tests. Its primary case is `D=2`, `T=0.05`,
`absolute=false`. The scientific-subdirectory option
`TEST_SCIENTIFIC_STAR_EXTENDED=ON` adds the other seven combinations of
`D={1,2}`, `T={0.05,0.2}`, and `absolute={false,true}` per frontend/backend.
That gives 48 license-free tests, or 64 including fresh initialization. The
additional entries have a `_D1_T005_absolute`-style name component and the
`scientific-extended` label. The option defaults OFF.
It is independent of `TEST_CHAIN_QUALIFICATION_EXTENDED` (also OFF), which
controls a different diagnostic domain outside this scientific frontend suite.

After configuring the desired scientific options and building the affected
targets, run from the repository root:

```sh
cmake --build build --target nrg nrgchain instantiate --parallel 2
ctest --test-dir build -R '^scientific_star_' --output-on-failure --no-tests=error
ctest --test-dir build -R '^scientific_star_.*_rkpw$' \
  -LE '^chain-generation$' --output-on-failure --no-tests=error
```

Each entry has its own build-tree work root and uses `SCIENTIFIC_TEST_ENV`,
including the build-tree library and single-thread numerical libraries.
Direct invocation must likewise prefer the build-tree library, for example
on Linux:

```sh
LD_LIBRARY_PATH="$PWD/build/c++:${LD_LIBRARY_PATH:-}" \
  python3 -B test/scientific/validate_star.py \
  --backend rkpw --lane instantiate --nrg build/c++/nrg \
  --instantiate build/tools/instantiate --D 2 --temperature 0.2 --mode absolute \
  --work-root build/star-validation
```

Use `--nrgchain build/tools/nrgchain` for the `nrgchain` lane, no extra
executable for `runtime`, or `--kernel /path/to/WolframKernel` for `nrginit`.
Do not run in source fixtures. Every invocation creates a fresh owned parent
and `run` directory; successful and failed artifacts remain there. Initial
fixture-integrity failures occur before creating work or launching tools.

## Physical Reference

`fixtures/asymmetric_star/star.json` is the sole physical input:

```text
E = {0.9, -0.7, 0.4, -0.2}       order: (+0,-0,+1,-1)
v = {0.06, 0.12, 0.12, 0.24}
epsilon_d = -0.17, U = B = 0
V = ||v|| = 0.3, sum(v^2) = 0.09
normalized amplitudes = {0.2, 0.4, 0.4, 0.8}
```

`finite_star.py` constructs the reference Jacobi matrix in physical units by
orthogonalizing the fixed diagonal star multiplication operator with two full
reorthogonalization passes. Neither this function nor its Green-function
oracle accepts candidate coefficients. The known first onsite is `-0.14` and
first hopping is `sqrt(0.1424)`. In particular, `xi[0]` is the **f0-f1 hopping**,
not the impurity hybridization `V`.

For the complete four-site bath, the primary Green-function reference is the
direct pole formula, not a candidate-chain resolvent:

```text
Delta(s) = sum_k v_k^2 / (s - E_k)
G_d(s) = 1 / (s - epsilon_d - Delta(s))
```

The legacy algorithms cannot represent exact terminal finite support. Legacy
lanes therefore request `Nmax=2`: three bath sites, a strict prefix of the
independently constructed four-site Jacobi matrix. The **reference prefix**
Green function uses that mathematical prefix's spectral measure. It is neither
the full-star Green function nor the result of dropping one physical pole.
RKPW requests `Nmax=3`: all four sites, with exact terminal `xi[3]=0`.
Both backends must supply `Nmax+1` hoppings and onsites, including the unused
extra hopping at the end of a strict prefix.

All lanes start at `Ninit=1`, the impurity plus two bath sites (64 physical
states). Final legacy and RKPW systems have 256 and 1024 states respectively.
`keep=1024`, full diagonalization and explicit `kept == total` sector checks
exclude many-body truncation. The seed report has index 0; extensions start at
1. Every prefix spectrum is compared with `ed_siam.solve`, including absolute
energies and gaps. The reference includes the bath constant `-sum(zeta)` from
the centered onsite terms `zeta*(n-1)`.

Final thermodynamics, identity and impurity occupancy are checked with ED.
Both real and imaginary components of the per-spin Green function are checked
at 64 Matsubara frequencies. The asymmetric real component is nonzero. Oracle
self-tests additionally compare the full star with direct single-particle
matrix inversion, full many-body occupation spectra, and finite-temperature
Lehmann Green functions of each prefix.

## Frontend Units

- Saved tables have positive energy magnitudes `abs(E)/D` and normalized `du`.
- `theta.dat` is the **physical** value `pi*0.09`, independent of `D`.
- The initializer hook sets `thetaCh=theta`, `df=thetaCh*du_pos^2`, and
  `dfminus=thetaCh*du_neg^2`. `Gamma=1/D` makes
  `theta0=D*Gamma*thetaCh=theta`, hence `sqrt(theta0/pi)=V`.
- `import_star.m` runs through `hook_pre_lanczosinit`; it overrides only star
  inputs. Actual normalization, reconstruction, model construction, seed
  diagonalization, operator generation and serialization still execute.
- Imported doubles are promoted exactly with `SetPrecision[..., Infinity]`
  before setting the working precision to 80 digits. No tolerance-based
  rationalization changes their values. Every fresh initialization first runs
  `check_star_import.m` against the real hook: `0.5` and its adjacent double
  `0.5000000000000001` must remain distinct and exactly equal to the imported
  binary values. The check uses an isolated `import-check/` directory and does
  not modify the physical input tables.
- `band=flat` is only an upstream placeholder. Saved-table mode and the hook
  replace the star; these tests do not qualify density discretization.
- Full instantiation explicitly selects non-`cpp` `tri`,
  `nrgchain_tables_load=true`, and `rescalexi=false`. Its data must contain
  only the `z` trailer, not a deferred `T` block.
- Prepared runtime seeds contain physical `z` entries through `Ninit`, plus
  the normalized `T` block. At `D=1` only the normalized star is re-expressed;
  the physical seed stays unchanged. Runtime stdout coefficients are printed
  **before** multiplication by `D`; the validator accounts for that factor.

All available coefficient tables are checked: `data/z`, `xi.dat`, `zeta.dat`,
and the instantiator's additional `xi1.dat`/`zeta1.dat` and theta copies, or
runtime stdout as appropriate. Coefficients use `atol=rtol=5e-15`; physical
outputs retain the existing scientific tolerances `atol=1e-11, rtol=1e-10`.
No fitted energy offsets or normalization factors are used.

## Fixture Generation

The new assets comprise one symbolic QS template (23 text files: `data.in`,
10 Hamiltonians, 12 impurity-operator matrices) and two separately generated
backend seeds. Each set has its own provenance. They are solver **inputs**,
not golden binary/HDF5 or physical-output references. No existing flat fixture
or provenance was regenerated.

The template uses ordinary `GENERATE_TEMPLATE`, `Ninit=1`, `SCALE=1`,
`data_has_rescaled_energies=false` and `I A_d n_d`. It has symbolic `eps` and
`U`; it does not embed a different interacting impurity Hamiltonian. In this
mode `tri=old` avoids numerical export guards, but `GENERATE_TEMPLATE` disables
chain reconstruction entirely. The template provenance's `backend=rkpw`
records the backend selected for its subsequent qualification, not a legacy
or RKPW reconstruction during symbolic generation.

`prepare_star.py` uses only repository `sneg.m`/`initial.m` entry points,
`PACKAGEPATH={NRGDIR}`, and isolated directories, avoiding user model/operator
overrides without changing license discovery. Provenance records kernel
version, revision, all initializer `.m` source hashes (excluding `mix*` paths),
the aggregate source hash, preparation/hook/entrypoint/import-check hashes,
physical-input file hashes, star specification hashes, precision and complete
parameters. The source aggregate hashes the canonical sorted JSON source-hash
map (`case_digest`), so it can be verified without reading historical sources.
Manifest checks fail on missing, changed, or undeclared template assets.
The verifier requires every claimed metadata field, validates source scope,
paths and SHA-256 formats, checks the aggregate against its map, and matches
the exact five input-file hashes to the canonical saved-star text at `D=2`.
It also parses the generation parameter record: seeds require `tri=cpp` and
the matching backend method; templates require the no-reconstruction
`GENERATE_TEMPLATE` configuration. Both require the fixed `Ninit=1`, `mMAX=1`,
`prec=80`, physical seed scaling, `bandrescale=2`, `Gamma=0.5`, and associated
model/unit controls. Missing or inconsistent records fail before tool launch.
Current source hashes need not match historical provenance: regeneration is
qualified by physics, not by eigenvector phases or byte equality.

To regenerate and qualify without changing source fixtures:

```sh
LD_LIBRARY_PATH="$PWD/build/c++:${LD_LIBRARY_PATH:-}" \
  python3 -B test/scientific/prepare_star.py \
  --backend rkpw --kernel /path/to/WolframKernel \
  --nrg build/c++/nrg --instantiate build/tools/instantiate \
  --work-root build/star-preparation --check
```

Add `--template` to generate the symbolic template instead of a `tri=cpp` seed;
choose `--backend legacy` for an independent legacy seed. Every candidate is
qualified at both temperatures in both execution modes before publication.
`--write` replaces `--check` to publish qualified star assets. It can update an
existing asset directory only if its file set exactly matches the generated
manifest; it never removes unexpected files. To prepare a replacement for
review without updating the fixture, use `--check` or a separate `--fixture`
directory containing the physical `star.json`.

`validation.json` records the selected backend/frontend, physical reference
kind, bath prefix length, omitted bath sites, state counts, numerical checks
and failure diagnostics. `generation.json` and `initialization.log` preserve
fresh initialization provenance. No artifacts are written into source fixtures
except by the explicit preparation action.

## Limits

This is a small real scalar finite-system qualification, not a continuum,
truncation-error, interacting nonflat, complex, matrix-chain or broadening
benchmark. Extended coverage changes units, temperature and solver scaling,
not the physical star. NumPy and the solver may share LAPACK, while the
physical construction and direct-pole Green-function oracle remain independent.
