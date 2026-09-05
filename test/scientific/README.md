# Scientific Validation

This opt-in suite compares untruncated, finite-chain single-impurity Anderson
model (SIAM) calculations with a separately constructed exact-diagonalization
(ED) reference. It is distinct from the saved-output regression comparisons
described in [`test/README.md`](../README.md).

## Running Tests

The suite requires Python >=3.10 and NumPy >=1.26,<3, in addition to the normal
[source-build dependencies](../../README.md#advanced-compilation-from-source).
`TEST_SCIENTIFIC` defaults to `OFF`; enabling it requires `Build_Tests=ON` and a
top-level build. Leaving it off adds no scientific-test Python or NumPy
dependency to the default build.

From the repository root, create an environment outside the checkout, install
the requirements, and configure the existing `build` tree with that exact
interpreter:

```sh
python3 -m venv "$HOME/nrg-scientific-venv"
"$HOME/nrg-scientific-venv/bin/python" -m pip install \
  -r test/scientific/requirements.txt
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/nrgljubljana" \
  -DBuild_Tests=ON -DTEST_SCIENTIFIC=ON \
  -DPython3_EXECUTABLE="$HOME/nrg-scientific-venv/bin/python"
cmake --build build --target nrg --parallel
ctest --test-dir build -L '^scientific$' --output-on-failure --no-tests=error
```

CMake checks NumPy using `Python3_EXECUTABLE`, and CTest uses the same
interpreter. Installing NumPy into a different Python environment will not
satisfy that check. CTest runs the build-tree `nrg` target and prefers the
build-tree library; installation is not required. For a multi-config generator,
select the same configuration with `cmake --build build --config Release` and
`ctest --test-dir build -C Release`, retaining the other options above.

All 13 CTest entries have the `scientific` label: one Python unit-test entry and
12 NRG comparisons. The latter are three fixtures times two temperatures
(`T=0.05`, `0.2`) times two execution modes (`rescaled`, `absolute`). Their exact
names are:

```text
scientific_unit
scientific_siam_qsz_T005_rescaled
scientific_siam_qsz_T005_absolute
scientific_siam_qsz_T02_rescaled
scientific_siam_qsz_T02_absolute
scientific_siam_qs_T005_rescaled
scientific_siam_qs_T005_absolute
scientific_siam_qs_T02_rescaled
scientific_siam_qs_T02_absolute
scientific_siam_u0_qsz_T005_rescaled
scientific_siam_u0_qsz_T005_absolute
scientific_siam_u0_qsz_T02_rescaled
scientific_siam_u0_qsz_T02_absolute
```

A focused run, or a direct invocation with an explicit work root, is:

```sh
ctest --test-dir build \
  -R '^scientific_siam_qsz_T005_rescaled$' --output-on-failure --no-tests=error
"$HOME/nrg-scientific-venv/bin/python" -B test/scientific/validate_siam.py \
  --fixture test/scientific/fixtures/siam_qsz \
  --nrg build/c++/nrg --temperature 0.05 --mode rescaled \
  --work-root build/scientific-validation
```

The direct commands here use the single-config executable path; adjust it for
multi-config builds. CTest also fixes numerical-library thread counts to one
and uses `LC_ALL=C`. The adapter applies these settings to each NRG subprocess.

## Reference Boundary

| File | Responsibility |
| --- | --- |
| [`ed_siam.py`](ed_siam.py) | Build the Hamiltonian and fermion operators from physical parameters; solve all number sectors; evaluate thermodynamics, expectations, and Lehmann Green functions. |
| [`test_ed.py`](test_ed.py) | Check fermion anticommutation, sector dimensions and eigenpairs, atomic and noninteracting limits, thermal limits, spectral moments, and invalid inputs. |
| [`test_validation.py`](test_validation.py) | Check parsers, symmetry multiplicities, fixture integrity, isolation, and failure handling without running real NRG. |
| [`validate_siam.py`](validate_siam.py) | Translate a fixture model into NRG parameters, inspect its seed and chain coefficients, run NRG, and compare parsed outputs with ED. |
| [`prepare_siam.py`](prepare_siam.py) | Generate candidate NRG inputs with the repository initializer and qualify them before writing or checking fixtures. |

The independent reference is the Hamiltonian/operator construction and the
full finite-temperature calculation from the model specification. ED does not
read NRG `data`, eigenvalues, eigenvectors, operator matrices, or saved reference
outputs to construct its answer. The adapter reads those input/output formats
only on the NRG side of the comparison.

ED uses `numpy.linalg.eigh`. NumPy and NRG may share a LAPACK implementation or
even the same eigensolver algorithm, so this is not a claim of an independent
diagonalization algorithm or independent numerical libraries.

## Hamiltonian And API

`AndersonModel(epsilon_d, U, V, zeta, t, B=0.0)` describes one impurity `d` and
bath sites `f0`, `f1`, ... . All parameters must be finite and real. Hopping is
spin-conserving, with no spin-flip, pairing, or complex hopping terms. There must
be 1-5 bath sites, at most 6 spatial orbitals including the impurity, and
`len(t) == len(zeta) - 1`.

The physical Hamiltonian is:

```text
H = epsilon_d * (n_d_up + n_d_down) + U * n_d_up * n_d_down
  + (B/2) * (n_d_up - n_d_down)
  + sum_j zeta_j * (n_j_up + n_j_down - 1)
  + V * sum_sigma (d_sigma^dagger f0_sigma + h.c.)
  + sum_j t_j * sum_sigma (fj_sigma^dagger f(j+1)_sigma + h.c.)
```

In particular, the impurity uses the unshifted `epsilon_d * n_d + U * n_d_ud`
convention, the bath onsite terms are centered at one electron, and positive
`B` raises the up-spin impurity level by `B/2`. Dropping the bath constant
`-sum_j zeta_j` changes absolute energies and free energies, even though it
does not change excitation gaps.

Spin orbitals are ordered `d_up, d_down, f0_up, f0_down, f1_up, ...`, starting at
bit zero of the integer Fock basis. For `L` spatial orbitals, the full space has
`4**L` states. `solve(model)` diagonalizes every `(N_up, N_down)` sector, whose
dimension is `binomial(L, N_up) * binomial(L, N_down)`, and retains every state.
The returned `AndersonSolution.energies` and `.vectors` are grouped by sector,
not globally sorted; `.sector_energies` exposes that grouping and
`.energies.min()` gives the ground energy.

The ensemble includes all particle-number sectors at chemical potential zero,
with `k_B=1` and finite, strictly positive `T`. It is not restricted to half
filling or a fixed particle number. `thermodynamics(T)` returns `E`, `F`, `S`,
and `C` for the entire finite chain, not impurity quantities obtained by
subtracting a free bath. `expectations(T)` returns `I`, `n_d`, and `n_d_ud`.
`greens(T, frequencies)` returns an array of shape `(2, len(frequencies))`,
ordered up then down, from the full thermal Lehmann sum. Supply actual complex
frequency arguments away from poles: it adds neither a factor of `i` nor
broadening implicitly.

The solver can be used without CMake, NRG, fixtures, or Mathematica. With the
Python environment above, this example runs from the repository root:

```sh
PYTHONPATH=test/scientific "$HOME/nrg-scientific-venv/bin/python" -B - <<'PY'
import numpy as np
from ed_siam import AndersonModel, solve

model = AndersonModel(-0.41, 0.93, -0.27, (0.14, -0.22), (0.31,), B=0.19)
solution = solve(model)
T = 0.2
omega = (2 * np.arange(64) + 1) * np.pi * T
print(solution.thermodynamics(T))
print(solution.expectations(T))
print(solution.greens(T, 1j * omega)[:, 0])
PY
```

Standalone self-checks likewise need only Python and NumPy:

```sh
"$HOME/nrg-scientific-venv/bin/python" -B -m unittest discover \
  -s test/scientific -p 'test_*.py'
```

## Fixture Models

The NRG adapter has a narrower, explicit model source schema than the general
`AndersonModel` API. Each fixture's `model.json` must contain exactly the
top-level fields `symmetry` and `model`. `symmetry` is `QS` or `QSZ`; `model`
must contain exactly these fields:

| Field | Meaning and constraint |
| --- | --- |
| `epsilon_d`, `U`, `B` | Finite real impurity parameters in the Hamiltonian above; `B` is mandatory in the JSON and must be zero for `QS`. |
| `Gamma` | Constant hybridization width, finite and nonnegative. |
| `D` | Bath half-bandwidth, finite and positive. |
| `Lambda` | Logarithmic discretization parameter, finite and greater than one. |
| `z` | Discretization twist; only `z=1` is supported. |
| `bath_sites` | Integer number of bath sites; 2-5 for the NRG adapter, which requires an extension beyond the `d+f0` seed. |

`flat_band_model(parameters)` constructs the reference chain from these
physical parameters rather than from the fixture's numerical coefficient
trailer. As a standalone helper it also permits one bath site and defaults
an omitted `B` to zero. For a flat band with `discretization=Z` at `z=1`, the
analytic physical coefficients are:

```text
zeta_j = 0
V = sqrt(2 * D * Gamma / pi)
t_n = D * (1 - Lambda^(-1)) / log(Lambda) * Lambda^(-n/2)
      * (1 - Lambda^(-(n+1)))
      / sqrt((1 - Lambda^(-(2*n+1))) * (1 - Lambda^(-(2*n+3))))
```

Here `log` is the natural logarithm, `n >= 0`, and `t_n` connects `f_n` to
`f_(n+1)`. `wilson_hopping(n, Lambda, D=1.0)` implements this formula without
iteration-dependent rescaling. This is the `Z` convention, not the usual
Wilson prefactor `(1 + Lambda^(-1))/2` from a different discretization.

The registered baseline uses `D=1`, `Lambda=2`, `z=1`, `B=0`, and three bath
sites. The JSON files remain the source of the full parameter sets:

| Fixture | Case |
| --- | --- |
| [`siam_qsz`](fixtures/siam_qsz/model.json) | Interacting, particle-hole-asymmetric SIAM with charge and spin-projection symmetry. |
| [`siam_qs`](fixtures/siam_qs/model.json) | The same interacting model with charge and SU(2) spin symmetry. |
| [`siam_u0_qsz`](fixtures/siam_u0_qsz/model.json) | Noninteracting (`U=0`) counterpart with charge and spin-projection symmetry. |

### Prefixes And Normalization

The initializer includes `d+f0` in the seed. Thus `Ninit=0`, `Nmax=2` means
the final system is `d+f0+f1+f2`, not an impurity with only two bath sites:

| Compared prefix | Physical state count |
| --- | ---: |
| Seed `d+f0` | 16 |
| First extension `d+f0+f1` | 64 |
| Final extension `d+f0+f1+f2` | 256 |

With `calc0=true`, the seed and first extension both have report index `N=0`;
the final extension has `N=1`. The validator checks this report sequence and
every prefix against the corresponding independently constructed ED chain.
It uses full diagonalization (`diag=dsyevd`, `diagratio=1`) and a keep limit
large enough for the full space (`keep=256` here). More importantly, it checks
every subspace's expected dimension and requires `kept == total`, so a run
cannot silently truncate states and still pass.

For each prefix with `L` spatial orbitals, charge is
`Q = N_up + N_down - L`. `QS` labels reduced multiplets by `(Q, SS)` with
`SS=2*S+1`; each multiplet energy is expanded over its `2*S+1` spin projections
before comparison with ED. `QSZ` labels physical sectors by `(Q, SSZ)` with
`SSZ=2*Sz+1=N_up-N_down+1`, **not** `2*Sz`. Its second label can be zero or
negative. The state counts above include all spin multiplicities.

Green functions are normalized per spin, with unit zeroth spectral moment and
`G_sigma(i*omega) ~ 1/(i*omega)` at large frequency. In `QS`, the single
`spec_FDMmats_dens_A_d-A_d.dat` file represents either identical spin component,
not a spin sum. In `QSZ`, the `-u` and `-d` suffixes select the up/down components
and both are compared. There is no extra factor of two or fitted normalization.

## Comparisons And Artifacts

The adapter checks the real legacy `#!9` seed in physical units (`# SCALE 1`,
`data_has_rescaled_energies=false`), including its ground-energy offset and
sector spectra. It checks every hopping and onsite coefficient in the chain
trailer against the analytic formula, including the final trailer hopping
that is not needed in the last finite-chain Hamiltonian. It then checks
absolute ground energies, all sector energies and excitation gaps at every
prefix, the no-truncation subspace counts, and final-chain observables:

- `customfdm`: identity, impurity occupancy, and impurity double occupancy.
- `tdfdm`: full-chain energy, heat capacity, free energy, and entropy.
- `spec_FDMmats_dens_A_d-A_d*.dat`: real and imaginary parts of the Green
  function at 64 positive fermionic Matsubara frequencies
  `omega_n=(2*n+1)*pi*T`, `n=0,...,63`.

`parameter_text()` in `validate_siam.py` is the source of the complete runtime
parameters. It selects direct FDM Matsubara evaluation (`fdmmats=true`) and
FDM expectations, with `fdm_cutoff=0`. No real-frequency binning, broadening,
or subsequent spectral-to-Matsubara transform enters this comparison.
`rescaled` and `absolute` select `absolute=false` and `absolute=true` execution,
respectively; both are compared in the same physical units, not by matching
one NRG run to the other.

Current numerical comparisons use the symmetric rule
`abs(actual - expected) <= atol + rtol * max(abs(actual), abs(expected))`:

| Quantity | `atol` | `rtol` |
| --- | ---: | ---: |
| Energies, gaps, thermodynamics, expectations, Green-function components | `1e-11` | `1e-10` |
| Each chain hopping and onsite coefficient | `5e-15` | `5e-15` |
| Matsubara frequency grid | `1e-14` | `1e-14` |

Shapes, sector sets, dimensions, and required table columns are checked
explicitly; empty comparisons and nonfinite numerical values fail. These
tolerances belong to the scientific adapter, not to `compare.pl` or `.tol`
regression sidecars.

Each invocation creates a fresh `run-*` directory beneath its `--work-root`
and requires a fresh `DONE` marker. CTest work roots are under
`build/test/scientific/runs/<config>/<test-name>/`. Successful and failed run
directories are preserved: inspect `param`, `data`, NRG stdout/stderr in `log`,
and the generated outputs when present. `validation.json` records status,
errors, the executable, NumPy version and build configuration, and per-quantity comparison counts
and maximum absolute errors. Initial schema or provenance failures can occur
before a run directory is created; `log` exists once NRG has been launched.
The validator prints the run path; use CTest's `-V` to see it on successful
runs as well.

## Fixture Preparation

Normal scientific tests consume prepared `data` and `provenance.json` files
and need no Mathematica kernel or license. Regeneration is a separate,
explicit maintainer action requiring a licensed Wolfram kernel and a built
`nrg`; it is not part of default configuration, CTest, or license-free CI.
From the repository root, replace the kernel path below with your executable:

```sh
"$HOME/nrg-scientific-venv/bin/python" -B test/scientific/prepare_siam.py \
  --fixture test/scientific/fixtures/siam_qsz \
  --kernel /path/to/WolframKernel \
  --nrg build/c++/nrg \
  --work-root build/scientific-preparation \
  --check
```

Choose exactly one action: `--check` as above, or replace it with `--write` to
update that fixture's `data` and `provenance.json` only after all four
temperature/mode validations pass. Repeat with the other fixture paths when
regenerating the complete suite. Candidate directories, `initialization.log`,
and validation runs are preserved under a fresh `prepare-*` directory.

`--source` defaults to the repository containing the script. Preparation calls
that source tree's `nrginit/sneg.m` and `nrginit/initial.m` directly, sets
`PACKAGEPATH={NRGDIR}` and `SYMTYPE="runtime"`, and resolves the actual symmetry
from the generated parameters. It does not use an installed `nrginit` wrapper
or home-directory overrides such as `~/sneg.m`. Both the candidate directory
and its parent are fresh, preventing optional modules in the caller's working
directory from overriding repository modules. `HOME` is not changed, so normal
kernel license discovery remains available.

Provenance records the model-case and `data` SHA-256 hashes, source revision,
aggregate initializer-source hash, preparation-script hash, kernel version,
`mMAX`, `prec`, and the full generation parameter text (`parameters`). The
generation parameters use `T=0.05` and `rescaled` mode; the candidate is then
validated in both modes at both temperatures.

`--check` regenerates a new candidate and revalidates it against ED. It also
inspects the stored fixture's seed and coefficients and verifies its case/data
hashes, without modifying source files. It does **not** require byte equality
between the candidate and stored `data`, nor does it compare their operator
matrices or eigenvector phases. Run the scientific CTest suite to validate the
stored fixtures end to end. Provenance hashes detect changed stored inputs;
they are not a byte-for-byte regeneration criterion.

### Convergence Checks

Preparation defaults to `mMAX=80` (`--mmax 80`) and `prec=1000`
(`--precision 1000`). These control the discretization/Lanczos sum cutoff and
initializer arithmetic precision, not the number of retained many-body states.
A higher-cutoff, higher-precision check is:

```sh
"$HOME/nrg-scientific-venv/bin/python" -B test/scientific/prepare_siam.py \
  --fixture test/scientific/fixtures/siam_qsz \
  --kernel /path/to/WolframKernel \
  --nrg build/c++/nrg \
  --work-root build/scientific-preparation \
  --mmax 120 --precision 1500 --check
```

For a convergence study, first vary one control at a time: use
`--mmax 120 --precision 1000`, then `--mmax 80 --precision 1500`, before the
combined check above. Repeat for each fixture and retain the provenance and
validation reports to compare numerical errors.

### Initial Qualification

On 2026-09-05, all 12 numerical cases and 50 Python self/contract checks passed
on Linux x86-64 with GCC 14.2, MKL 2024.0 ILP64, Python 3.13.1, and NumPy 2.3.1.
The 13-entry scientific CTest suite took about five seconds with two concurrent
tests. The largest absolute discrepancy among the checked numerical values was
approximately `3.1e-14`; each quantity was judged against its own stated
tolerance, without fitted offsets or normalization.

The fixtures were prepared with Wolfram Mathematica 13.3. The primary QSZ case
also passed all four temperature/mode checks after separate regenerations at
`mMAX=120, prec=1000` and `mMAX=80, prec=1500`. These checks address the short
chain's coefficient preparation, not continuum convergence. The OpenBLAS and
MKL LP64 configurations enabled in GitHub CI were not run during this local
qualification.

## Validation Limits

This is a finite-chain correctness test, not a validation of general production
NRG accuracy. In particular:

- No states are truncated. The tests do not qualify truncation errors or the
  FDM discarded-kept/kept-discarded (`DK`/`KD`) contributions that require
  discarded states at earlier shells.
- Three bath sites at `D=1`, one `Lambda`, and `z=1` do not establish continuum,
  chain-length, discretization, or z-averaging convergence. The schema's wider
  parameter range is not a qualification of every accepted model.
- The NRG fixtures use `B=0` and a flat band. General bath onsite terms, signed
  hoppings, and finite fields have standalone ED checks, not the same complete
  NRG fixture coverage. Complex, spin-mixing, and superconducting models are
  outside this reference's scope.
- Direct Matsubara checks do not validate real-frequency binning, broadening,
  analytic continuation, or the other spectral algorithms.
- The adapter uses serial CPU `dsyevd` diagonalization and BLAS multiplication.
  Passing it does not qualify GPU/CUDA, distributed, or alternative solver
  paths, even if those features were enabled when building NRG.
- NumPy-backed ED can share LAPACK with NRG. The independent physical
  construction and analytic-limit self-checks do not eliminate all common
  numerical-library failure modes.
