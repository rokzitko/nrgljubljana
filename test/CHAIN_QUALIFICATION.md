# Scalar Mapping Qualification

These tests separate reconstruction error on an identical finite star from
finite-cutoff error and from physical frontend/seed errors. They do not change
production defaults or call the other production backend as an oracle.

**Current qualification policy:** compact qualification passes for both backends.
The known extended RKPW hard-gap hopping miss documented below is a non-blocking
advisory, not a numerical fix. Successful tests neither prove a universal
`2e-12` bound nor recommend a production-default switch.

## Commands And Gates

With ordinary unit tests enabled, build and run the compact numerical suite:

```sh
cmake --build build --target scalar_chain_qualification --parallel 2
ctest --test-dir build -R '^scalar_chain_qualification_(legacy|rkpw|report_(legacy|rkpw))$' \
  --output-on-failure --no-tests=error
```

Backend registration follows `TEST_CHAIN_LEGACY` and `TEST_CHAIN_RKPW`.
`TEST_CHAIN_QUALIFICATION_EXTENDED` defaults OFF and is intentionally independent
of `TEST_LONG`. To run the extended acceptance gate:

```sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX="$HOME/nrgljubljana/" \
  -DTEST_CHAIN_QUALIFICATION_EXTENDED=ON
cmake --build build --target scalar_chain_qualification --parallel 2
ctest --test-dir build -R '^scalar_chain_qualification_extended_(report_)?(legacy|rkpw)$' \
  --output-on-failure --no-tests=error
```

Under this policy, the known RKPW miss measured on GCC 14.2/x86-64 is advisory.
No `WILL_FAIL` property inverts the result. Different platform rounding may meet
the unchanged target and produce no advisory. Disable the extended option again
for the ordinary compact suite.

Tests have `chain-qualification` and their independent `chain-legacy` or
`chain-rkpw` labels; extended numerical entries additionally have
`chain-qualification-extended`. Each backend/tier has a separate workdir under
`build/test/unit/nrgchain/`. Neither needs another backend's outputs.

The scientific finite-star frontend suite is separately opt-in through
`TEST_SCIENTIFIC`. Its extra units/temperature/solver-scaling combinations use
`TEST_SCIENTIFIC_STAR_EXTENDED`, not the numerical acceptance-gate option above.
See [asymmetric finite-star qualification](scientific/STAR.md).

## Same-Star Contract

`test/unit/nrgchain/scalar_chain_qualification.cpp` builds analytic scalar
measures on a fixed logarithmic mesh, computes their exact-form shell masses
and Z representative energies in high precision, then rounds each input
energy and amplitude once to double. The same rounded values are sent to:

- the real `nrgchain` loaded-table API;
- runtime `Tridiag<double>`;
- runtime `Tridiag<std::complex<double>>` with real input.

The tool adapter uses round-trip text precision. All three paths request a
strict prefix shorter than effective support, since legacy reconstruction
does not have an exact-terminal-zero contract. They use physical coefficients
with `bandrescale=1` and `rescalexi=false`; separate scientific tests cover
nonunit bandwidth and impurity-coupling conventions.

The oracle is independent diagonal-star Krylov construction with Rayleigh
quotients, residual norms and two full reorthogonalization passes. It uses
Boost's header-only multiprecision only in the test reference, at explicitly
specified **binary** precisions:

| Contract | Compact | Extended |
| --- | --- | --- |
| Oracle precision pair | 256 / 512 bits | 512 / 1024 bits |
| Relative hopping / local onsite oracle agreement | `<1e-35` | `<1e-35` |
| Maximum basis orthogonality defect | `<1e-50` | `<1e-50` |
| Production hopping relative-error target | `<=2e-12` | `<=2e-12` |
| Production onsite local-scale error | `<=2e-12` | `<=2e-12` |

The unchanged `2e-12` is a provisional engineering target reused from earlier
tests, not a theoretical or paper-derived bound or a demonstrated physics
requirement. Only RKPW hop index 19 of the exact extended `gap_large_lambda`
tuple below is exempt from its fatal relative-error gate. All other hops and
all onsites, finite/positive checks (including that hop), and independent
reference convergence remain fatal gates. Legacy is unchanged.

Precision agreement is checked before narrowing the oracle. An unconverged
reference fails rather than becoming the expected result. Onsite errors use
`max(abs(zeta[n]), xi[n], xi[n-1])` from the reference, omitting the preceding
link at site zero. There is no fixed bandwidth-sized absolute floor on tail
hoppings. Legacy precision is explicitly set to a conservative test budget,
`512 + ceil(count^2 log2(Lambda))` GMP bits; RKPW does not use that parameter.

The compact matrix has eight stars: a flat analytic anchor, power laws with
exponents `0.5`, `1`, and `2`, a smooth asymmetric density, a near-unity
`Lambda=1.05` case, and symmetric/asymmetric hard gaps. Positive/negative
weights include ratios up to 20:1, twists `0.05`, `0.37`, and `1`, and 16--24
requested coefficients. The extended matrix adds ten cases, reaching 100
coefficients and `Lambda=12`. The complete immutable tuples live in the
driver's `matrix()` function and are recorded in each report.

Power measures are `rho_+(E)=E^r`, `rho_-(-E)=Cminus*E^r`. Smooth measures are
`rho_+(E)=1+0.8E`, `rho_-(-E)=Cminus*(1-0.4E)`. Hard gaps use constant density
on `boundary<abs(E)<1`, obtained by an affine transformation of the flat measure.
Stars remain in outer-to-inner sign-paired order. This study does not claim
accuracy for arbitrary insertion orders or unresolved gap-edge poles.
Analytic construction isolates reconstruction and cutoff effects; it does not
measure `adapt`/FSOL tabulation error or general DOS-integration precision.
Those need separate input-resolution studies, not a larger `preccpp` value.

## Cutoff And Finite Resolvent

The cutoff study fixes `Lambda=1.05`, `z=0.37`, `r=2`, `Cminus=0.05` and a
16-coefficient prefix. Extended mode uses `mMAX=80,160,320,640`; compact mode
uses the final two. Each finite star first receives its own precision-converged
reference and production check. Reference prefixes are then compared with the
`mMAX=640` anchor, with a `1e-10` acceptance target only on the final pair.
Coarse differences are diagnostics; convergence is not assumed monotone.
Analytic omitted-weight fractions are reported separately.

On the tested build, the final-pair hopping drift is about `4.91e-19`. The
coarse `mMAX=80` prefix differs by about `7.93e-4` even though reconstruction
of that finite star is accurate. This is omitted-shell error, not a failure
of the tridiagonalizer, and illustrates why the two checks are separate.

An additional eight-pole test compares the boundary-bath resolvent with its
direct spectral sum. Both backends produce `rank-1` coefficients; the last
diagonal is closed by trace invariance, with an exact zero terminal link.
This is an explicitly identified test adapter, not a claim that production
legacy supports terminal hopping generation. The independent oracle also
reconstructs the entire finite matrix and checks that identity separately.
The scientific suite instead validates a full physical impurity Green function
and treats legacy's truncated prefix as a different physical system.

## Reports

Each numerical invocation writes
`scalar_chain_qualification_<backend>_<compact|extended>.tsv` and a companion
`.tsv.status.json` in its own workdir. Rows contain the case parameters,
precision budgets, frontend, error maxima and worst indices, omitted weight,
gap-edge distance, and gate status. Metric families are `reference_precision`,
`reconstruction`, `known_gap_hop`, `analytic_anchor`, `cutoff`, and `bath_resolvent`.

The reporting contract splits only the exempt RKPW hop into `known_gap_hop`,
with `gated=0` and status `advisory` when the target is exceeded, `pass` otherwise.
The `reconstruction` row retains all remaining hops and all onsites as fatal
checks. An exceedance prints `ADVISORY` to stdout; the TSV retains the measured
error, index, and target. The `gated` column also distinguishes coarse-cutoff
diagnostics from acceptance checks.

Always check the process exit status and completion JSON as well as TSV rows:
a case-construction or oracle exception may prevent a row from being emitted.
Completion retains `running`, `passed`, `failed`, `not-run`, and `incomplete`
and the actual execution/failure/skip counts. It adds `advisories`, counting
exceeding frontend rows, including repeated executions. `passed-with-advisories`
requires at least one advisory and completed real execution with no fatal
failures or skips; no-op runs cannot receive this status. With zero advisories,
otherwise successful completion remains `passed`.
Listing tests, an empty filter, or zero repetitions cannot report a successful
qualification. `scalar_chain_qualification_report_<backend>` checks these
contracts, including inherited `GTEST_REPEAT=0` and repeated real executions.
The opt-in `scalar_chain_qualification_extended_report_<backend>` tests also
check the known-hop rows, retained blocking checks, stdout warnings, and
advisory counts against a fresh extended matrix run.

## Known Hard-Gap Limit

The extended `gap_large_lambda` tuple is:

```text
Measure=gap, Lambda=4, z=1, r=0, Cminus=0.05, boundary=0.1, mMAX=13, count=20
```

It remains well inside the test's geometric representability checks, but the
inverse reconstruction is sensitive. On GCC 14.2/x86-64, hop index 19 is:

```text
reference: 4.3519025361456866e-6
RKPW:      4.3519025361646171e-6
relative error: approximately 4.35e-12 (target 2e-12)
absolute error: approximately 1.89e-17
```

The 512/1024-bit independent references agree to roughly `8.5e-151` in relative
hopping error. All three C++ frontends reproduce the miss, while legacy passes.
Additional scratch investigations found that one adjacent-double perturbation
of an input energy can change this hopping by about `7.56e-12` relatively.
FMA or local arithmetic rearrangements that improve this single case do not
consistently qualify its neighborhood. These diagnostics do not justify a
universal relaxed gap tolerance.

By explicit policy, this exact case remains enabled, with only the hop-19
relative-error miss treated as advisory. There is no global tolerance
relaxation, algorithm patch, higher-precision production fallback, or
production-default change. Broader accuracy claims or a default-switch
recommendation need separate qualification.
