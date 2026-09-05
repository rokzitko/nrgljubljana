# Regression comparison tools

This directory contains the regression-test runners, reference results, and
comparison tools used by NRG Ljubljana. The central comparison path is:

```text
test runner
  -> clean-physical-output.pl
  -> program under test
  -> compare.pl
       -> mycomp.pl for text references
       -> built-in semantic validators for binary and structural outputs
```

The comparison tools are intended for numerical regression testing. They do
not perform byte-for-byte comparison unless a test runner explicitly invokes
`cmp` or `diff` instead.

## Scientific Validation

The separate [`scientific/` suite](scientific/README.md) compares untruncated
finite-chain SIAM calculations with a model-driven NumPy ED reference, rather
than saved NRG outputs. It is opt-in with `-DTEST_SCIENTIFIC=ON` and requires
Python >=3.10 and NumPy >=1.26,<3 in CMake's selected interpreter. Prepared
fixtures need no Mathematica license; regeneration is a separate action.

After following the suite's setup instructions, run its unit-test entry and
12 temperature/mode comparisons with:

```sh
ctest --test-dir build -L '^scientific$' --output-on-failure --no-tests=error
```

The suite guide documents model conventions, exact test names, tolerances,
preserved run artifacts, and validation limits. The regression policies below
do not apply to its ED comparisons.

## Components

| File | Purpose |
| --- | --- |
| `mycomp.pl` | Compare two nonempty text files field by field, with absolute and relative numerical tolerances. |
| `compare.pl` | Compare a reference directory with an actual-results directory and dispatch semantic validators. |
| `mycomp-ignoresigns.pl` | Run `mycomp.pl` with sign-insensitive comparison enabled. |
| `compare-ignoresigns.pl` | Run `compare.pl` with sign-insensitive comparison enabled. |
| `PhysicalOutput.pm` | Define which result filenames are physical or spectral outputs. |
| `clean-physical-output.pl` | Remove stale classified outputs and `DONE` before a regression run. |
| `scripts/comparators.t` | Exercise the comparison tools as an executable contract test. |

## Regression workflow

The main NRG test runners use the following lifecycle:

1. Remove stale physical outputs and the previous `DONE` marker.
2. Run the executable in the test's build-directory working tree.
3. Require the executable to produce a fresh `DONE` marker.
4. Compare generated files with the source test's `ref/` directory using
   `compare.pl --strict`.

Tool tests usually invoke `compare.pl` without `--strict`. Initializer-only
tests use `compare-ignoresigns.pl` where basis-vector gauge choices can change
signs. Some specialized tests use direct `mycomp.pl`, `cmp`, or `diff` calls.

## Text comparison

### Invocation

```sh
test/mycomp.pl [OPTIONS] EXPECTED ACTUAL
```

| Option | Meaning |
| --- | --- |
| `--atol VALUE` | Set the absolute tolerance for ordinary fields. Default: `1e-12`. |
| `--rtol VALUE` | Set the relative tolerance for ordinary fields. Default: `1e-5`. |
| `--first-atol VALUE` | Set a separate absolute tolerance for the first whitespace-separated field. By default it inherits `--atol`. |
| `--first-rtol VALUE` | Set a separate relative tolerance for the first whitespace-separated field. By default it inherits `--rtol`. |
| `--check-comments` | Compare comment records instead of ignoring them. |
| `--ignore-signs` | Permit the sign relaxations described below. |

All tolerance values must be finite, representable, and nonnegative.

### Records and fields

Blank lines are ignored. A line whose first non-whitespace character is `#`
is ignored unless `--check-comments` is active. Even when comments are checked,
each input must contain at least one noncomment data record.

Comparable records are paired in order. Each record is split into
whitespace-separated fields. The files must have the same number of records
and the paired records must have the same number of fields.

The supported field forms are:

| Form | Comparison |
| --- | --- |
| `1.25`, `-3e-8` | Finite real-number comparison. |
| `(1.25,-3e-8)` | Componentwise finite complex-number comparison. |
| `label=1.25` | The label must match exactly; the payload is compared numerically. |
| `label=(1,2)` | The label must match exactly; both components are compared numerically. |
| Other text | Exact text comparison after whitespace tokenization. |

NaN, infinity, overflow, nonzero values that silently underflow to zero, and
recognized malformed numeric or complex forms are input errors. They never
compare equal, even if the same invalid spelling occurs in both files. Tokens
that do not match a recognized numeric or invalid-numeric form are treated as
text and compared exactly.

### Numerical rule

For expected value `e` and actual value `a`, define

```text
difference = abs(a - e)
scale      = max(abs(a), abs(e))
limit      = atol + rtol * scale
```

The values match when `difference <= limit`. The rule is symmetric: exchanging
the expected and actual operands does not change the result. Complex values use
the same rule independently for their real and imaginary components.

The absolute term handles zero and near-zero values. If the expected value is
exactly zero and `rtol < 1`, the largest accepted actual magnitude is

```text
atol / (1 - rtol)
```

This is slightly larger than `atol` because the relative term uses the larger
operand. For example, the generic defaults accept approximately
`1.0000100001e-12` against zero.

There is no separate rule saying that any two individually small values are
equal. Their difference must still satisfy the formula. Thus `8e-13` and
`-8e-13` fail with the generic defaults because their difference is `1.6e-12`.
References for quantities that mathematically must vanish should use a literal
zero when possible.

### Sign-insensitive comparison

With `--ignore-signs`, real fields are compared using their absolute values.
For a complex field, either the original value or one whole-complex sign flip
may match:

```text
(real, imaginary) ~= (actual_real, actual_imaginary)
(real, imaginary) ~= (-actual_real, -actual_imaginary)
```

Independent component sign changes and arbitrary phase changes are not
accepted. Labels and textual fields still match exactly.

### Diagnostics and status

`mycomp.pl` reports at most the first ten differing fields, while retaining the
overall difference count.

| Status | Meaning |
| ---: | --- |
| `0` | The files match. |
| `1` | Valid inputs contain one or more differences. |
| `2` | Usage, file, or numeric-input error. |

## Directory comparison

### Invocation

```sh
test/compare.pl [OPTIONS] REF_DIR
```

| Option | Meaning |
| --- | --- |
| `--actual DIR` | Read actual results from `DIR`. Default: current directory. |
| `--exclude FILE` | Exclude one exact top-level filename. May be repeated. This is not a glob. |
| `--ignore-signs` | Forward sign-insensitive behavior to ordinary text comparisons. |
| `--strict` | Enable physical-output classification, strict comments, class-specific tolerances, and detection of unreferenced physical outputs. |

Reference directories are flat. Their entries must be regular files; nested
directories are rejected. Hidden files and filenames containing spaces or
shell metacharacters are supported.

Each ordinary reference file must have a same-named regular file in the actual
directory. Ordinary text references are passed to `mycomp.pl`. Golden `.bin`
and `.h5` references are prohibited because such outputs must use semantic
validation instead.

Without `--strict`, additional actual-side files are ignored. With `--strict`,
an additional actual-side file fails only if `PhysicalOutput.pm` classifies its
name as a physical output. Inputs, logs, and other unclassified files remain
outside that census.

### Built-in tolerance policies

`compare.pl` selects the following text-comparison policy before applying an
optional per-reference override:

| Mode | Fields | `atol` | `rtol` |
| --- | --- | ---: | ---: |
| Generic or non-strict | All | `1e-12` | `1e-5` |
| Strict, unclassified file | All | `1e-12` | `1e-5` |
| Strict, nonspectral physical output | All | `1e-9` | `1e-5` |
| Strict spectral text output | First field | `1e-12` | `1e-5` |
| Strict spectral text output | Remaining fields | `1e-12` | `2e-2` |

For spectral tables, the first field is normally the frequency or temperature
grid. It remains tightly checked while dependent values may vary by roughly two
percent. "First" means the first whitespace-separated field, not the first
field recognized as numeric.

Strict mode is structurally stricter but not numerically tighter in every
case. Compared with generic mode, it uses a larger absolute tolerance for
nonspectral physical outputs and a larger relative tolerance for spectral
dependent values.

## Per-reference tolerance files

An ordinary text reference may have a companion named `FILENAME.tol` in the
same reference directory:

```text
ref/FILENAME
ref/FILENAME.tol
```

The companion is reference metadata. It is not compared with an actual-side
file, and no `FILENAME.tol` should be generated in the actual directory.
Direct `mycomp.pl` calls do not discover companions; discovery belongs to
`compare.pl`.

### Format

A single number overrides only the absolute tolerance:

```text
# Permit platform-dependent residue around zero.
1e-9
```

Two numbers override the absolute and relative tolerances, in that order:

```text
# Absolute tolerance, then relative tolerance.
1e-11 5e-4
```

The two values may also occur on separate lines:

```text
# Absolute
1e-11

# Relative
5e-4
```

Blank lines and full-line comments whose first non-whitespace character is `#`
are ignored. Inline comments are not supported. After removing blank and
comment lines, the file must contain exactly one or two whitespace-separated
numeric tokens. The values obey the same finite, representable, nonnegative
requirements as command-line tolerances.

### Precedence

Tolerance policy is resolved in this order:

1. Start with the generic `mycomp.pl` defaults.
2. Apply strict physical or spectral defaults when applicable.
3. Replace `atol` with the first companion value.
4. Replace `rtol` only when the companion supplies a second value.

For a strict spectral text output, the companion applies to fields after the
first field. The first-field grid policy remains fixed at `atol=1e-12` and
`rtol=1e-5` and cannot be loosened by `FILENAME.tol`.

For generic and strict nonspectral physical comparisons, the effective
companion tolerances apply to every field. Sign-insensitive comparison and
strict comment checking remain active independently of tolerance overrides.

### Association and errors

The lowercase `.tol` suffix is reserved for comparison metadata in reference
directories. Every active companion must name an existing ordinary text
reference after removing one `.tol` suffix.

The following are configuration errors with status `2`:

- an empty or malformed companion;
- more than two values;
- invalid, negative, nonfinite, overflowing, or underflowing values;
- an orphan companion without a matching ordinary reference;
- a companion targeting an output declared in `.physical-outputs`;
- a companion targeting a `.bin` or `.h5` reference;
- recursive metadata names such as `FILENAME.tol.tol`.

An actual-side `.tol` file is never consulted and cannot change comparison
policy.

`--exclude FILENAME` excludes both the reference and its companion. Explicitly
excluding only `FILENAME.tol` disables the override while leaving `FILENAME` to
be compared with the normal policy. Excluded companions are not parsed.

## Physical-output classification

`PhysicalOutput.pm` recognizes these exact basenames:

```text
absolute_energies.dat
annotated.dat
custom
customfdm
energies.nrg
raw-dm.h5
raw.h5
report.nrg
states.nrg
subspaces.dat
td
tdfdm
```

It also recognizes case-sensitive spectral basenames with one of these
prefixes and a `.bin` or `.dat` suffix:

```text
chit_ corr_ gt_ i1t_ i2t_ orbspin_ spec_ specq_ spin_
```

The classifier controls four behaviors:

1. Strict text-comparison tolerances and comment checking.
2. Detection of unreferenced physical outputs under `--strict`.
3. Which filenames may be declared in `.physical-outputs`.
4. Which stale files `clean-physical-output.pl` removes.

When a new generated physical-output filename is introduced, update the
classifier and its contract tests rather than duplicating filename knowledge in
individual runners.

## Semantic output validation

Some outputs are unsuitable for checked-in golden bytes. A reference directory
can declare these in `.physical-outputs`. Each nonblank, noncomment line has the
form:

```text
VALIDATOR FILENAME
```

Example:

```text
# Native real spectral records
binary-real spec_FDM_dens_d-d.bin

# HDF5 density-matrix output
hdf5 raw.h5
```

Manifest filenames must be classified physical-output basenames. A declaration
cannot duplicate an ordinary reference or another declaration. Excluded
filenames are not validated. Manifest processing is active with or without
`--strict`.

Semantic validators check that generated data is structurally usable; they do
not compare its numerical values against a reference dataset.

| Validator | Compatible output | Validation performed |
| --- | --- | --- |
| `binary-real` | Spectral `.bin` | Nonempty native records of two `double` values; complete records and finite values. |
| `binary-complex` | Spectral `.bin` | Nonempty native records of three `double` values; complete records and finite values. |
| `hdf5` | `.h5` | `h5dump` can decode the file, at least one dataset exists, and the dump contains no NaN or infinity tokens. |
| `subspaces` | `subspaces.dat` | Iteration grammar, increasing iterations, row counts, invariant arity and uniqueness, and `kept <= total`. |
| `states` | `states.nrg` | Iteration and subspace grammar, finite energies and vectors, dimensions, vector counts, coefficient kind, and orthonormality. |

The state validator requires the printed `norm-1` magnitude to be at most
`1e-8`. Recomputed squared norms may differ from one by at most `1e-4`, and
pairwise inner-product magnitudes may be at most `1e-4`.

When an ordinary `energies.nrg` reference is active, the state validator also
requires the actual state dump to cover the same actual energy-dump iterations,
subspaces, eigenvalue counts, ordering, and canonicalized energies. The
ordinary `energies.nrg` file is independently compared with its reference using
the selected text tolerances. Excluding `energies.nrg` disables this additional
state-to-energy correlation.

The HDF5 validator requires an executable `h5dump` in `PATH`. Missing validator
tooling is a configuration/environment error, not a numerical mismatch.

## Stale-output cleanup

```sh
test/clean-physical-output.pl
```

The cleaner accepts no arguments and operates in the current directory. It
removes `DONE` and regular entries classified by `PhysicalOutput.pm`. It does
not remove `data`, `param`, logs, reference metadata, or arbitrary unclassified
files. It refuses to recursively remove an output-like directory.

Running cleanup before the executable prevents stale results or a stale `DONE`
marker from making a failed calculation appear successful.

## Wrapper scripts

`mycomp-ignoresigns.pl` and `compare-ignoresigns.pl` are thin wrappers. They add
`--ignore-signs`, forward all remaining arguments, and replace themselves with
the underlying comparator. They do not implement a separate parser or
tolerance policy.

## Adding a regression reference

For an ordinary deterministic text output:

1. Add the expected file under the case's `ref/` directory.
2. Use the standard tolerance policy unless the output has a demonstrated,
   domain-specific numerical requirement.
3. If needed, add `FILENAME.tol` with the narrowest suitable absolute and
   optional relative tolerance.
4. Prefer a literal zero in the reference for quantities that mathematically
   vanish.
5. Run the individual regression and the comparator contract test.

For binary, HDF5, or nondeterministic structural output, add an appropriate
`.physical-outputs` declaration instead of committing a golden binary file.
Add a new validator only when the existing validators cannot express the
required invariant.

When adding a new physical-output basename or spectral naming family, update
`PhysicalOutput.pm`, cleanup expectations, strict extra-output tests, and any
manifest compatibility rules together.

## Testing the comparison tools

The primary contract test is registered with CTest as `comparator_contract`.
From the repository root, run:

```sh
ctest --test-dir build -R '^comparator_contract$' --output-on-failure
```

It can also be run directly without a configured build tree:

```sh
perl test/scripts/comparators.t test
```

The contract covers text parsing, tolerance boundaries, nonfinite and
unrepresentable inputs, sign handling, companion tolerance files, strict-mode
policies, special filenames, exclusions, semantic manifests, state-vector
consistency, and stale-output cleanup.

Useful syntax checks are:

```sh
perl -c test/mycomp.pl
perl -c test/compare.pl
perl -c test/scripts/comparators.t
```

## Suite-local comparison helpers

`tools/integ/check_scalar.pl` is a specialized integration-test helper. It
requires explicit absolute and relative tolerances and exactly one finite
numeric scalar in each input. It uses the same combined tolerance formula but
does not provide the record grammar, first-field policy, sidecar discovery, or
directory semantics of the central comparators.

Tests that invoke system `cmp` or `diff` intentionally require exact bytes or
text and are not governed by the tolerance policies in this document.
