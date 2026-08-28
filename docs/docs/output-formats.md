# Output Format Reference

This page describes files emitted by the current release of `nrginit` and
`nrg`. It is a description of current behavior, not a promise that every layout
will remain byte-for-byte compatible across releases. In particular, native
binary files, HDF5 trees, diagnostic dumps, Mathematica expressions, console
text, and temporary workdir archives are implementation-dependent.

For text tables, treat whitespace as the delimiter and header labels as
authoritative. Column padding, trailing spaces, and the textual representation
of floating-point values are not field boundaries.

## Locations And Lifecycle

Both programs use the process current directory for inputs and persistent
results:

- `nrginit` reads `param` and normally writes `data` and `mmalog`.
- `nrg` reads `param` and `data` and writes result files beside them.
- `nrg -w DIR` and `NRG_WORKDIR=DIR` select only the parent of a uniquely named
  temporary workdir. They do not relocate persistent results.

Most result streams truncate a file when that output is opened. Files from
features that are no longer requested are not removed, so a reused calculation
directory can contain stale spectra, diagnostics, HDF5 files, or completion
markers. Final result files are generally written in place rather than through
an atomic temporary-file rename; an interrupted run can leave partial output.

Use the program's exit status as the primary success signal. With `done=true`,
`nrg` creates an empty `DONE` file after the requested phases return. It does
not remove an existing `DONE` at startup, so a marker in a reused directory
does not prove that the latest invocation succeeded.

## Artifact Summary

| Artifact | Created when | Current role |
| --- | --- | --- |
| `data` | Normal `nrginit` run | Generated, same-release input for `nrg`. |
| `mmalog` | `nrginit` loads the initializer | Volatile Mathematica diagnostic log. |
| `data.in` | Template generation is selected | Input to the `instantiate` tool, not directly to `nrg`. |
| `td` | Normal NRG measurements produce rows | Shell thermodynamics text table. |
| `custom` | Main NRG phase starts | Shell expectation-value text table. |
| `tdfdm` | An FDM workflow computes `rhoFDM` | Fixed-temperature FDM thermodynamics. |
| `customfdm` | `fdmexpv=true` | FDM expectation values. |
| `*_dens_*.dat` | A requested spectrum or response is saved | Real-frequency, Matsubara, or temperature-dependent text table. |
| `*_dens_*.bin` | `savebins=true` for a real-frequency spectrum | Unversioned native binary bins. |
| `annotated.dat` | `dumpannotated>0` | RG-flow diagnostic. |
| `energies.nrg` | `dumpenergies=true` | Eigenspectrum diagnostic. |
| `states.nrg` | `dumpstates=true` | Energy/vector diagnostic. |
| `report.nrg` | `reportdiagonal>0` | Per-sector diagonal-observable report. |
| `absolute_energies.dat` | `dumpabsenergies=true` | Absolute-energy diagnostic. |
| `subspaces.dat` | `dumpsubspaces=true` | Per-subspace dimension diagnostic. |
| `raw.h5` | `h5raw=true` during the NRG phase | Implementation-dependent HDF5 state. |
| `raw-dm.h5` | `h5raw=true` during a second phase | Implementation-dependent HDF5 DM/FDM state. |
| `DONE` | `done=true` and requested phases return | Empty completion marker with stale-file caveats. |
| `unitaryN`, `rhoN`, `rhofdmN` | Required by selected multiphase workflows | Temporary Boost archives inside the generated workdir. |

C++ and Mathematica extension modules can create additional model-specific
files that are outside this inventory.

## Thermodynamics: `td`

`td` is a whitespace-separated text table. It has one `#` header row followed
by one row per reported shell:

```text
#               T           <Sz^2>              <Q>            <Q^2>              <E>            <E^2>                C                F                S
       1.020139447     0.2751638928   -0.03678911966     0.8696561051     0.4349139256     0.3007523743     0.1116022517     -2.278648456      2.713562382
```

The base column order is:

```text
T, symmetry-specific fields, <E>, <E^2>, C, F, S
```

The fields have these conventions:

| Field | Meaning |
| --- | --- |
| `T` | Effective shell temperature `energyscale/betabar`, in bandwidth units. |
| symmetry fields | Thermal moments of the conserved quantum numbers, named by the header. |
| `<E>` | Dimensionless `beta * <H>`. |
| `<E^2>` | Dimensionless `beta^2 * <H^2>`. |
| `C` | Heat capacity in units of `k_B`. |
| `F` | `F/(k_B*T)`, equal to `-ln(Z)`. |
| `S` | Entropy in units of `k_B`. |

Symmetry fields are inserted immediately after `T` in this order:

| Symmetry | Fields |
| --- | --- |
| `NONE`, `P`, `PP` | none |
| `U1` | `<Q>`, `<Q^2>` |
| `SU2` | `<Q^2>` |
| `SL` | `<Q>`, `<Q^2>`, `<sQ^2>` |
| `SL3` | `<Q1>`, `<Q1^2>`, `<sQ1^2>`, `<Q2>`, `<Q2^2>`, `<sQ2^2>`, `<Q3>`, `<Q3^2>`, `<sQ3^2>` |
| `QS`, `QSLR`, `QSC3` | `<Sz^2>`, `<Q>`, `<Q^2>` |
| `QSZ`, `QSZLR` | `<Sz^2>`, `<Sz>`, `<Q>`, `<Q^2>` |
| `QST`, `QSTZ` | `<Sz^2>`, `<Tz^2>`, `<Q>`, `<Q^2>` |
| `QSZTZ` | `<Sz>`, `<Sz^2>`, `<Tz>`, `<Tz^2>`, `<Q>`, `<Q^2>` |
| `QJ` | `<Jz^2>`, `<Q>`, `<Q^2>` |
| `SPSU2`, `SPSU2LR`, `SPSU2C3` | `<Sz^2>` |
| `SPSU2T` | `<Sz^2>`, `<Tz^2>` |
| `SPU1`, `SPU1LR` | `<Sz^2>`, `<Sz>` |
| `ISO`, `ISO2`, `ISOLR`, `ISO2LR` | `<Sz^2>`, `<Q^2>` |
| `ISOSZ`, `ISOSZLR` | `<Sz^2>`, `<Sz>`, `<Q^2>` |
| `DBLSU2` | `<Q1^2>`, `<Q2^2>` |
| `DBLISOSZ` | `<Sz^2>`, `<Sz>`, `<Q1^2>`, `<Q2^2>` |
| `DBLQSZ` | `<Sz^2>`, `<Sz>`, `<Q1>`, `<Q1^2>`, `<Q2>`, `<Q2^2>` |

With the default `calc0=true`, the first row describes the seed problem and is
followed by the iterative-shell rows. Do not infer a shell index from the line
number alone.

`width_td` is a minimum display width and `prec_td` is the number of
general-format significant digits. Values can exceed the requested width.

## Expectation Values: `custom`

`custom` has two `#` header rows. The first numbers the columns; the second
names them:

```text
#               1                2
#               T              n_d
       1.020139447     0.9630254109
```

The first column is the same shell temperature convention as `td`. Remaining
columns are expectation values for local and global singlet operators encoded
in `data`. Local operators are ordered lexically by name, followed by global
operators in lexical order; this need not match the order written in `ops`.

For complex calculations, the file stores only the real component of each
expectation value. `width_custom` and `prec_custom` control minimum width and
significant digits.

`custom` is opened for the main phase even when no scalar operators are
available. In that case it contains the header and temperature column only.

## FDM Tables

### `tdfdm`

`tdfdm` has one header and normally one fixed-temperature row:

```text
#               T            E_fdm            C_fdm            F_fdm            S_fdm
              0.1     -2.041644463      1.843302768     -2.316274192       2.74629729
```

`T`, `E_fdm`, and `F_fdm` are in bandwidth units. `C_fdm` and `S_fdm` are in
units of `k_B`. Formatting uses `width_td` and `prec_td`.

### `customfdm`

`customfdm` uses the same two-header grammar and operator ordering as `custom`.
It is created by `fdmexpv=true` and receives a row when the run reaches
`fdmexpvn`. If that iteration is outside the executed range, the file can
contain headers only.

## Spectral And Response Files

### Filename Grammar

The current basename grammar is:

```text
<prefix>_<algorithm>_dens_<operator1>-<operator2>[-u|-d].<suffix>
```

The optional `-u` and `-d` select polarized components. Operator names are
inserted without filename escaping, so use filesystem-safe operator names in
automation.

| Prefix | Quantity selected by the parameter file |
| --- | --- |
| `spec` | Doublet spectral function. |
| `corr` | Singlet correlator. |
| `spin` | Triplet response. |
| `orbspin` | Orbital-triplet response. |
| `specq` | Quadruplet spectral function. |
| `gt` | Temperature-dependent conductance. |
| `i1t` | First temperature-dependent transport moment. |
| `i2t` | Second temperature-dependent transport moment. |
| `chit` | Temperature-weighted susceptibility `k_B*T*chi(T)`. |

Current algorithm identifiers include `FT`, `FTmats`, `DMNRG`, `DMNRGmats`,
`CFS`, `CFSgt`, `CFSls`, `FDM`, `FDMgt`, `FDMls`, `FDMmats`, `GT`, `I1T`,
`I2T`, and `CHIT`.

### Broadened Real-Frequency `.dat`

When `broaden=true`, real-frequency spectra are headerless text:

```text
omega real_density
```

With `reim=true`, a third column is appended:

```text
omega real_density imag_density
```

Frequencies are in bandwidth units. Negative mesh points are written first in
increasing numerical order, followed by positive mesh points in increasing
order. The logarithmic mesh has no zero-frequency row. `prec_xy` controls
significant digits. These columns are broadened densities or response values;
the unbroadened integrated bin weights are stored only in `.bin` output.

### Matsubara `.dat`

Matsubara files are headerless and always contain three columns:

```text
omega_n real imag
```

Frequencies are ascending. Fermionic frequencies are `(2*n+1)*pi*T`; bosonic
frequencies are `2*n*pi*T`. The operator type selects the applicable sequence.

### Temperature-Dependent `.dat`

`gt`, `i1t`, `i2t`, and `chit` outputs are headerless curves sorted by
increasing temperature:

```text
T value
```

With `reim=true`, an imaginary-value column is appended.

For the `chit` prefix, `value` is `k_B*T*chi(T)` (`chi/beta`), not `chi(T)`.
Recovering `chi(T)` therefore requires division by `k_B*T` in the units used by
the calculation.

### Unbroadened `.bin`

With `savebins=true`, each record consists of native host `double` values:

```text
omega, real_weight
```

With `reim=true`, each record contains three doubles. Positive bins are written
first from low to high frequency, followed by negative bins from near zero
toward more-negative frequency. The file has no magic bytes, version, record
count, scalar declaration, or byte-order declaration. It is an
implementation-dependent same-platform artifact, not a portable interchange
format.

## Diagnostic Text Files

Diagnostic layouts are useful for inspection but are not stable interchange
formats.

### `annotated.dat`

Each iteration is a block terminated by a blank line. There is no header or
explicit iteration number. With the default `dumpgroups=true`, a row has this
shape:

```text
energy (quantum numbers) [total multiplicity]
```

With `dumpgroups=false`, it is:

```text
energy quantum-number-1 quantum-number-2 ...
```

Entries are sorted by energy. At least `dumpannotated` entries are retained,
and a degeneracy group crossing the cutoff is completed. `dumpscaled`,
`dumpabs`, `dumpEscale`, `dumpprecision`, and `grouptol` control energy units
and presentation.

### `energies.nrg` And `states.nrg`

Both use blocks beginning with:

```text
===== Iteration number: N
```

`energies.nrg` then contains `Subspace:` sections and raw eigenvalues. Optional
`corr=` and `crit=` records are enabled by `dumpcorr` and `dumpcrit` when
`dumpenergies=true`. `states.nrg` adds the vector information still available
at the point where the diagnostic is written. These files use diagnostic
formatting and should not be treated as complete restart data.

Iteration labels are zero-based. With `calc0=true` and `Ninit=0`, a seed block
and the first iterative block can both carry label `0`.

### `report.nrg`

The report is organized by iteration and invariant sector:

```text
=== Report N=3
Sector I=0 1
I=0 1 n=0 E=0 <n_d>=0.946016949
```

State index `n` is zero-based. Local operators are emitted first in lexical
order, followed by global operators in lexical order. A missing operator block
can make rows differ in field count. Complex diagonal values are represented by
their real component.

### `absolute_energies.dat`

This blank-line-separated diagnostic contains iteration and `Subspace:`
headers followed by energies relative to the final absolute ground state in
bandwidth units. Its displayed iteration label is one-based, unlike the labels
in `energies.nrg`.

### `subspaces.dat`

The current block grammar is:

```text
Iteration N
len_dm=<count>
I=<quantum numbers> kept=<count> total=<Hamiltonian dimension>
```

`N` is zero-based. `total` is the full subspace dimension and can exceed the
number of eigenpairs computed by a partial eigensolver.

## HDF5 Output

`h5raw=true` opens `raw.h5` for the main phase and `raw-dm.h5` for a requested
second phase in overwrite mode. The files have no root schema-version marker.
Their hierarchy depends on `h5all`, `h5last`, `h5ham`, `h5ops`, `h5vectors`,
`h5U`, and `h5struct`.

Current content can include:

```text
params/<keyword>
<iteration>/eigen/<invariant>/...
<iteration>/hamiltonian/<invariant>/matrix
<iteration>/structure/<invariant>/ancestors
<iteration>/U/<invariant>/...
<iteration>/<operator-type>/<name>/<I1>/<I2>/matrix
stats/...
expv/<name>
store/...
fdm/...
```

Parameter and scalar values are commonly stored as one-element datasets rather
than HDF5 scalar datasets. A complex matrix is split into a real dataset at the
base path and an imaginary dataset with the suffix `-imag`. Iteration indexing
is not uniform for every seed and runtime group. With `h5all=true`, seed
operators and eigenspectrum vectors are stored independently of `h5ops` and
`h5vectors`; those parameters control later iteration saves. Consumers must
inspect the actual file produced by their release and parameter set.

## `nrginit` Hand-Off Files

### `data`

`data` is whitespace-separated generated input. The current producer writes a
`#!9` marker followed by symmetry, channel, chain, seed-eigenspectrum,
operator, and coefficient blocks. Complex files encode values textually as
`(real,imag)` and identify themselves in the comment header.

The format is coupled to the runtime reader and contains inferred dimensions
and symmetry-specific blocks. Treat it as a same-release hand-off: do not edit
it, parse it as a result format, or assume that generated files are portable
between incompatible releases. Keep the original `param` and regenerate
`data` when generation-locked settings change.

### `data.in` And Symbolic Artifacts

Template options emit `data.in`, which can contain `DIAG` and external matrix
placeholders. It must be processed by `instantiate` before use as `data`.
Optional `basis`, `ham_*`, `op.*`, and `opf*` files contain unversioned
Mathematica syntax whose ordering and expression form can vary with the kernel
version.

### `mmalog`

`mmalog` records initializer progress, paths, Mathematica messages, generated
expressions, and timings after the log is opened. Early launcher and parser
messages can remain console-only. The content is intended for troubleshooting,
not machine parsing.

## Temporary Workdir Files

Multiphase calculations can store `unitaryN`, `rhoN`, and `rhofdmN` in the
generated workdir. They are Boost binary archives used only by the running
solver. They are not portable across Boost versions, architectures, scalar
types, or releases.

Individual archives are written through temporary files and renamed after
serialization. During normal shutdown, the workdir is removed recursively;
`removefiles=false` can delay per-file deletion but does not preserve the whole
workdir after a normal run. A crash can leave it behind.

## Standard Output, Standard Error, And `log`

`nrg` writes build information, effective parameters, unused settings,
iteration progress, diagnostics, spectral summaries, total energy, timings,
and memory information primarily to standard output. Warnings and fatal
preambles are not consistently separated onto standard error. ANSI styling is
used when output is a terminal and disabled when output is redirected.

Console text is diagnostic and can change. Even familiar lines such as `Total
energy:` should not be treated as a general structured-output interface.

`nrg` does not create a file named `log`. Repository test wrappers do so by
combining standard output and standard error through `tee`. Users who want the
same behavior can run:

```sh
nrg 2>&1 | tee log
```

With a pipeline, enable shell `pipefail` or inspect the `nrg` status explicitly
so that a successful `tee` process does not hide solver failure.

See the [parameter reference](parameter-reference.md) for every switch that
enables or formats these artifacts.
