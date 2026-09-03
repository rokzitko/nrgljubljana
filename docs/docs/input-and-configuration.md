# Input and Configuration

This page describes the startup side of the codebase: how `param` and `data` are interpreted, how symmetry is selected, and how the initial runtime objects are constructed. For the complete user-facing key list, see the [parameter reference](parameter-reference.md).

## Runtime Inputs

The native executable expects two files in the working directory:

- `param`: user-facing initialization and runtime configuration
- `data`: generated seed data for the first NRG step

The split is intentional:

- `param` controls both how `nrginit` generates the model and how `nrg` runs
- `data` contains the precomputed initial Hamiltonian basis, operator blocks, and coefficient tables

`data` is a generated same-release hand-off rather than a user-edited result
format. Regenerate it after changing model, symmetry, discretization, operator,
or Wilson-chain settings.

## End-User Syntax And Ownership

A portable file uses named blocks and one case-sensitive `key=value`
assignment per line:

```text
[param]
symtype=QS
Lambda=2.0

[extra]
# Model-specific values
```

Both parsers ignore blank lines and whole-line `#` comments when indentation
uses ordinary ASCII spaces. Do not use tabs because the `nrginit` parser does
not strip them. Inline comments, quoted strings, and duplicate keys should not
be used. Initializer-only keys such as `model`, `band`, and `Nmax` remain in the
file when `nrg` starts and are therefore printed as unused by the C++ parser.
This is expected, but the same report also exposes misspelled runtime keys.

The [parameter reference](parameter-reference.md#portable-file-syntax) defines
the common syntax, boolean values, stage ownership, defaults, and constraints.

## `param` Parsing

`Params` is defined in `c++/params.hpp`.

Important characteristics:

- user-facing parameters are stored as many `param<T>` members
- `Params` also holds derived runtime quantities such as `Nmax`, `Nlen`, and channel/combinatorics metadata
- the constructor parses a named block from `param`, normally `[param]`
- if an `[extra]` block is present, it is parsed into `extra_params`

Startup flow in the constructor is roughly:

1. parse the requested block from `param`
2. apply recognized values to the registered `param<T>` members
3. print any leftover keys as unused settings
4. optionally parse `[extra]`
5. validate the resulting configuration

The C++ parser rejects duplicate keys. Unrecognized keys are reported rather
than rejected so that the shared file can contain initializer settings.

## `data` Parsing

The `data` file is read by `InputData<S>` in `c++/read-input.hpp`.

The `InputData` constructor performs these steps:

1. read the file header and version information
2. extract symmetry name, channel count, `Nmax`, and the number of initial subspaces
3. instantiate the correct `Symmetry<S>` backend through `set_symmetry(...)`
4. load the seed eigenspectrum into `DiagInfo<S>`
5. load the seed `f` operators into `Opch<S>`
6. read the remaining operator blocks and coefficient tables
7. finalize chain-length-related derived parameters in `Params`

After `InputData` has finalized `Nmax`, `NRG_calculation` scans the exact
checkpoint directory when `resume=true`. Deferring discovery until this point
ensures the scan covers the actual runtime iteration range. Every archive in
the discovered prefix is then validated before existing result files are
replaced.

The marker in `data` versions the producer/consumer hand-off, but the layout is
not a general interchange contract. Its current role is described in the
[output format reference](output-formats.md#nrginit-hand-off-files).

## Symmetry Selection

Symmetry objects are created through `set_symmetry(...)` in `c++/mk_sym.hpp`.

That function:

1. resets the global `Invar` metadata
2. stores the selected `symtype` into `Params`
3. derives channel and combination counts
4. constructs the matching concrete symmetry backend
5. loads backend-specific quantum number tables

This means the symmetry object is a structural dependency for almost everything that follows:

- interpreting invariant subspaces
- generating new subspaces
- building Hamiltonians
- recalculating operators

## Coefficients And Derived Lengths

`Coef<S>` is populated while reading the `data` blocks. After that, `determine_Nmax_Nlen(...)` in `read-input.hpp` finalizes the effective number of iteration steps, taking `substeps` into account.

This is an important boundary: before `data` is loaded, `Params` is not fully initialized for iteration; after `data` is loaded, the runtime has enough information to start the NRG loop.

## Ownership After Initialization

Once `NRG_calculation` is constructed in `c++/nrg-general.hpp`, the startup objects are held like this:

```text
NRG_calculation
  ├── Params P
  ├── InputData<S> input
  │    ├── shared_ptr<Symmetry<S>> Sym
  │    ├── DiagInfo<S> diag
  │    ├── Operators<S> operators
  │    └── Coef<S> coef
  ├── Stats<S> stats
  ├── ThermoStore<S> store
  └── BackiterStore store_all
```

From this point on, the runtime has everything it needs to start `run_phase(RUNTYPE::NRG, ...)`.

## Why This Matters For Contributors

If you change any of the following, this page's area of the code is where you should start:

- adding a new parameter
- changing parameter validation
- altering the `data` file format
- selecting or configuring symmetry backends
- modifying chain-length or coefficient initialization rules

Any user-visible change in this area also requires an update to the [parameter
reference](parameter-reference.md).

For the next stage after startup, see [Iteration engine](iteration-engine.md).
