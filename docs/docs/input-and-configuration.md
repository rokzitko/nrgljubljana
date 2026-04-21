# Input and Configuration

This page describes the startup side of the codebase: how `param` and `data` are interpreted, how symmetry is selected, and how the initial runtime objects are constructed.

## Runtime Inputs

The native executable expects two files in the working directory:

- `param`: user-facing runtime configuration
- `data`: generated seed data for the first NRG step

The split is intentional:

- `param` controls how the run should proceed
- `data` contains the precomputed initial Hamiltonian basis, operator blocks, and coefficient tables

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
6. initialize resume-related bookkeeping

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

For the next stage after startup, see [Iteration engine](iteration-engine.md).
