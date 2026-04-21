# Tools

The `tools/` directory is the repository's standalone utility layer. These executables are separate from the main `nrg` runtime, but they are part of the same workflow ecosystem: some prepare discretization-related inputs, some transform intermediate files, and many postprocess spectral or tabulated output.

This page documents the tools as a subsystem rather than as a flat command list.

## Design Pattern

Most tools follow a simple structure:

1. a tiny `main()` in `tools/<name>/<name>.cc`
2. the actual implementation in `tools/<name>/<name>.hpp` or `.cc`
3. shared parsing or I/O helpers from `tools/common/` and `c++/`
4. regression coverage in `test/tools/<name>/`

This keeps each tool easy to build and test in isolation while still reusing the project's numerical and file-format utilities.

Representative examples:

- `tools/adapt/adapt.cc` sets up the sign and parameter filename, then delegates to `NRG::Adapt::Adapt`
- `tools/broaden/broaden.cc` is a thin wrapper around `NRG::Broaden::Broaden`
- `tools/h5write/h5write.cc` and `tools/unitary/unitary.cc` follow the same pattern

The main exception is `matrix`, which is parser-driven and includes generated scanner/parser sources.

## Build Integration

All tool executables are declared centrally in `tools/CMakeLists.txt`.

The build logic is intentionally uniform:

- each tool is an `add_executable(...)`
- all tools link against the same warning target and core dependency stack
- the tools include the project C++ headers directly
- the tools use the build-tree runtime path for the core library

This means tool code is not a separate mini-project. It is a thin executable layer built on the same numerics, matrix, HDF5, and helper infrastructure as the main solver.

## Shared Tool Infrastructure

The most important reusable helpers live in `tools/common/`:

- `tools/common/parser.hpp`
  - block lookup in parameter-like files
  - simple key/value parsing
  - `get_or_default(...)` helpers for typed access
- `tools/common/io.hpp`
  - small file-opening wrappers
  - next-data-line handling for comment-heavy text files
- `tools/common/tabulated.hpp`
  - parsing and writing of tabulated `(x, y)` data
  - interval reporting
  - sign-restricted loaders used by several spectral tools

Many tools also reuse common code from `c++/`, especially:

- `misc.hpp`
- `basicio.hpp`
- `io.hpp`
- `h5.hpp`
- `traits.hpp`

The practical consequence is that changes to low-level parsing or I/O helpers can affect both the main runtime and multiple standalone tools.

## Tool Families

## Input And Chain Preparation

### `adapt`

Purpose:

- solve the discretization ODE and generate an adapted mesh for NRG input preparation

Typical inputs:

- `param`
- discretization-related tabulated input files

Implementation notes:

- `adapt` is a front-end tool with its own parameter parser and calculation object
- it leans on `tools/common/parser.hpp` and tabulated-data helpers
- it sits close to the start of the workflow, before the main NRG iteration

### `nrgchain`

Purpose:

- generate or inspect Wilson-chain-related data derived from the discretization setup

Implementation notes:

- `nrgchain` has its own lightweight parser logic in `tools/nrgchain/parser.h`
- it is conceptually paired with `adapt`: one computes the discretization mesh, the other translates that information into chain data used downstream

These tools are the closest thing in `tools/` to the front-end preparation phase of the project.

## Spectral Postprocessing

### Tools In This Family

- `broaden`
- `bw`
- `hilb`
- `kk`
- `resample`
- `specmoments`

### Shared Role

These tools operate on tabulated spectral or response data after the main `nrg` run. They are used to:

- broaden discrete spectral output into smoother curves
- perform Kramers-Kronig or Hilbert-transform-style conversions
- resample onto a new grid
- compute moments or summaries of spectral weight

### Shared Structural Traits

- heavy use of tabulated text-file parsing
- strong reliance on `tools/common/tabulated.hpp`
- a lot of command-line option parsing and validation
- frequent interaction with GSL-based interpolation or integration routines

### Notable Distinctions

- `broaden` and `bw` are both broadening-oriented, but they represent different broadening strategies and input assumptions
- `hilb` and `kk` are transform tools rather than smoothing tools
- `resample` is focused on interpolation and grid conversion
- `specmoments` is summary-oriented rather than transformation-oriented

If you are changing common spectral-file parsing or interpolation behavior, this is the family most likely to be affected broadly.

## Matrix And Linear-Algebra Utilities

### Tools In This Family

- `diag`
- `matrix`
- `unitary`

### Shared Role

These tools operate closer to matrix-level structure than to full NRG workflows. They are useful for:

- inspecting matrices
- diagonalizing text-based matrix input
- composing or transforming matrices
- debugging or validating intermediate algebraic steps

### `diag`

`diag` is the simplest representative: it consumes matrix input and exposes basic diagonalization-related functionality.

### `unitary`

`unitary` applies chained unitary-style transformations to matrix inputs. It is used both as a real utility and as a testable linear-algebra helper.

### `matrix`

`matrix` is structurally different from the other tools because it includes generated parser sources. That makes it the most parser-like tool in the tree and the one most sensitive to lexer/parser regeneration or parser-adjacent code changes.

## Averaging And Integration Utilities

### Tools In This Family

- `binavg`
- `intavg`
- `integ`
- `mats`
- `tdavg`

### Shared Role

These tools work on already-produced numerical data and compute derived tables or aggregates:

- averages over bins or multiple runs
- integration over tabulated values
- Matsubara-oriented summaries
- temperature-dependent averaging helpers

### Common Characteristics

- line-oriented parsing of simple text or binary data
- small command-line interfaces
- straightforward failure modes that are heavily covered by regression tests

These are often the easiest tools to read when you want to understand the project's utility conventions because their data flow is linear and self-contained.

## File And Format Utilities

### `h5write`

Purpose:

- write scalar or matrix content into an HDF5 file at a named dataset path

Why it matters architecturally:

- it is the clearest standalone bridge between external text input and the repository's HDF5 helpers
- it reuses the same HDF5 support used by the main runtime

### `matrix`

Although grouped above as a matrix utility, `matrix` also belongs here conceptually because it acts as a parser-driven format conversion and inspection tool.

## Relationship To The Main Runtime

The tools are not wrappers around `nrg`. Most are independent executables.

Still, they are tightly coupled to the rest of the project in three ways:

1. they share the same low-level helper code and numerical libraries
2. they operate on files generated by the main workflow
3. they are validated through repository-local regression tests that encode expected behavior and error messages

This means tool behavior is part of the project's compatibility surface.

## Testing

Tool tests live in `test/tools/`.

The usual pattern is:

1. copy minimal inputs into a test working directory
2. run the tool via a small `run` script
3. check produced output or expected failure text

That test style is especially useful for catching:

- argument parsing regressions
- error-message changes
- file-format handling mistakes
- differences in external-library behavior

## Suggested Reading Order For Tool Work

If you are about to change tool code, a good reading order is:

1. `tools/CMakeLists.txt`
2. the specific tool's `main()` file
3. the tool implementation header/source
4. the shared helpers in `tools/common/`
5. the corresponding tests under `test/tools/`

If the tool touches HDF5, matrix parsing, or numerical helpers, also inspect the reused code under `c++/`.

## Extending The Tool Layer

When adding a new tool, keep it aligned with the existing pattern:

- one directory per tool under `tools/`
- thin `main()`
- reusable parsing and I/O helpers where possible
- regression coverage under `test/tools/`
- minimal dependency footprint beyond what the repository already builds against

This keeps the tool layer coherent and prevents it from drifting into a collection of unrelated one-off programs.
