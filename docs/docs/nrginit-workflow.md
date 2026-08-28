# nrginit Workflow

The `nrginit/` part of the repository is the Mathematica-side initialization layer that prepares the seed input for the C++ runtime.

## Why It Exists

The C++ executable does not build the model Hamiltonian from symbolic expressions. Instead, that work is done before the iterative run starts.

The initialization layer is responsible for:

- parsing model parameters
- constructing the initial basis
- defining operators
- diagonalizing the initial Hamiltonian
- generating the `data` file consumed by the C++ runtime

## Entry Script

The installed entry point is the `nrginit` shell launcher. Run it from a
calculation directory containing `param`:

```sh
nrginit
```

It runs the command `math` by default. To select another kernel executable,
pass it as the only positional argument:

```sh
nrginit /path/to/WolframKernel
```

`nrginit -h` prints the short launcher usage. The launcher does not accept
model settings as command-line options; those belong in `param`.

The launcher locates and loads `nrginit/nrginit.m` from the installation.

Its top-level flow is compact and helpful to keep in mind:

1. define package search paths
2. load `sneg.m`
3. load `initialparse.m`
4. parse `param`
5. load `initial.m`
6. call `makedata["data"]`

In other words, `nrginit.m` is primarily the entrypoint and dispatcher; most of the actual initialization logic lives in the auxiliary Mathematica files it loads.

## Working Directory And Results

`nrginit` preserves the caller's current directory. Its normal persistent
artifacts are:

| Artifact | Purpose |
| --- | --- |
| `data` | Generated seed input consumed by `nrg`. |
| `mmalog` | Mathematica messages, progress, generated expressions, and timings. |
| `data.in` | Template form generated instead of `data` by template options. |

Advanced `WRITE` and template options can also produce `basis`, `ham_*`,
`op.*`, and `opf*` Mathematica expressions. The `writedir` parameter affects
those expressions only; it does not relocate `data`, `data.in`, or `mmalog`.

A normal successful invocation exits with status zero and prints `Success!`.
Initialization errors normally print `Aborting.` and exit nonzero. Always check
the process status before running `nrg`: a failed rerun does not remove an old
`data` file. `mmalog` is diagnostic output and can vary with the Mathematica
version.

## Hand-Off To The C++ Runtime

The key artifact produced by `nrginit` is `data`.

From the C++ point of view, `data` contains:

- header information such as symmetry and channel count
- seed eigenspectra
- seed operator blocks
- chain coefficient tables

That file is then read by `InputData<S>` in `c++/read-input.hpp`. Although the
file begins with a format marker, it is a generated same-release hand-off, not
a supported user-editable result format. Regenerate it after changing model,
symmetry, discretization, operator, or Wilson-chain settings.

## Relationship To `param`

Both sides of the project look at `param`, but for different purposes:

- Mathematica initialization uses it to define the model and generate `data`
- the C++ runtime uses it to decide how to run the iterative solver and postprocessing phases

This is why changes to parameter semantics often need to be checked in both `nrginit/` and `c++/params.hpp`.

The [parameter reference](parameter-reference.md) identifies initializer-only,
shared, and runtime settings. Initializer-only keys appearing as `Unused
settings` in the later `nrg` startup report are expected.

## Run A Minimal Model

[Getting started](getting-started.md#first-end-to-end-calculation) contains a
complete SIAM `param` file and commands that run `nrginit` followed by `nrg`.
Use that example to verify a new installation before adding custom models,
band files, or `[extra]` parameters.

## Repository Role

From a contributor perspective, `nrginit/` is the bridge between symbolic/high-level model specification and the low-level iterative runtime.

If you are changing:

- model definitions
- initialization-stage operator content
- `data` file layout
- Mathematica-generated metadata needed by the runtime

then this directory is part of the implementation, not just a frontend wrapper.

## Recommended Companion Reading

- `nrginit/nrginit.m`
- `c++/read-input.hpp`
- `c++/params.hpp`
- [Input and configuration](input-and-configuration.md)
- [Parameter reference](parameter-reference.md)
- [Output format reference](output-formats.md)
- [Runtime flow](runtime-flow.md)
