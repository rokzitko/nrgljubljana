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

The main script is `nrginit/nrginit.m`.

Its top-level flow is compact and helpful to keep in mind:

1. define package search paths
2. load `sneg.m`
3. load `initialparse.m`
4. parse `param`
5. load `initial.m`
6. call `makedata["data"]`

In other words, `nrginit.m` is primarily the entrypoint and dispatcher; most of the actual initialization logic lives in the auxiliary Mathematica files it loads.

## Hand-Off To The C++ Runtime

The key artifact produced by `nrginit` is `data`.

From the C++ point of view, `data` contains:

- header information such as symmetry and channel count
- seed eigenspectra
- seed operator blocks
- chain coefficient tables

That file is then read by `InputData<S>` in `c++/read-input.hpp`.

## Relationship To `param`

Both sides of the project look at `param`, but for different purposes:

- Mathematica initialization uses it to define the model and generate `data`
- the C++ runtime uses it to decide how to run the iterative solver and postprocessing phases

This is why changes to parameter semantics often need to be checked in both `nrginit/` and `c++/params.hpp`.

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
- [Runtime flow](runtime-flow.md)
