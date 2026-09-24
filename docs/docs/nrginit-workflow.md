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

## Wilson-Chain Reconstruction

The default remains `tri=old`, the existing high-precision recursion in
`wilson.m`. The runtime method default remains `tridiag_method=lanczos`.
For scalar normal-state chains, select the unsquared
Rutishauser/Gragg-Harrod (RKPW) backend in the `[param]` block:

```ini
tri=rkpw
```

Then run `nrginit` normally. This generates the full requested coefficient
table in Mathematica using machine arithmetic for reconstruction only.
There is no C++ reconstruction of that table in this mode.

For C++ runtime reconstruction, select both settings:

```ini
tri=cpp
tridiag_method=rkpw
```

Here Mathematica uses machine RKPW to generate the coefficients needed by
the initial cluster, through `Ninit`; the C++ runtime reconstructs its chain
from the exported star data. With `tri=cpp` and the default
`tridiag_method=lanczos`, the initializer keeps its old high-precision seed
recursion. `tri=none` also honors `tridiag_method=rkpw` for the `Ninit` seed,
without changing its existing coefficient-table output policy.
`tridiag_method` does not override an explicit `tri=old`, `orth`, or `rkpw`.
Unknown `tri` or `tridiag_method` strings are initialization errors, even when
the runtime-method setting would otherwise be unused.

### Precision and support

RKPW does **not** eliminate upstream arbitrary-precision work. `prec` still
controls discretization energies, hybridization integrals, and normalized
initial amplitudes. Its default is 30 for `tri=rkpw`, as for `tri=cpp` and
`tri=none`; the `tri=old` default stays 1000. Increase upstream precision when
the band integration requires it. Increasing `prec` does not turn RKPW
reconstruction into an arbitrary-precision algorithm.

The scalar reconstruction inserts signed energies in shell order,
`+de[m], -deminus[m]`, from the outside inward. It discards exactly zero
amplitudes and merges exactly identical **machine** energies, preserving
their first occurrence and combining amplitudes with a scaled Euclidean
norm. Before merging or reconstruction, all amplitudes receive one common
power-of-two scale that puts the largest amplitude in `[1,2)`. This preserves
common-scale invariance even for the least subnormal input and prevents
amplitude-norm overflow. A nonzero amplitude that rounds to zero under this
normalization is an error, not discarded support. The binary preprocessing
does not change the machine-arithmetic recurrence or the exported `du/dv[0]`.
The backend neither groups all positive energies first nor removes small
coefficients using a bandwidth-relative threshold.

For `count=DISCNMAX+1`, the effective unique nonzero support must have at
least `count` poles. When support equals `count`, the output has `count`
onsite entries, `count-1` positive hoppings, and a terminal zero hopping.
Requests beyond support fail instead of padding a longer chain. A short
requested prefix still incorporates every pole; only coefficient storage and
the insertion sweeps are capped. For `tri=cpp`/`none`, the initializer uses
`DISCNMAX=Ninit` for this check.

Nonfinite or nonrepresentable machine inputs, a nonzero input rounded to
zero, and numerical breakdown are errors. Scaled norms avoid unnecessarily
squaring tiny amplitudes or tail hoppings, but representability still limits
how deep a chain can go. RKPW does not silently retry in arbitrary precision.
Independent scalar channels are supported; matrix, rung, and superconducting
chains are not. The existing `sc` and `sc2` algorithms are unchanged.

For RKPW only, `bandrescale` must be a finite positive machine real. The
initializer also checks the final rescaled coefficient tables, including
onsite adjustments from `gap`, `shift0`, and bulk fields. Nonfinite results
or any nonzero coefficient rounded to zero are errors; representable
subnormals and the exact terminal zero hopping are allowed. A very small
positive bandwidth can therefore fail even when the unscaled reconstruction
succeeded. These checks also apply to RKPW seeds with `tri=cpp`/`none`, but
do not change legacy-backend scaling behavior.

The normalized high-precision `du/dv[0]` amplitudes remain available for star
output. With `disccheck` present, RKPW reports the initial-state normalization,
first-moment error, and square-root-of-variance error for each channel. It
does not reconstruct higher Lanczos vectors for orthogonality diagnostics.

### Ownership and testing

`tri` is an initializer setting. `tridiag_method` is also read by the
initializer to select the `cpp`/`none` seed and by the runtime to select its
reconstruction method. Treat these settings, along with `prec`, `Ninit`,
`nrxi`, and discretization inputs, as generation-locked: changing them
requires regenerating `data` before running the solver. See the
[parameter reference](parameter-reference.md) for the ownership inventory.

The focused Mathematica-only test needs no model generation or C++ build:

```sh
sh test/nrginit/test_rkpw "$PWD" "$(command -v math)"
```

It checks an independent high-precision asymmetric-star reference, spectral
nodes and weights, the analytic flat-band chain through a long tail, support
and representability boundaries, least-subnormal amplitude scale invariance,
final physical-coefficient scaling, and initializer dispatch/seed behavior.
When the Mathematica test suite is configured (Mathematica detected and
`SYM_ALL` enabled), its CTest name is `nrginit_rkpw`.
The `nrginit_rkpw_pipeline` test additionally generates fresh `data` with both
RKPW execution paths, runs the solver against existing physical references,
and checks a nontrivial initial cluster. Build `nrg` before running it.

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
- [Floquet model construction](floquet-nrginit.md)
- [Output format reference](output-formats.md)
- [Runtime flow](runtime-flow.md)
