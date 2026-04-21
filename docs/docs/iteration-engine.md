# Iteration Engine

This page breaks down the main forward iteration in terms of the actual coordination functions in `c++/core.hpp`.

## The Main Loop

The outer loop is `nrg_loop(...)` in `c++/core.hpp`.

At a high level it does:

```text
for each Step:
  iterate(...)
```

`Step` packages the current iteration index, energy scale, runtime mode, and convenience helpers used throughout the loop.

## One Iteration: `iterate(...)`

`iterate(...)` is the main coordinator for a single NRG step.

Its structure is:

1. construct `SubspaceStructure` from the previous `DiagInfo`
2. build a `TaskList` from the subspace structure
3. optionally persist structure information to HDF5
4. run diagonalization via `do_diag(...)`
5. run post-diagonalization processing via `after_diag(...)`
6. trim operator matrices and drop dense eigenvectors no longer needed
7. print timing/memory summary

## Diagonalization: `do_diag(...)`

`do_diag(...)` is responsible for the part of the loop that may restart if too few states were computed.

The flow is:

1. print iteration information through `step.infostring()`
2. report Wilson-chain coefficients via `Sym->show_coefficients(...)`
3. attempt to obtain the eigenspectra with `load_or_compute_diag(...)`
4. initialize the energy reference with `initialize_diag_energy_reference(...)`
5. prepare truncation metadata with `prepare_diag_for_truncation(...)`
6. catch `NotEnough` and optionally retry with a larger `diagratio`

This is where the eigensolver policy and restart behavior meet.

## Hamiltonian Construction

The actual block Hamiltonian for one invariant subspace is built by `hamiltonian(...)`.

That function:

1. computes the ancestor list for the target `Invar`
2. constructs `SubspaceDimensions`
3. initializes the diagonal from the previous shell's corrected eigenvalues
4. delegates symmetry-specific matrix filling to `Sym->make_matrix(...)`

This is the point where the generic iteration engine hands off to the symmetry backend.

## Post-Diagonalization: `after_diag(...)`

`after_diag(...)` is the second half of an iteration. It is where the freshly computed eigenspectra are turned into data products and next-step state.

The main branches are:

1. if this is an NRG step:
   - handle Floquet-specific postprocessing or split eigenvectors into blocks
   - finalize metadata such as scale and total-energy shifts
   - persist outputs for the current iteration
2. optionally recalculate operators and measure before truncation
3. perform actual truncation with `diag.truncate_perform()`
4. archive iteration state into `ThermoStore` and `BackiterStore`
5. if not at the last step, recalculate irreducible operators for the next step
6. optionally recalculate operators and measure after truncation
7. optionally check operator sum rules and write operator HDF5 output

## Measurement And Recalculation Boundary

The helper `recalculate_and_measure(...)` in `core.hpp` is a useful conceptual boundary. It does exactly two things:

1. `Oprecalc::recalculate_operators(...)`
2. `calculate_spectral_and_expv(...)`

That means operator recalculation and measurement are treated as a paired operation in the main loop.

## `docalc0(...)`

Before the first real NRG step, `run_phase(...)` may call `docalc0(...)` if `P.calc0` is enabled.

This phase performs measurements using the seed data from `data` without advancing the iteration loop. It is effectively a pre-iteration measurement stage.

## Function Map

The most useful function-level reading order in `core.hpp` is:

1. `nrg_loop(...)`
2. `iterate(...)`
3. `do_diag(...)`
4. `hamiltonian(...)`
5. `after_diag(...)`
6. `recalculate_and_measure(...)`
7. `archive_iteration_state(...)`

## Related Pages

- [Runtime flow](runtime-flow.md)
- [Data structures](data-structures.md)
- [State and persistence](state-and-persistence.md)
- [Symmetry framework](symmetry-framework.md)
