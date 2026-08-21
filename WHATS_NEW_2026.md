# NRG Ljubljana 2026.09: Major New Release

This release summarizes development since `2024.12`, with major improvements in GPU accelerator support, performance, memory efficiency, physical capabilities, and workflows.

## CUDA And Numerical Performance

- NVIDIA CUDA support has been added for numerical linear algebra operations.
- CUDA-enabled builds can select `diag=cuda` to use cuSOLVER diagonalisation routines (real and complex).
- Selected operator-recalculation matrix products can be offloaded with `mult=cuda`.
- Eigenvector blocks and intermediate accumulators remain cached on the GPU during recalculation.
- BLAS/LAPACK integration has been improved, including support for ILP64 numerical libraries for really large problems.
- Threaded BLAS, MPI scheduling, and oversubscription diagnostics are handled more consistently.
- CPU diagonalization now defaults to the divide-and-conquer LAPACK routines `dsyevd` and `zheevd`.

## Performance And Memory

- All unnecessary matrix, eigenspace, and density-matrix copies have been eliminated.
- Seed operators, diagonalization data, and obsolete eigenvector representations are released earlier.
- FDM and DMNRG spectral kernels now reuse weights, contractions, and existing operator-block information.
- FDM back-iteration exploits the diagonal structure of discarded-state density matrices.
- Memory requirements for diagonalization workspaces are now reported more clearly.
- Completed Wilson shells no longer redundantly retain full eigenvector and operator-block matrices.
- Thermodynamic histories and density-matrix back-iteration data now use separate, compact storage.

## New Physics Capabilities

- Multiple local phonon modes are supported, including independent bosonic cutoffs.
- Tensor-product construction of phonon bases and mode-resolved operators has been generalized.
- Superconducting Wilson chains can be instantiated using Nambu onsite and hopping coefficients.
- Symmetry triangle inequalities are enforced before constructing reduced operator matrix elements.
- Experimental Floquet-NRG support introduces Floquet bases and quasienergy-aware truncation.
- Three-channel `QST` calculations received important low-energy Hamiltonian and recalculation fixes.
- Orbital-triplet operator generation and several SNEG symbolic-algebra operations were corrected.

## Spectra And Observables

- `report.nrg` can list low-lying states together with diagonal observables, aiding fixed-point and sub-gap-state identification.
- `broaden` accepts arbitrary user-provided output-frequency meshes.
- `adapt --flat Gamma` directly supports constant hybridization functions.
- Thermal Fermi and Bose kernels are stable at extreme energies and near the Bose pole.
- Large Matsubara meshes no longer suffer from a small-integer index limitation.
- FDM partition-function accumulators consistently retain high numerical precision.
- Level-flow energies can be reported in user-selected or physical energy units.
- Raw HDF5 output and broadening sum-rule diagnostics received correctness fixes.

## Workflows And Reliability

- Parameter files, spectral meshes, and truncated inputs now receive substantially stricter validation.
- The new `instantiate`/`nrgspawn` workflow can run prepared model templates without invoking Mathematica for every parameter point.
- Basis, Hamiltonian, and operator blocks can be saved and reused in parameter sweeps.
- Checkpoints and density matrices are written atomically and loaded transactionally.
- The `DONE` marker is created only after all requested calculation stages complete successfully.
- Conda packaging and build coverage now include broader Linux, macOS, ARM, BLAS, and compiler configurations.

## Installation

NRG Ljubljana is available from conda-forge:

    conda install -c conda-forge nrgljubljana

Packages are available for Linux x86-64, Linux aarch64, macOS x86-64,
and macOS arm64. Source builds remain appropriate when CUDA support or
nonstandard numerical-library configurations are required.

## Requirements And Compatibility

- The C++ runtime requires a C++20 compiler.
- CPU diagonalization now defaults to the divide-and-conquer LAPACK routines dsyevd and zheevd Numerical results should remain equivalent within floating-point tolerances.
- BLAS/LAPACK and OpenMP runtime selection is handled more explicitly. Users combining MPI with threaded numerical libraries should review their rank and thread settings to avoid oversubscription.
- CUDA support is optional and must be enabled explicitly in source builds. It is not currently included in the standard conda-forge packages.
