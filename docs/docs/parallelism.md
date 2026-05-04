# Parallelism

NRG Ljubljana's default parallelism philosophy is simple: BLAS/LAPACK owns numerical threading.

The expensive dense diagonalisation kernels should run inside a threaded BLAS/LAPACK implementation such as MKL or OpenBLAS. Application-level OpenMP regions are disabled by default so the executable does not accidentally link a second OpenMP runtime into the same process. Mixing runtimes such as GNU `libgomp`, Intel `libiomp5`, and LLVM `libomp` is unsafe and can crash as soon as threaded numerical kernels are used.

## Default Model

- Build with `NRGLJUBLJANA_ENABLE_APP_OPENMP=OFF`, which is the default.
- Use threaded MKL or threaded OpenBLAS for BLAS/LAPACK.
- Control numerical threads with environment variables such as `MKL_NUM_THREADS`, `OPENBLAS_NUM_THREADS`, and `OMP_NUM_THREADS`.
- Use MPI for parallelising independent diagonalisation work across ranks when desired.
- Size MPI jobs as `mpi_ranks * blas_threads <= allocated_cpus` unless oversubscription is intentional and controlled by the scheduler.

In this mode the code has serial application-level scheduling, while each LAPACK call may use many threads internally.

## Expert Application OpenMP

`NRGLJUBLJANA_ENABLE_APP_OPENMP=ON` enables OpenMP regions in NRG Ljubljana itself. This affects simultaneous diagonalisation scheduling through `diag_mode=OpenMP` and `diagth`, plus a few non-BLAS loops.

Use this only intentionally. If `diagth>1` and BLAS/LAPACK also uses more than one thread, the program is nested-parallel: several diagonalisation tasks can run at once, and each task can also create BLAS/LAPACK worker threads. This can oversubscribe CPUs and can expose incompatible OpenMP runtimes.

The CMake configuration checks the visible BLAS/LAPACK and application OpenMP link line and fails if more than one OpenMP runtime family is detected. It also rejects explicitly sequential MKL selections. Opaque dispatcher libraries such as `mkl_rt` may not reveal the final runtime at configure time, so the executable prints runtime diagnostics at startup.

## Startup Diagnostics

The `nrg` executable reports parallel configuration at startup on rank 0, including:

- selected BLAS/LAPACK vendor from the build
- relevant thread-control environment variables
- application OpenMP build status
- loaded OpenMP runtime libraries from `/proc/self/maps` when available
- MKL version, threading layer, max threads, BLAS-domain max threads, and dynamic mode when MKL service symbols are visible
- OpenBLAS config, core, thread count, and threading model when OpenBLAS reporting symbols are visible
- MPI rank count times BLAS/LAPACK thread count, with an oversubscription warning when this exceeds online CPUs

If startup reports multiple OpenMP runtime families, rebuild or adjust the BLAS/LAPACK selection before trusting threaded runs.
