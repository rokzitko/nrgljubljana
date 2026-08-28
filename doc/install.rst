.. highlight:: bash

.. _install:

Compiling NRG Ljubljana from source
===================================


Installation steps
------------------

#. Clone the ``rokzitko/NRGLjubljana`` repository from GitHub. This checks out
   the current development branch; select a release tag as described below
   when a released version is required::

     $ git clone https://github.com/rokzitko/NRGLjubljana nrgljubljana.src

#. Create and move to a new directory where you will compile the code::

     $ mkdir nrgljubljana.build && cd nrgljubljana.build

#. In the build directory call cmake, including any additional custom CMake options, see below::

     $ cmake -DCMAKE_INSTALL_PREFIX=path_to_install_dir ../nrgljubljana.src

#. Compile the code, run the tests and install the application::

     $ make
     $ make test
     $ make install

Versions
--------

Release tags use dotted numeric names. Go into the source directory and list
them in descending version order::

     $ cd nrgljubljana.src
     $ git tag --list '[0-9]*' --sort=-version:refname

Select one of the listed tags in detached mode. For example::

     $ git switch --detach 2026.06.1

Then follow steps 2 to 4 above to compile the code. The current set of
published versions is also available on the `GitHub releases page
<https://github.com/rokzitko/nrgljubljana/releases>`_ and the `tags page
<https://github.com/rokzitko/nrgljubljana/tags>`_.

Custom CMake options
--------------------

The compilation of ``NRG Ljubljana`` can be configured using CMake options::

    cmake ../nrgljubljana.src -DOPTION1=value1 -DOPTION2=value2 ...

Common options include:

* Installation path: ``-DCMAKE_INSTALL_PREFIX=path``
* Debug build: ``-DCMAKE_BUILD_TYPE=Debug``
* Disable testing: ``-DBuild_Tests=OFF``
* Build the documentation: ``-DBuild_Documentation=ON``
* Enable application-level OpenMP regions:
  ``-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON``
* Select the MKL threading layer for ``mkl_rt`` builds:
  ``-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU``

BLAS/LAPACK integer ABI
-----------------------

The default ``NRGLJUBLJANA_BLAS_ILP64=OFF`` configuration uses the
32-bit-integer LP64 interface. ``NRGLJUBLJANA_BLAS_ILP64=ON`` selects the
64-bit-integer ILP64 interface. This setting is independent of pointer size.
The selected BLAS and LAPACK libraries and NRG Ljubljana must use the same
integer ABI and symbol convention; mixing LP64 and ILP64 may link but is
unsafe.

For MKL, use ``-DBLA_VENDOR=Intel10_64lp`` with LP64 or
``-DBLA_VENDOR=Intel10_64ilp`` with ILP64. OpenBLAS is an ILP64 option only
when the installed library was built for 64-bit integers and exports compatible
symbols. Apple Accelerate and the standard Conda variants use LP64. Prefer the
direct MKL ILP64 interface over ``mkl_rt``; a dispatcher build also requires a
matching runtime interface layer, which NRG Ljubljana does not select.

Parallelism
-----------

The default parallelism model is BLAS/LAPACK-internal threading. Use a threaded MKL or OpenBLAS build and control numerical kernel threads with variables such as ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS`` and ``OMP_NUM_THREADS``. Application-level OpenMP regions are disabled by default to avoid linking multiple OpenMP runtimes into the same executable.

For MKL builds using the ``mkl_rt`` dispatcher, ``-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU`` (or ``INTEL`` or ``LLVM``) selects the expected MKL backend. This links the compiler OpenMP runtime through CMake's ``OpenMP::OpenMP_CXX`` target without enabling NRG Ljubljana's own OpenMP regions.

The option ``-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON`` is intended only for expert runs that intentionally use simultaneous diagonalisation scheduling through ``diag_mode=OpenMP`` and ``diagth``. If BLAS/LAPACK is also threaded, this is nested parallelism. CMake checks the visible link line for mixed GNU/Intel/LLVM OpenMP runtime families, and the executable reports the detected MKL/OpenBLAS/OpenMP/MPI threading configuration at startup.
