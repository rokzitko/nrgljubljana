.. highlight:: bash

.. _install:

Compiling NRG Ljubljana from source
===================================


Installation steps
------------------

#. Download the source code of the latest stable version by cloning the ``rokzitko/NRGLjubljana`` repository from GitHub::

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

To use a particular version, go into the directory with the sources, and look at all available versions::

     $ cd nrgljubljana.src && git tag

Checkout the version of the code that you want::

     $ git checkout 2019/12

and follow steps 2 to 4 above to compile the code.

Custom CMake options
--------------------

The compilation of ``NRG Ljubljana`` can be configured using CMake-options::

    cmake ../nrgljubljana.src -DOPTION1=value1 -DOPTION2=value2 ...

+-----------------------------------------------------------------+-----------------------------------------------+
| Options                                                         | Syntax                                        |
+=================================================================+===============================================+
| Specify an installation path                                    | -DCMAKE_INSTALL_PREFIX=path                   |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build in Debugging Mode                                         | -DCMAKE_BUILD_TYPE=Debug                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Disable testing (not recommended)                               | -DBuild_Tests=OFF                             |
+-----------------------------------------------------------------+-----------------------------------------------+
| Build the documentation                                         | -DBuild_Documentation=ON                      |
+-----------------------------------------------------------------+-----------------------------------------------+
| Enable application-level OpenMP regions                         | -DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON          |
+-----------------------------------------------------------------+-----------------------------------------------+
| Select MKL threading layer for ``mkl_rt`` builds                 | -DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU       |
+-----------------------------------------------------------------+-----------------------------------------------+

Parallelism
-----------

The default parallelism model is BLAS/LAPACK-internal threading. Use a threaded MKL or OpenBLAS build and control numerical kernel threads with variables such as ``MKL_NUM_THREADS``, ``OPENBLAS_NUM_THREADS`` and ``OMP_NUM_THREADS``. Application-level OpenMP regions are disabled by default to avoid linking multiple OpenMP runtimes into the same executable.

For MKL builds using the ``mkl_rt`` dispatcher, ``-DNRGLJUBLJANA_MKL_THREADING_LAYER=GNU`` (or ``INTEL`` or ``LLVM``) selects the expected MKL backend. This links the compiler OpenMP runtime through CMake's ``OpenMP::OpenMP_CXX`` target without enabling NRG Ljubljana's own OpenMP regions.

The option ``-DNRGLJUBLJANA_ENABLE_APP_OPENMP=ON`` is intended only for expert runs that intentionally use simultaneous diagonalisation scheduling through ``diag_mode=OpenMP`` and ``diagth``. If BLAS/LAPACK is also threaded, this is nested parallelism. CMake checks the visible link line for mixed GNU/Intel/LLVM OpenMP runtime families, and the executable reports the detected MKL/OpenBLAS/OpenMP/MPI threading configuration at startup.
