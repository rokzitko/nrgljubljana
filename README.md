[![Build Status](https://travis-ci.org/rokzitko/nrgljubljana.svg?branch=master)](https://travis-ci.org/rokzitko/nrgljubljana)

"NRG Ljubljana" is a flexible framework for performing large-scale
numerical renormalization group (NRG) calculations for quantum
impurity problems. It is highly extensible without sacrificing
numerical efficiency.

*Copyright (C) 2006-2019 Rok Zitko*

The framework "NRG Ljubljana" is a set of interrelated computer codes
for performing numerical renormalization group (NRG) calculations for
quantum impurity problems, described by models such as the Kondo
exchange (s-d) model or the Anderson single impurity model, and their
multi-impurity and multi-channel generalizations. It also contains a
number of tools for analyzing the results (thermodynamic properties,
such as magnetic and charge susceptibility, entropy and heat capacity;
expectation values of arbitrary operators; spectral functions). It is
user-friendly, in the sense that it is easy to set up new types of
problems (Hamiltonians, perturbation terms, etc.) and the output is
formatted and annotated for easy interpretation, parsing and plotting.
It efficiently handles problems with different symmetries, such as
spin SU(2) symmetry, charge SU(2) symmetry, Z_2 reflection symmetry
(parity), etc.

To achieve a high degree of flexibility without sacrificing numerical
efficiency, "NRG Ljubljana" is composed of a hierarchy of modules:
high level modules are written in a mixture of functional and
procedural Mathematica code, while the low level numerically intensive
parts are programmed in the object oriented approach in the C++
language. The foundation of the framework is a Mathematica package for
performing calculations with non-commutative second quantisation
operators, SNEG. The next layer is a Mathematica program which defines
the Hamiltonian, the basis of states, and the physical operators of
interest: with the help of SNEG, Hamiltonian and operators can be
defined using the familiar second-quantization expressions. This
program performs the diagonalization of the initial Hamiltonian and
prepares the input for the NRG iteration proper.

For efficiency, NRG iteration is performed by a separate C++ program:
for a typical problem, most of the time (90%) is spent in the LAPACK
dsyev routine which solves the eigenvalue problem. There is very
little housekeeping overhead due to the tasks required by the NRG
iteration; "NRG Ljubljana" is thus suitable for performing large scale
NRG calculations on computer clusters.

1. Features

   - all parameters, model definitions and observables configurable at run-time
   - support for a large number of different symmetry types
   - flexible truncation schemes (energy cut-off truncation, avoidance
     of trunction within gaps, etc.)
   - density-matrix NRG (DM-NRG), complete Fock space (CFS) and full-density matrix
     (FD-NRG) for all symmetry types
   - FDM calculation of expectation values at finite temperatures
   - various spectral-function broadening schemes & stand-alone tools
   - self-energy trick calculations
   - arbitrary number of channels
   - calculations with real and complex numbers
   - support for superconducting bands, spin-polarized bands, etc.
   - support for global operators (i.e., operators defined on the
     Wilson chain sites)
   - calculation of temperature-dependent conductance, G(T)
   - automatic exact diagonalisation of the initial Hamiltonian with
     automagic generation of the basis
   - flexibility in the choice of operators whose thermodynamical averages
     are computed, they can be expressed using operators of second
     quantization
   - dynamic spin susceptibility, dynamic charge susceptibility, etc.
     It is possible to compute arbitrary spectral functions for any
     pair of local operators.
   - multiple logarithmic discretization schemes (Wilson/Krishnamurthy,
     Yoshida/Whitaker/Oliveira, Campo/Oliveira, ODE scheme)
   - support for non-flat bands (i.e. cosine band that arrises from
     tight-binding description of leads in quantum transport problems,
     or arbitrary hybridization as needed in DMFT)
   - object-oriented code for easy maintenance and expandability
   - very high numeric efficiency with optimized-for loops in the most
     numerically demanding parts of the code (chiefly the recalculation
     of irreducible matrix elements of operators)
   - configurable verbosity level
   - automatic recording of time elapsed in various parts of the program
   - monitoring of memory usage
   - formated output files for easy interpretation of the results
   - internal consistency checks, assertions, parameter compatibility
     and reasonableness checks that reduce the possibility of undetected
     bugs


2. License

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

   The full text of the GPL General Public License can be found
   in file LICENSE.

3. Compiling and dependencies

   NRG Ljubljana is very portable and it should work without any modification
   on any modern Linux distribution and, with some tweaking, on any Unix or
   Unix-like operating system with a good standards-compliant C++ compiler. It
   has been reported to me that it can also be compiled under Windows.

   The following libraries are required to compile the C++ part of the
   NRG code and the tools:

    * LAPACK and BLAS linear algebra libraries
    * Boost C++ libraries
    * GNU Scientific library (GSL)
    * GNU MP Bignum Library (GMP) for arbitrary-precision numerics

   Due to the heavy use of template metaprogramming in Boost libraries, a
   high-quality standards-compliant C++ compiler must be used. Tested to
   work with GCC, Clang and Intel C++ compiler. As of 2019, the code
   is written in C++14.

   Wolfram Reasearch Mathematica must be installed for running the
   Mathematica part of the NRG code. Versions 5 through 12 have been
   tested. Mathematica is only required for the initialization of the
   problem (basis construction, diagonalisation of the initial
   Hamiltonian, transformations of the operator matrices, etc.) which is
   relatively fast. When "NRG Ljubljana" is used on a cluster, it is
   therefore sufficient to have Mathematica installed on a single

   computer (for example on the cluster host computer), while the
   numerically demanding (C++) part of the program can be ran on the
   cluster nodes.

   Since Nov 2019, the code uses cmake for the configuration stage.
   The compilation thus consists of the following steps:

   mkdir build
   cd build
   cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/nrgljubljana/
   make
   make install

   For debugging, add -DCMAKE_BUILD_TYPE=Debug to cmake.


4. Contributing to "NRG Ljubljana"

   If you make improvements to "NRG Ljubljana", you are encouraged to
   share them with other users. Bug reports (and fixes) are very welcome
   as well.  The contact information is in the next section.


5. Contact information:

   "NRG Ljubljana" home-page: http://nrgljubljana.ijs.si/

   Rok Zitko
   "Jozef Stefan" Institute
   F1 - Theoretical physics
   Jamova 39
   SI-1000 Ljubljana
   Slovenia

   rok.zitko@ijs.si


6. Acknowledgements

   The development of the "NRG Ljubljana" framework started during
   author's PhD studies at the Faculty for mathematics and physics of the
   University of Ljubljana, and the "Jozef Stefan" Institute, Ljubljana,
   Slovenia. Discussions and collaboration with prof. Janez Bonca, prof.
   Anton Ramsak, dr. Jernej Mravlje and dr. Tomaz Rejec from the F1,
   Theoretical Physics department are acknowledged. I'm also grateful to
   prof. Thomas Pruschke, Robert Peters and Oliver Bodensiek from the
   University in Goettingen for many very fruitful discussions. I thank
   Marcus Greger from the University in Augsburg for contributing
   optimized routines for the spectral function calculation. Nils
   Wentzell helped me make the switch to cmake build system.
