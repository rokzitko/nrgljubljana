/*
 "NRG Ljubljana" - Numerical renormalization group for multiple
 impurities and an arbitrary number of channels

 Copyright (C) 2005-2020 Rok Zitko

   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the
   Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

   Contact information:
   Rok Zitko
   F1 - Theoretical physics
   "Jozef Stefan" Institute
   Jamova 39
   SI-1000 Ljubljana
   Slovenia

   rok.zitko@ijs.si
*/

#ifndef _nrg_general_h_
#define _nrg_general_h_

#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <complex>
#include <numeric>
#include <limits>
#include <memory>
#include <string>
using namespace std::string_literals;
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <deque>
#include <set>
#include <stdexcept>

// C headers
#include <cassert>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstring>
#include <unistd.h>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/optional.hpp>

// ublas matrix & vector containers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
using namespace boost::numeric;
using namespace boost::numeric::ublas; // keep this!

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>
namespace atlas = boost::numeric::bindings::atlas;

// Serialization support (used for storing to files and for MPI)
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>

// MPI support
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
using namespace fmt::literals;

#include <range/v3/all.hpp>

// Support for compiler dependant optimizations

#ifdef __GNUC__
#define PUREFNC __attribute__((pure))
#define CONSTFNC __attribute__((const))
#else
#define PUREFNC
#define CONSTFNC
#endif

#include "nrg-general.h"
#include "nrg-lib.h" // exposed interfaces for wrapping into a library
#include "portabil.h"
#include "debug.h"
#include "misc.h"
#include "openmp.h"
#include "mp.h"
#include "traits.h"
#include "workdir.h"
#include "params.h"
#include "numerics.h"
#include "io.h"
#include "time_mem.h"
#include "outfield.h"

// This is included in the library only. Should not be used if a cblas library is available.
#ifdef CBLAS_WORKAROUND
 #define ADD_
 #include "cblas_globals.c"
 #include "cblas_dgemm.c"
 #include "cblas_zgemm.c"
 #include "cblas_xerbla.c"
#endif

inline const size_t MAX_NDX = 1000; // max index number
inline const double WEIGHT_TOL = 1e-8; // where to switch to l'Hospital rule form

#include "invar.h"
#include "eigen.h"
#include "operators.h"
#include "subspaces.h"
#include "store.h"
#include "step.h"
#include "stats.h"
#include "spectral.h"
#include "coef.h"
#include "tridiag.h"
#include "diag.h"
#include "symmetry.h"
#include "matrix.h"
#include "recalc.h"
#include "read-input.h"
#include "spectrum.h"
#include "algo.h"
#include "dmnrg.h"
#include "splitting.h"
#include "output.h"
#include "oprecalc.h"
#include "measurements.h"
#include "truncation.h"
#include "core.h"
#include "mk_sym.h"

template <typename S> class NRG_calculation {
private:
  MPI mpi;
  Params P;
  Stats<S> stats;
  MemTime mt; // memory and timing statistics
public:
  auto run_nrg(Step &step, IterInfo<S> &iterinfo, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diag0,
               AllSteps<S> &dm, std::shared_ptr<Symmetry<S>> Sym) {
    diag0.states_report(Sym->multfnc());
    auto oprecalc = Oprecalc<S>(step.runtype, iterinfo, Sym, mt, P);
    auto output = Output<S>(step.runtype, iterinfo, stats, P);
    // If calc0=true, a calculation of TD quantities is performed before starting the NRG iteration.
    if (step.nrg() && P.calc0 && !P.ZBW)
      docalc0(step, iterinfo, diag0, stats, output, oprecalc, Sym, mt, P);
    auto diag = P.ZBW ? nrg_ZBW(step, iterinfo, stats, diag0, output, dm, oprecalc, Sym, mt, P)
      : nrg_loop(step, iterinfo, coef, stats, diag0, output, dm, oprecalc, Sym, mpi, mt, P);
    fmt::print(fmt::emphasis::bold | fg(fmt::color::red), FMT_STRING("\nTotal energy: {:.18}\n"), stats.total_energy);
    stats.GS_energy = stats.total_energy;
    if (step.nrg() && P.dumpsubspaces) dm.dump_subspaces();
    fmt::print("\n** Iteration completed.\n\n");
    return diag;
  }
  NRG_calculation(MPI &mpi, const Workdir &workdir, const bool embedded) : mpi(mpi), P("param", "param", workdir, embedded), stats(P) {
    auto [diag0, iterinfo, coef, Sym] = read_data<S>(P, stats);
    Step step{P, RUNTYPE::NRG};
    AllSteps<S> dm(P.Ninit, P.Nlen);
    auto diag = run_nrg(step, iterinfo, coef, stats, diag0, dm, Sym);
    if (std::string(P.stopafter) == "nrg") exit1("*** Stopped after the first sweep.");
    dm.shift_abs_energies(stats.GS_energy); // we call this here, to enable a file dump
    if (P.dumpabsenergies)
      dm.dump_all_absolute_energies();
    if (P.dm) {
      if (P.need_rho()) {
        auto rho = init_rho(step, diag, Sym);
        rho.save(step.lastndx(), P, fn_rho);
        if (!P.ZBW) calc_densitymatrix(rho, dm, Sym, mt, P);
      }
      if (P.need_rhoFDM()) {
        calc_ZnD(dm, stats, Sym, P.T);
        if (P.logletter('w')) 
          report_ZnD(stats, P);
        fdm_thermodynamics(dm, stats, Sym, P.T);
        auto rhoFDM = init_rho_FDM(step.lastndx(), dm, stats, Sym, P.T);
        rhoFDM.save(step.lastndx(), P, fn_rhoFDM);
        if (!P.ZBW) calc_fulldensitymatrix(step, rhoFDM, dm, stats, Sym, mt, P);
      }
      if (std::string(P.stopafter) == "rho") exit1("*** Stopped after the DM calculation.");
      auto [diag0_dm, iterinfo_dm, coef_dm, Sym_dm] = read_data<S>(P, stats);
      Step step_dmnrg{P, RUNTYPE::DMNRG};
      run_nrg(step_dmnrg, iterinfo_dm, coef_dm, stats, diag0_dm, dm, Sym_dm);
      my_assert(num_equal(stats.GS_energy, stats.total_energy));
    }
  }
  ~NRG_calculation() {
    if (!P.embedded) mt.report(); // only when running as a stand-alone application
    if (P.done) { std::ofstream D("DONE"); } // Indicate completion by creating a flag file
  }
};

#endif
