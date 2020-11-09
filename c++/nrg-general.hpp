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

#ifndef _nrg_general_hpp_
#define _nrg_general_hpp_

#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <complex>
#include <numeric>
#include <limits>
#include <memory>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <list>
#include <deque>
#include <set>
#include <stdexcept>

// C headers
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstring>
#include <unistd.h>

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/optional.hpp>

// ublas matrix & vector containers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>

// This is included in the library only. Should not be used if a cblas library is available.
#ifdef CBLAS_WORKAROUND
 #define ADD_
 #include "cblas_globals.c"
 #include "cblas_dgemm.c"
 #include "cblas_zgemm.c"
 #include "cblas_xerbla.c"
#endif

// Serialization support (used for storing to files and for MPI)
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

// MPI support
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>

#include <range/v3/all.hpp>

#include "nrg-lib.hpp" // exposed interfaces for wrapping into a library
#include "portabil.hpp"
#include "debug.hpp"
#include "basicio.hpp"
#include "misc.hpp"
#include "openmp.hpp"
#include "mp.hpp"
#include "traits.hpp"
#include "workdir.hpp"
#include "params.hpp"
#include "numerics.hpp"
#include "io.hpp"
#include "time_mem.hpp"
#include "outfield.hpp"
#include "core.hpp"
#include "mk_sym.hpp"

namespace NRG {

using namespace std::string_literals;
using namespace fmt::literals;
using namespace boost::numeric;
using namespace boost::numeric::ublas; // keep this!
namespace atlas = boost::numeric::bindings::atlas;

template <typename S> class NRG_calculation {
private:
  MPI_diag mpi;
  Params P;
  Stats<S> stats;
  MemTime mt; // memory and timing statistics
public:
  auto run_nrg(Step &step, Operators<S> &operators, const Coef<S> &coef, Stats<S> &stats, const DiagInfo<S> &diag0,
               Store<S> &store, std::shared_ptr<Symmetry<S>> Sym) {
    diag0.states_report(Sym->multfnc());
    auto oprecalc = Oprecalc<S>(step.get_runtype(), operators, Sym, mt, P);
    auto output = Output<S>(step.get_runtype(), operators, stats, P);
    if (output.h5raw) {
       diag0.h5save(output.h5raw.value(), std::to_string(step.ndx()) + "/eigen/", step.dmnrg());
       operators.h5save(output.h5raw.value(), std::to_string(step.ndx()));
    }
    // If calc0=true, a calculation of TD quantities is performed before starting the NRG iteration.
    if (step.nrg() && P.calc0 && !P.ZBW)
      docalc0(step, operators, diag0, stats, output, oprecalc, Sym.get(), mt, P);
    auto diag = P.ZBW ? nrg_ZBW(step, operators, stats, diag0, output, store, oprecalc, Sym.get(), mt, P)
                      : nrg_loop(step, operators, coef, stats, diag0, output, store, oprecalc, Sym.get(), mpi, mt, P);
    fmt::print(fmt::emphasis::bold | fg(fmt::color::red), FMT_STRING("\nTotal energy: {:.18}\n"), stats.total_energy);
    stats.GS_energy = stats.total_energy;
    if (step.nrg() && P.dumpsubspaces) store.dump_subspaces();
    if (output.h5raw) store.h5save(output.h5raw.value(), "/store");
    fmt::print("\n** Iteration completed.\n\n");
    return diag;
  }
  NRG_calculation(MPI_diag &mpi, const Workdir &workdir, const bool embedded) : mpi(mpi), P("param", "param", workdir, embedded), stats(P) {
    auto [diag0, operators, coef, Sym] = read_data<S>(P, stats);
    Step step{P, RUNTYPE::NRG};
    Store<S> store(P.Ninit, P.Nlen);
    auto diag = run_nrg(step, operators, coef, stats, diag0, store, Sym);
    if (std::string(P.stopafter) == "nrg") exit1("*** Stopped after the first sweep.");
    store.shift_abs_energies(stats.GS_energy); // we call this here, to enable a file dump
    if (P.dumpabsenergies)
      store.dump_all_absolute_energies();
    if (P.dm) {
      if (P.need_rho()) {
        auto rho = init_rho(step, diag, Sym->multfnc());
        rho.save(step.lastndx(), P, fn_rho);
        if (!P.ZBW) calc_densitymatrix(rho, store, Sym.get(), mt, P);
      }
      if (P.need_rhoFDM()) {
        calc_ZnD(store, stats, Sym.get(), P.T);
        if (P.logletter('w'))
          report_ZnD(stats, P);
        fdm_thermodynamics(store, stats, Sym.get(), P.T);
        auto rhoFDM = init_rho_FDM(step.lastndx(), store, stats, Sym->multfnc(), P.T);
        rhoFDM.save(step.lastndx(), P, fn_rhoFDM);
        if (!P.ZBW) calc_fulldensitymatrix(step, rhoFDM, store, stats, Sym.get(), mt, P);
      }
      if (std::string(P.stopafter) == "rho") exit1("*** Stopped after the DM calculation.");
      auto [diag0_dm, operators_dm, coef_dm, Sym_dm] = read_data<S>(P, stats);
      Step step_dmnrg{P, RUNTYPE::DMNRG};
      run_nrg(step_dmnrg, operators_dm, coef_dm, stats, diag0_dm, store, Sym_dm);
      my_assert(num_equal(stats.GS_energy, stats.total_energy));
    }
  }
  NRG_calculation(const NRG_calculation &) = delete;
  NRG_calculation(NRG_calculation &&) = delete;
  NRG_calculation & operator=(const NRG_calculation &) = delete;
  NRG_calculation & operator=(const NRG_calculation &&) = delete;
  ~NRG_calculation() {
    if (!P.embedded) mt.report(); // only when running as a stand-alone application
    if (P.done) { std::ofstream D("DONE"); } // Indicate completion by creating a flag file
  }
};

} // namespace

#endif
