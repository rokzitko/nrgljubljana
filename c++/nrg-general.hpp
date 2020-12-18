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

// Explicit request to stop at some calculation stage
inline void exit1(const std::string &message) {
  std::cout << std::endl << message << std::endl;
  exit(1);
}

template <scalar S> class NRG_calculation {
private:
  MPI_diag mpi;
  Params P;
  InputData<S> input;
  std::shared_ptr<Symmetry<S>> Sym;
  Stats<S> stats;
  Store<S> store;
  MemTime mt; // memory and timing statistics
public:
  auto run_nrg(const RUNTYPE runtype, Operators<S> operators, const Coef<S> &coef, DiagInfo<S> diag0) {
    Step step{P, runtype};
    diag0.states_report(Sym->multfnc());
    auto oprecalc = Oprecalc<S>(step.get_runtype(), operators, Sym, mt, P);
    auto output = Output<S>(step.get_runtype(), operators, stats, P);
    if (P.h5raw) {
       diag0.h5save(*output.h5raw, std::to_string(step.ndx()) + "/eigen/");
       operators.h5save(*output.h5raw, std::to_string(step.ndx()));
    }
    // If calc0=true, a calculation of TD quantities is performed before starting the NRG iteration.
    if (step.nrg() && P.calc0 && !P.ZBW())
      docalc0(step, operators, diag0, stats, output, oprecalc, Sym.get(), mt, P);
    auto diag = P.ZBW() ? nrg_ZBW(step, operators, stats, diag0, output, store, oprecalc, Sym.get(), mt, P)
                        : nrg_loop(step, operators, coef, stats, diag0, output, store, oprecalc, Sym.get(), mpi, mt, P);
    fmt::print(fmt::emphasis::bold | fg(fmt::color::red), FMT_STRING("\nTotal energy: {:.18}\n"), stats.total_energy);
    stats.GS_energy = stats.total_energy;
    if (step.nrg()) {
      store.shift_abs_energies(stats.GS_energy);
      if (P.dumpabsenergies) store.dump_all_absolute_energies();
      if (P.dumpsubspaces) store.dump_subspaces();
    }
    if (P.h5raw) store.h5save(*output.h5raw, "/store");
    fmt::print("\n** Iteration completed.\n\n");
    return diag;
  }
  void calc_rho(const DiagInfo<S> &diag) {
    Step step{P, RUNTYPE::NRG};
    step.set_last();
    auto rho = init_rho(step, diag, Sym->multfnc());
    rho.save(step.lastndx(), P, fn_rho);
    if (!P.ZBW()) calc_densitymatrix(rho, store, Sym.get(), mt, P);
  }
  void calc_rhoFDM() {
    Step step{P, RUNTYPE::NRG};
    step.set_last();
    calc_ZnD(store, stats, Sym.get(), P.T);
    if (P.logletter('w'))
      report_ZnD(stats, P);
    fdm_thermodynamics(store, stats, Sym.get(), P.T);
    auto rhoFDM = init_rho_FDM(step.lastndx(), store, stats, Sym->multfnc(), P.T);
    rhoFDM.save(step.lastndx(), P, fn_rhoFDM);
    if (!P.ZBW()) calc_fulldensitymatrix(step, rhoFDM, store, stats, Sym.get(), mt, P);
  }
  NRG_calculation(MPI_diag &mpi, std::unique_ptr<Workdir> workdir, const bool embedded) : 
    mpi(mpi), P("param", "param", std::move(workdir), embedded), input(P, "data"), Sym(input.Sym),
    stats(P, Sym->get_td_fields(), input.GS_energy), store(P.Ninit, P.Nlen)
  {
    auto diag = run_nrg(RUNTYPE::NRG, input.operators, input.coef, input.diag);
    if (P.dm) {
      if (P.need_rho()) calc_rho(diag);
      if (P.need_rhoFDM()) calc_rhoFDM();
      run_nrg(RUNTYPE::DMNRG, input.operators, input.coef, input.diag);
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
