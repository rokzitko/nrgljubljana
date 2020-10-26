#ifndef _step_hpp_
#define _step_hpp_

#include <algorithm>
#include "params.hpp"
#include "numerics.hpp"

#include <fmt/format.h>

namespace NRG {

class Step {
 private:
   // N denotes the order of the Hamiltonian. N=0 corresponds to H_0, i.e. the initial Hamiltonian
   int trueN = 0; // "true N", sets the energy scale, it may be negative, trueN <= ndxN
   size_t ndxN = 0; // "index N", iteration step, used as an array index, ndxN >= 0
   const Params &P; // reference to parameters (beta, T)
   const RUNTYPE runtype; // NRG vs. DM-NRG run
 public:
   void set(const int newN) {
     trueN = newN;
     ndxN = std::max(newN, 0);
   }
   void init() { set(P.Ninit); }
   Step(const Params &P_, const RUNTYPE runtype_) : P(P_), runtype(runtype_) { init(); }
   void next() { trueN++; ndxN++; }
   [[nodiscard]] auto N() const { return ndxN; }
   [[nodiscard]] auto ndx() const { return ndxN; }
   [[nodiscard]] auto energyscale() const { return P.SCALE(trueN+1); } // current energy scale in units of bandwidth D
   [[nodiscard]] auto scale() const { // scale factor as used in the calculation
     return P.absolute ? 1.0 : energyscale();
   }
   [[nodiscard]] auto unscale() const { // 'unscale' parameter for dimensionless quantities
     return P.absolute ? energyscale() : 1.0;
   }
   [[nodiscard]] auto Teff() const { return energyscale()/P.betabar; }  // effective temperature for thermodynamic calculations
   [[nodiscard]] auto TD_factor() const { return P.betabar / unscale(); }
   [[nodiscard]] auto scT() const { return scale()/P.T; } // scT = scale*P.T, scaled physical temperature that appears in the exponents in spectral function calculations (Boltzmann weights)
   [[nodiscard]] std::pair<size_t, size_t> NM() const {
     const size_t N = ndxN / P.channels;
     const size_t M = ndxN - N*P.channels; // M ranges 0..channels-1
     return {N, M};
   }
   void infostring() const {
     auto info = fmt::format(" ***** [{}] Iteration {}/{} (scale {}) ***** ", runtype == RUNTYPE::NRG ? "NRG"s : "DM"s, 
                             ndxN+1, int(P.Nmax), energyscale());
     info += P.substeps ? fmt::format(" step {} substep {}", NM().first+1, NM().second+1) : "";
     fmt::print(fmt::emphasis::bold, "\n{}\n", info);
   }
   void set_ZBW() {
     trueN = P.Ninit - 1; // if Ninit=0, trueN will be -1 (this is the only exceptional case)
     ndxN = P.Ninit;
   }
   // Return true if the spectral-function merging is to be performed at the current step
   [[nodiscard]] auto N_for_merging() const {
     if (P.NN1) return true;
     if (P.NN2avg) return true;
     return P.NN2even ? is_even(ndxN) : is_odd(ndxN);
   }
   [[nodiscard]] size_t firstndx() const { return P.Ninit; }
   [[nodiscard]] size_t lastndx() const { return P.ZBW ? P.Ninit : P.Nmax-1; }
   // Return true if this is the first step of the NRG iteration
   [[nodiscard]] auto first() const { return ndxN == firstndx(); }
   // Return true if N is the last step of the NRG iteration
   [[nodiscard]] auto last(int N) const {
     return N == lastndx() || (P.ZBW && N == firstndx()); // special case!
   }
   [[nodiscard]] auto last() const { return last(ndxN); }
   [[nodiscard]] auto end() const { return ndxN >= P.Nmax; } // ndxN is outside the allowed range
   // NOTE: for ZBWcalculations, Ninit=0 and Nmax=0, so that first() == true and last() == true for ndxN=0.
   [[nodiscard]] auto nrg() const { return runtype == RUNTYPE::NRG; }
   [[nodiscard]] auto dmnrg() const { return runtype == RUNTYPE::DMNRG; }
   // Index 'n' of the last site in the existing chain, f_n (at iteration 'N'). The site being added is f_{n+1}. This
   // is the value that we use in building the matrix.
   [[nodiscard]] auto getnn() const { return ndxN; }
   [[nodiscard]] auto get_runtype() const { return runtype; }
};

} // namespace

#endif
