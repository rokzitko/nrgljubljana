#ifndef _store_hpp_
#define _store_hpp_

#include "invar.hpp"
#include "eigen.hpp"
#include "subspaces.hpp"

namespace NRG {

// Information about the number of states, kept and discarded, rmax, and eigenenergies. Required for the
// density-matrix construction.
template<typename S>
struct DimSub {
  size_t kept  = 0;
  size_t total = 0;
  Rmaxvals rmax;
  Eigen<S> eig;
  bool is_last = false;
  auto min() const { return is_last ? 0 : kept; } // min(), max() return the range of D states to be summed over in FDM
  auto max() const { return total; }
  auto all() const { return boost::irange(min(), max()); }
};

// Full information about the number of states and matrix dimensions
// Example: dm[N].rmax[I] etc.
template<typename S>
using Subs = std::map<Invar, DimSub<S>>;

template<typename S>
class AllSteps : public std::vector<Subs<S>> {
 public:
   const size_t Nbegin, Nend; // range of valid indexes
   AllSteps(const size_t Nbegin, const size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend ? Nend : 1); } // at least 1 for ZBW
   auto Nall() const { return boost::irange(Nbegin, Nend); }
   void dump_absenergyG(std::ostream &F) const {
     for (const auto N : Nall()) {
       F << std::endl << "===== Iteration number: " << N << std::endl;
       for (const auto &[I, ds]: this->at(N))
         F << "Subspace: " << I << std::endl << ds.eig.absenergyG << std::endl;
     }
   }
   void dump_all_absolute_energies(const std::string &filename = "absolute_energies.dat"s) {
     std::ofstream F(filename);
     this->dump_absenergyG(F);
   }
   // Save a dump of all subspaces, with dimension info, etc.
   void dump_subspaces(const std::string &filename = "subspaces.dat"s) const {
     std::ofstream O(filename);
     for (const auto N : Nall()) {
       O << "Iteration " << N << std::endl;
       O << "len_dm=" << this->at(N).size() << std::endl;
       for (const auto &[I, DS] : this->at(N))
         O << "I=" << I << " kept=" << DS.kept << " total=" << DS.total << std::endl;
       O << std::endl;
     }
   }
   void shift_abs_energies(const double GS_energy) {
     for (const auto N : Nall())
       for (auto &ds : this->at(N) | boost::adaptors::map_values)
         ds.eig.subtract_GS_energy(GS_energy);
   }
   void store(const size_t ndx, const DiagInfo<S> &diag, const QSrmax &qsrmax, const bool last) {
     my_assert(Nbegin <= ndx && ndx < Nend);
     for (const auto &[I, eig]: diag)
       (*this)[ndx][I] = { eig.getnrkept(), eig.getdim(), qsrmax.at_or_null(I), eig, last };
   }
};

} // namespace

#endif
