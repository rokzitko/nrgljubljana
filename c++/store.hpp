#ifndef _store_hpp_
#define _store_hpp_

#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>
#include "invar.hpp"
#include "eigen.hpp"
#include "subspaces.hpp"

namespace NRG {

// Container for all information which needs to be gathered in each invariant subspace.
// Required for the density-matrix construction.
template<typename S>
struct Sub {
  Eigen<S> eig;
  SubspaceDimensions rmax;
  bool is_last = false;
  [[nodiscard]] auto kept() const { return eig.getnrkept(); }
  [[nodiscard]] auto total() const { return eig.getdim(); }
  [[nodiscard]] auto min() const { return is_last ? 0 : kept(); } // min(), max() return the range of D states to be summed over in FDM
  [[nodiscard]] auto max() const { return total(); }
  [[nodiscard]] auto all() const { return boost::irange(min(), max()); }
};

template<typename S> 
class Subs : public std::map<Invar, Sub<S>> {
 public:
   Subs() = default;
   Subs(const DiagInfo<S> &diag, const SubspaceStructure &substruct, const bool last) {
     for (const auto &[I, eig]: diag)
       (*this)[I] = { eig, substruct.at_or_null(I), last };
   }
};

template<typename S>
class Store : public std::vector<Subs<S>> {
 public:
   const size_t Nbegin, Nend; // range of valid indexes
   Store(const size_t Nbegin, const size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend ? Nend : 1); } // at least 1 for ZBW
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
       for (const auto &[I, sub] : this->at(N))
         O << "I=" << I << " kept=" << sub.kept() << " total=" << sub.total() << std::endl;
       O << std::endl;
     }
   }
   void shift_abs_energies(const double GS_energy) {
     for (const auto N : Nall())
       for (auto &ds : this->at(N) | boost::adaptors::map_values)
         ds.eig.subtract_GS_energy(GS_energy);
   }
};

} // namespace

#endif
