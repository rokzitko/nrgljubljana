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
#include "h5.hpp"

namespace NRG {

template<scalar S>
struct ThermoSub {
  StoredEigen<S> eig;
  [[nodiscard]] auto kept() const { return eig.kept(); }
  [[nodiscard]] auto total() const { return eig.total(); }
  [[nodiscard]] auto min() const { return eig.min(); }
  [[nodiscard]] auto max() const { return eig.max(); }
  [[nodiscard]] auto all() const { return eig.all(); }
  void h5save(H5Easy::File &fd, const std::string &name) const {
    h5_dump_scalar(fd, name + "/kept", kept());
    h5_dump_scalar(fd, name + "/total", total());
    h5_dump_scalar(fd, name + "/min", min());
    h5_dump_scalar(fd, name + "/max", max());
  }
};

template<scalar S>
class ThermoSubs : public std::map<Invar, ThermoSub<S>> {
 public:
    ThermoSubs() = default;
    explicit ThermoSubs(const DiagInfo<S> &diag, const bool last) {
      for (const auto &[I, eig]: diag)
        (*this)[I] = { StoredEigen<S>(eig, last) };
    }
    void h5save(H5Easy::File &fd, const std::string &name) const {
      const std::vector<int> dummy = {1};
      for (const auto &I : *this | boost::adaptors::map_keys)
        H5Easy::dump(fd, name+"/list/" + I.name(), dummy);
     for (const auto &[I, sub]: *this)
       sub.h5save(fd, name + "/data/" + I.name());
    }
};

template<scalar S>
class ThermoStore : public std::vector<ThermoSubs<S>> {
 public:
    const size_t Nbegin, Nend; // range of valid indexes
    ThermoStore(const size_t Nbegin, const size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend); }
    auto Nall() const { return boost::irange(Nbegin, Nend); }
    void dump_abs_G(std::ostream &F) const {
      for (const auto N : Nall()) {
        F << std::endl << "===== Iteration number: " << N+1 << std::endl; // mind the shift by 1
        for (const auto &[I, ds]: this->at(N))
         F << "Subspace: " << I << std::endl << (ds.eig.values.all_abs_G() | ranges::to_vector) << std::endl;
     }
   }
   void dump_all_absolute_energies(const std::string &filename = "absolute_energies.dat"s) {
     std::ofstream F(filename);
     this->dump_abs_G(F);
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
   void h5save(H5Easy::File &fd, const std::string &name) const {
     const std::vector range = {Nbegin, Nend};
     H5Easy::dump(fd, name + "/range", range);
     for (const auto N : Nall())
        this->at(N).h5save(fd, name + "/" + std::to_string(N+1)); // note the shift by 1
    }
};

struct BackiterSub {
  SubspaceDimensions rmax;
  size_t nrkept = 0;
  size_t dim = 0;
  [[nodiscard]] auto kept() const { return nrkept; }
  [[nodiscard]] auto total() const { return dim; }
};

class BackiterSubs : public std::map<Invar, BackiterSub> {
 public:
   BackiterSubs() = default;
   template<scalar S>
   BackiterSubs(const DiagInfo<S> &diag, const SubspaceStructure &substruct) {
     for (const auto &[I, eig] : diag)
       (*this)[I] = { substruct.at_or_null(I), eig.getnrkept(), eig.getdim() };
   }
};

class BackiterStore : public std::vector<BackiterSubs> {
 public:
   const size_t Nbegin, Nend;
   BackiterStore(const size_t Nbegin, const size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend); }
   auto Nall() const { return boost::irange(Nbegin, Nend); }
};

} // namespace

#endif
