#ifndef _subspaces_hpp_
#define _subspaces_hpp_

#include <memory>
#include <vector>
#include <h5cpp/all>
#include <range/v3/all.hpp>
#include "eigen.hpp"
#include "invar.hpp"

namespace NRG {

template<typename S> class Symmetry;

// Dimensions of the invariant subspaces |r,1>, |r,2>, |r,3>, etc. 
class SubspaceDimensions {
 private:
   std::vector<size_t> dims;
   std::vector<Invar> ancestors;
 public:
   SubspaceDimensions() = default;
   template<typename S>
     SubspaceDimensions(const Invar &I, const InvarVec &ancestors, const DiagInfo<S> &diagprev, std::shared_ptr<Symmetry<S>> Sym);
   [[nodiscard]] auto combs() const { return dims.size(); } // number of subspaces
   [[nodiscard]] auto rmax(const size_t i) const { // subspace dimension
     my_assert(i < combs());
     return dims[i];
   }
   [[nodiscard]] auto operator[](const size_t i) const { return rmax(i); }
   [[nodiscard]] auto exists(const size_t i) const {
     my_assert(i < combs());
     return dims[i] > 0;
   }
   [[nodiscard]] auto offset(const size_t i) const { // offset in the Hamiltonian matrix
     my_assert(i < combs());
     return ranges::accumulate(std::begin(dims), std::begin(dims) + i, size_t{0});
   }
   [[nodiscard]] auto chunk(const size_t i1) const {
     return std::make_pair(offset(i1-1), rmax(i1-1));
   }
   [[nodiscard]] auto view(const size_t i) const { // index range in the Hamiltonian matrix
     return boost::irange(offset(i), offset(i)+rmax(i));
   }
   [[nodiscard]] auto view_mma(const size_t i) const {
     return view(i-1); // Mathematica uses 1-based indexing
   }
   [[nodiscard]] auto uboost_view(const size_t i) const {
     return ublas::range(offset(i), offset(i)+rmax(i));
   }
   [[nodiscard]] auto uboost_view_mma(const size_t i) const {
     return uboost_view(i-1); // Mathematica uses 1-based indexing
   }
   [[nodiscard]] auto total() const { return ranges::accumulate(dims, 0); } // total number of states
   // *** Mathematica interfacing: i1,j1 are 1-based
   [[nodiscard]] bool offdiag_contributes(const size_t i1, const size_t j1) const { // i,j are 1-based (Mathematica interface)
     my_assert(1 <= i1 && i1 <= combs() && 1 <= j1 && j1 <= combs());
     my_assert(i1 != j1);
     return exists(i1-1) && exists(j1-1); // shift by 1
   }
   [[nodiscard]] Invar ancestor(const size_t i) const { return ancestors[i]; }
   void h5save(h5::fd_t &fd, const std::string &name) const {
     h5::write(fd, name + "/dims", dims);
     std::vector<std::string> ancestor_names;
     std::transform(ancestors.begin(), ancestors.end(), std::back_inserter(ancestor_names), [](const auto &I){ return I.name(); });
     h5::write(fd, name + "/ancestors", ancestor_names);
   }
 private:
   friend std::ostream &operator<<(std::ostream &os, const SubspaceDimensions &rmax) {
     for (const auto &x : rmax.dims) os << x << ' ';
     return os;
   }
   template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &dims; }
   friend class boost::serialization::access;
};

class SubspaceStructure : public std::map<Invar, SubspaceDimensions> {
 public:
   SubspaceStructure() = default;
   template<typename S> SubspaceStructure(const DiagInfo<S> &, std::shared_ptr<Symmetry<S>>);
   // List of invariant subspaces in which diagonalisations need to be performed
   [[nodiscard]] std::vector<Invar> task_list() const {
     std::vector<std::pair<size_t, Invar>> tasks_with_sizes;
     for (const auto &[I, rm] : *this)
       if (rm.total())
         tasks_with_sizes.emplace_back(rm.total(), I);
     ranges::sort(tasks_with_sizes, std::greater<>()); // sort in the *decreasing* order!
     auto nr       = tasks_with_sizes.size();
     auto min_size = tasks_with_sizes.back().first;
     auto max_size = tasks_with_sizes.front().first;
     std::cout << "Stats: nr=" << nr << " min=" << min_size << " max=" << max_size << std::endl;
     return tasks_with_sizes | ranges::views::transform( [](const auto &p) { return p.second; } ) | ranges::to<std::vector>();
   }
   void dump() const {
     for(const auto &[I, rm]: *this)
       std::cout << "rmaxvals(" << I << ")=" << rm << " total=" << rm.total() << std::endl;
   }
   [[nodiscard]] auto at_or_null(const Invar &I) const {
     const auto i = this->find(I);
     return i == this->cend() ? SubspaceDimensions() : i->second;
   }
   void h5save(h5::fd_t &fd, const std::string &name) const {
     for (const auto &[I, rm]: *this)
       rm.h5save(fd, name + "/" + I.name());
   }
};

} // namespace

#endif
