#ifndef _subspaces_hpp_
#define _subspaces_hpp_

#include <memory>
#include <vector>

#include <range/v3/all.hpp>

#include "eigen.hpp"
#include "invar.hpp"
#include "h5.hpp"

namespace NRG {

template<scalar S> class Symmetry;

// Dimensions of the invariant subspaces |r,1>, |r,2>, |r,3>, etc. 
class SubspaceDimensions {
 private:
   std::vector<size_t> dims;
   std::vector<Invar> ancestors;
 public:
   SubspaceDimensions() = default;
   template<scalar S>
     SubspaceDimensions(const Invar &I, const InvarVec &ancestors, const DiagInfo<S> &diagprev, 
                        const Symmetry<S> *Sym, const bool ignore_inequality = false);
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
   [[nodiscard]] auto view(const size_t i) const { // index range in the Hamiltonian matrix // XXX: rename to range
     return boost::irange(offset(i), offset(i)+rmax(i));
   }
   [[nodiscard]] auto view_mma(const size_t i) const {
     return view(i-1); // Mathematica uses 1-based indexing
   }
   [[nodiscard]] auto part(const size_t i) const {
     return std::make_pair(offset(i), offset(i)+rmax(i));
   }
   [[nodiscard]] auto part_mma(const size_t i1) const {
     return part(i1-1); // Mathematica uses 1-based indexing
   }
   [[nodiscard]] auto total() const { return ranges::accumulate(dims, 0); } // total number of states
   // *** Mathematica interfacing: i1,j1 are 1-based
   [[nodiscard]] bool offdiag_contributes(const size_t i1, const size_t j1) const { // i,j are 1-based (Mathematica interface)
     my_assert(1 <= i1 && i1 <= combs() && 1 <= j1 && j1 <= combs());
     my_assert(i1 != j1);
     return exists(i1-1) && exists(j1-1); // shift by 1
   }
   [[nodiscard]] Invar ancestor(const size_t i) const { return ancestors[i]; }
   void dump(std::ostream &F = std::cout) const {
     my_assert(dims.size() == ancestors.size());
     for (int i = 0; i < dims.size(); i++)
       F << "[" << ancestors[i] << "] total=" << dims[i] << std::endl;
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
     std::vector<std::string> ancestor_names;
     for (const auto i : range0(combs()))
       if (dims[i]) ancestor_names.push_back(ancestors[i].name()); // only true ancestors with dim>0
     H5Easy::dump(fd, name + "/ancestors", ancestor_names);
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
   template<scalar S> SubspaceStructure(const DiagInfo<S> &, const Symmetry<S> *);
   void dump(std::ostream &F = std::cout) const {
     for(const auto &[I, rm]: *this)
       F << "rmaxvals(" << I << ")=" << rm << " total=" << rm.total() << std::endl;
   }
   [[nodiscard]] auto at_or_null(const Invar &I) const {
     const auto i = this->find(I);
     return i == this->cend() ? SubspaceDimensions() : i->second;
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
     for (const auto &[I, rm]: *this)
       rm.h5save(fd, name + "/" + I.name());
   }
};

// List of invariant subspaces in which diagonalisations need to be performed
class TaskList {
 private:
   std::vector<std::pair<size_t, Invar>> tasks_with_sizes;
   std::vector<Invar> tasks;
 public:
   void stats(std::ostream &F) {
     auto nr       = tasks_with_sizes.size();
     auto min_size = tasks_with_sizes.back().first;
     auto max_size = tasks_with_sizes.front().first;
     F << "Stats: nr=" << nr << " min=" << min_size << " max=" << max_size << std::endl;
   }
   explicit TaskList(const SubspaceStructure &structure, std::ostream &F = std::cout) {
     for (const auto &[I, rm] : structure)
       if (rm.total())
         tasks_with_sizes.emplace_back(rm.total(), I);
     ranges::sort(tasks_with_sizes, std::greater<>()); // sort in the *decreasing* order!
     stats(F);
     tasks = tasks_with_sizes | ranges::views::transform( [](const auto &p) { return p.second; } ) | ranges::to<std::vector>();
   }
   [[nodiscard]] std::vector<Invar> get() const { return tasks; }
};

} // namespace

#endif
