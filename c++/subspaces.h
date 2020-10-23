#ifndef _subspaces_h_
#define _subspaces_h_

#include <memory>
#include <vector>
#include "eigen.h"
#include <range/v3/all.hpp>

template<typename S> class Symmetry;

// Dimensions of the invariant subspaces |r,1>, |r,2>, |r,3>, etc. The name "rmax" comes from the maximal value of
// the index "r" which ranges from 1 through rmax.

class Rmaxvals {
 private:
   std::vector<size_t> values;
 public:
   Rmaxvals() = default;
   template<typename S>
     Rmaxvals(const Invar &I, const InvarVec &In, const DiagInfo<S> &diagprev, std::shared_ptr<Symmetry<S>> Sym);
   auto combs() const { return values.size(); }
   auto rmax(const size_t i) const {
     my_assert(i < combs());
     return values[i];
   }
   auto exists(const size_t i) const {
     my_assert(i < combs());
     return values[i] > 0;
   }
   auto offset(const size_t i) const {
     my_assert(i < combs());
     return ranges::accumulate(std::begin(values), std::begin(values) + i, size_t{0});
   }
   auto view(const size_t i) const {
     return boost::irange(offset(i), offset(i)+rmax(i));
   }
   auto view_mma(const size_t i) const {
     return view(i-1); // Mathematica uses 1-based indexing
   }
   auto ubview(const size_t i) const {
     return ublas::range(offset(i), offset(i)+rmax(i));
   }
   auto ubview_mma(const size_t i) const {
     return ubview(i-1); // Mathematica uses 1-based indexing
   }
   auto operator[](const size_t i) const { return rmax(i); }
   auto total() const { return ranges::accumulate(values, 0); } // total number of states
   // *** Mathematica interfacing: i1,j1 are 1-based
   bool offdiag_contributes(const size_t i1, const size_t j1) const { // i,j are 1-based (Mathematica interface)
     my_assert(1 <= i1 && i1 <= combs() && 1 <= j1 && j1 <= combs());
     my_assert(i1 != j1);
     return exists(i1-1) && exists(j1-1); // shift by 1
   }
   auto chunk(const size_t i1) const {
     return std::make_pair(offset(i1-1), rmax(i1-1));
   }
 private:
   friend std::ostream &operator<<(std::ostream &os, const Rmaxvals &rmax) {
     for (const auto &x : rmax.values) os << x << ' ';
     return os;
   }
   template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &values; }
   friend class boost::serialization::access;
};

class QSrmax : public std::map<Invar, Rmaxvals> {
 public:
   QSrmax() {}
   template<typename S> QSrmax(const DiagInfo<S> &, std::shared_ptr<Symmetry<S>>);
   // List of invariant subspaces in which diagonalisations need to be performed
   std::vector<Invar> task_list() const {
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
   auto at_or_null(const Invar &I) const {
     const auto i = this->find(I);
     return i == this->cend() ? Rmaxvals() : i->second;
   }
};

#endif
