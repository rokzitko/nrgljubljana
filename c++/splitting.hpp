// Code for correcting floating-point roundoff errors
// Rok Zitko, rok.zitko@ijs.si

#ifndef _splitting_hpp_
#define _splitting_hpp_

#include <iostream>
#include <unordered_map>
#include "portabil.hpp"
#include "traits.hpp"
#include "eigen.hpp"

namespace NRG {

template<typename T>
inline void cluster_show(const T &i0, const T &i1) {
  std::cout << "[";
  for (auto j = i0; j != i1; ++j) { std::cout << HIGHPREC(*j) << " "; }
  std::cout << "]" << std::endl;
}

// Returns true if not all the states have the same energy.
template<typename T>
inline bool cluster_splitting(const T &i0, const T &i1) {
  my_assert(i0 != i1); // non-empty set
  // We need to compare all distinct pairs.
  for (auto i = i0; i != i1; ++i)
    for (auto j = i + 1; j != i1; ++j)
      if (*i != *j) return true;
  return false;
  }

template<scalar S, typename t_eigen = eigen_traits<S>>
class Clusters {
 private:
   const Params &P;
 public:
   std::unordered_map<t_eigen, t_eigen> cluster_mapping;
   // Fix splittings of eigenvalues.
   void fix_it(DiagInfo<S> &diag) {
     for(auto &[I, eig]: diag) {
       auto v = eig.values.all_rel_zero() | ranges::to_vector;
       for (auto &r : v)
         if (auto m = cluster_mapping.find(r); m != cluster_mapping.cend())
           r = m->second;
       eig.values.set_corr(std::move(v)); // ATTENTION: these are now zero-offset values!
       if (!P.floquet)
         eig.values.crit_copy_corr(); // update truncation criterion with corrected eigenvalues; zero-offset!
     }
   }
   // Find clusters of values which differ by at most 'epsilon'
    Clusters(DiagInfo<S> &diag, const double epsilon, const Params &_P, bool fix = true) : P(_P) {
      const auto energies = diag.sorted_energies_rel_zero();
      my_assert(energies.size());
      auto e0 = energies.front();  // energy of the lower boundary of the cluster, [e0:e1]
      auto i0 = energies.cbegin(); // iterator to the lower boundary of the cluster, [i0:i1]
      const auto process_cluster = [this, &i0](const auto i1) {
        if (std::distance(i0, i1) > 1 && cluster_splitting(i0, i1)) {
          const auto replace_with = *i0;
          for (auto j = i0 + 1; j != i1; ++j)
            if (*j != *i0) cluster_mapping.insert({*j, replace_with});
        }
      };
      for (auto i = energies.begin() + 1; i != energies.end(); ++i) {
        if ((*i - e0) >= epsilon) {
          process_cluster(i);
          e0   = *i;
          i0   = i;
        }
      }
      process_cluster(energies.cend());
      if (fix) fix_it(diag);
    }
};

} // namespace

#endif
