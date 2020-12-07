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
 public:
   std::unordered_map<t_eigen, t_eigen> cluster_mapping;
   // Fix splittings of eigenvalues.
   void fix_it(DiagInfo<S> &diag) {
     for(auto &[I, eig]: diag) {
       for (auto &r : eig.value_zero)
         if (auto m = cluster_mapping.find(r); m != cluster_mapping.cend())
           r = m->second;
     }
   }
   // Find clusters of values which differ by at most 'epsilon'
   Clusters(DiagInfo<S> &diag, const double epsilon, bool fix = true) {
     const auto energies = diag.sorted_energies_rel_zero();
     my_assert(energies.size());
     auto e0 = energies.front();  // energy of the lower boundary of the cluster, [e0:e1]
     auto i0 = energies.cbegin(); // iterator to the lower boundary of the cluster, [i0:i1]
     int size = 1;                // number of states in the current cluster
     for (auto i = energies.begin(); i != energies.end(); ++i) {
       if ((*i - e0) < epsilon) { // in the cluster
         size++;
       } else { // end of cluster detected
         auto i1 = i;
         if (size > 1) {            // is this a real cluster?
           if (cluster_splitting(i0, i1)) { // are the states actually split?
             auto replace_with = *i0;    // use the lowest eigenvalue of the cluster
             for (auto j = (i0 + 1); j != i1; ++j) // skip 1st
               if (*j != *i0) cluster_mapping.insert({*j, replace_with});
           }
         }
         e0   = *i;
         i0   = i;
         size = 1;
       }
     }
     if (fix) fix_it(diag);
   }
};

} // namespace

#endif
