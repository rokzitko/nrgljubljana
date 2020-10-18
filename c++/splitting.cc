// Code for correcting floating-point roundoff errors
// Rok Zitko, rok.zitko@ijs.si

using mapdd = std::unordered_map<t_eigen, t_eigen>;

// Fix splittings of eigenvalues. Returns true if any changes had been made.
template<typename S>
void fix_splittings(DiagInfo_tmpl<S> &diag, const mapdd &cluster_mapping) {
  for(auto &[I, eig]: diag) { 
    for (auto &r : eig.value_zero) 
      if (auto m = cluster_mapping.find(r); m != cluster_mapping.cend())
        r = m->second;
  }
}

template<typename T>
  void cluster_show(const T &i0, const T &i1) {
    std::cout << "[";
    for (auto j = i0; j != i1; ++j) { std::cout << HIGHPREC(*j) << " "; }
    std::cout << "]" << std::endl;
  }

// Returns true if not all the states have the same energy.
template<typename T>
  bool cluster_splitting(const T &i0, const T &i1) {
    my_assert(i0 != i1); // non-empty set
    // We need to compare all distinct pairs.
    for (auto i = i0; i != i1; ++i)
      for (auto j = i + 1; j != i1; ++j)
        if (*i != *j) return true;
    return false;
  }

// Find clusters of values which differ by at most 'epsilon'
mapdd find_clusters(const std::vector<t_eigen> &energies, double epsilon) {
  my_assert(energies.size());
  mapdd cluster_mapping;
  auto e0 = energies[0];      // energy of the lower boundary of the cluster, [e0:e1]
  auto i0 = cbegin(energies); // iterator to the lower boundary of the cluster, [i0:i1]
  int size = 1;                // number of states in the current cluster
  for (auto i = begin(energies); i != end(energies); ++i) {
    if ((*i - e0) < epsilon) { // in the cluster
      size++;
    } else { // end of cluster detected
      auto i1 = i;
      if (size > 1) {            // is this a real cluster?
        if (cluster_splitting(i0, i1)) { // are the states actually split?
          auto replace_with = *i0;    // use the lowest eigenvalue of the cluster
          for (auto j = (i0 + 1); j != i1; ++j) // skip 1st
            if (*j != *i0) cluster_mapping.insert(make_pair(*j, replace_with));
        }
      }
      e0   = *i;
      i0   = i;
      size = 1;
    }
  }
  return cluster_mapping;
}
