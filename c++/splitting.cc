// Code for correcting floating-point roundoff errors
// Rok Zitko, rok.zitko@ijs.si

// Fix splittings of eigenvalues. Returns true if any changes had been made.
bool fix_splittings(DiagInfo &diag) {
  bool changes_made = false;
  LOOP(diag, is)
  for (auto &r : EIGEN(is).value) {
    const auto m = STAT::cluster_mapping.find(r);
    if (m != end(STAT::cluster_mapping)) {
      r            = m->second;
      changes_made = true;
    }
  }
  return changes_made;
}

// Iterator over eigenvalues
using svdi = STDEVEC::iterator;

void cluster_show(const svdi &i0, const svdi &i1) {
  cout << "[";
  for (svdi j = i0; j != i1; ++j) { cout << HIGHPREC(*j) << " "; }
  cout << "]" << endl;
}

// Returns true if not all the states have the same energy.
bool cluster_splitting(const svdi &i0, const svdi &i1) {
  my_assert(i0 != i1);
  // We need to compare all distinct pairs.
  for (svdi i = i0; i != i1; ++i)
    for (auto j = i + 1; j != i1; ++j)
      if (*i != *j) return true;
  return false;
}

// Find clusters of values which differ by at most 'epsilon'
void find_clusters(STDEVEC &energies, double epsilon, mapdd &cluster_mapping) {
  my_assert(energies.size() > 0);
  t_eigen e0       = energies[0];     // energy of the lower boundary of the cluster, [e0:e1]
  auto i0          = begin(energies); // iterator to the lower boundary of the cluster, [i0:i1]
  int cluster_size = 1;               // number of states in the current cluster

  for (auto i = begin(energies); i != end(energies); ++i) {
    if ((*i - e0) < epsilon) { // in the cluster
      cluster_size++;
    } else { // end of cluster detected
      auto i1 = i;
      if (logletter('X')) cluster_show(i0, i1);
      if (cluster_size > 1) {            // is this a real cluster?
        if (cluster_splitting(i0, i1)) { // are the states actually split?
          t_eigen replace_with = *i0;    // use the lowest eigenvalue of the cluster
          if (logletter('X')) cout << " -> " << setprecision(std::numeric_limits<double>::max_digits10) << replace_with << endl;
          for (auto j = (i0 + 1); j != i1; ++j) // skip 1st
            if (*j != *i0) cluster_mapping.insert(make_pair(*j, replace_with));
        }
      }
      e0           = *i;
      i0           = i;
      cluster_size = 1;
    }
  }
}
