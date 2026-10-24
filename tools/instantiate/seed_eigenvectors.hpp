#ifndef _instantiate_seed_eigenvectors_hpp_
#define _instantiate_seed_eigenvectors_hpp_

#include <cmath>
#include <traits.hpp>

namespace NRG::Spawn {

inline void canonicalize_seed_eigenvectors(Matrix_traits<double> &vectors) {
  if (vectors.cols() == 0) return;
  for (::Eigen::Index row = 0; row < vectors.rows(); ++row) {
    ::Eigen::Index pivot = 0;
    for (::Eigen::Index col = 1; col < vectors.cols(); ++col)
      if (std::abs(vectors(row, col)) > std::abs(vectors(row, pivot))) pivot = col;
    // LAPACK leaves each row's overall sign arbitrary; the first largest component fixes it.
    const bool flip = vectors(row, pivot) < 0.0;
    for (::Eigen::Index col = 0; col < vectors.cols(); ++col) {
      auto &value = vectors(row, col);
      if (value == 0.0) value = 0.0;
      else if (flip) value = -value;
    }
  }
}

} // namespace NRG::Spawn

#endif
