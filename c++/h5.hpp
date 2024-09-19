#ifndef _h5_hpp_
#define _h5_hpp_

#include <iostream>
#include <vector>
#include <string>
#include <complex>
#include <type_traits>

#include "traits.hpp"

#define H5_USE_EIGEN
#include <highfive/highfive.hpp>
#include <highfive/eigen.hpp>
#include <highfive/H5Easy.hpp>

namespace NRG {
   template<typename T>
   void h5_dump_scalar(H5Easy::File &file, const std::string &path, const T x) {
     if constexpr (std::is_same<T, bool>::value) {
       std::vector<int> vec = {x ? 1 : 0}; // workaround for bool
       H5Easy::dump(file, path, vec);
     } else {
       std::vector<T> vec = {x};
       H5Easy::dump(file, path, vec);
     }
   }

   template<typename T>
   void h5_dump_vector(H5Easy::File &file, const std::string &path, const std::vector<T> &vec) {
     H5Easy::dump(file, path, vec);
   }

   template <real_Eigen_matrix REM>
   void h5_dump_matrix(H5Easy::File &file, const std::string &path, const REM &m) {
     H5Easy::dump(file, path, m);
   }

  template <complex_Eigen_matrix CEM>
  void h5_dump_matrix(H5Easy::File &file, const std::string &path, const CEM &m) {
    EigenMatrix<double> mr = m.real();
    h5_dump_matrix(file, path, mr);
    EigenMatrix<double> mi = m.imag();
    h5_dump_matrix(file, path + "-imag", mi);
  }
}

#endif
