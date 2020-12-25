#ifndef _h5_hpp_
#define _h5_hpp_

#include <iostream>
#include <string>
#include <complex>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include "traits.hpp"

#ifdef INCL_UBLAS
#define H5_USE_BOOST // for H5Easy
#endif
#ifdef INCL_EIGEN
#define H5_USE_EIGEN // for H5Easy
#endif
#include <highfive/H5Easy.hpp>

namespace NRG {
   using namespace boost::numeric;

   template<typename T>
   void h5_dump_scalar(H5Easy::File &file, const std::string &path, const T x) {
     std::vector<T> vec = {x};
     H5Easy::dump(file, path, vec);
   }

   template<typename T>
   void h5_dump_vector(H5Easy::File &file, const std::string &path, const std::vector<T> &vec) {
     H5Easy::dump(file, path, vec);
   }

#ifdef INCL_UBLAS
   template<real_ublas_matrix RUM>
   void h5_dump_matrix(H5Easy::File &file, const std::string &path, const RUM &m) {
     H5Easy::detail::createGroupsToDataSet(file, path);
     HighFive::DataSet dataset = file.createDataSet<double>(path, HighFive::DataSpace::From(m));
     dataset.write(m);
   }

   template<complex_ublas_matrix CUM>
   void h5_dump_matrix(H5Easy::File &file, const std::string &path, const CUM &m) {
     ublas::matrix<double> mr = ublas::real(m);
     h5_dump_matrix(file, path, mr);
   }
#endif

#ifdef INCL_EIGEN
   template <real_Eigen_matrix REM>
   void h5_dump_matrix(H5Easy::File &file, const std::string &path, const REM &m) {
      H5Easy::detail::createGroupsToDataSet(file, path);
      HighFive::DataSet dataset = file.createDataSet<double>(path, HighFive::DataSpace::From(m));
      dataset.write(m);
    }

  template <complex_Eigen_matrix CEM>
  void h5_dump_matrix(H5Easy::File &file, const std::string &path, const CEM &m) {
     EigenMatrix<double> mr = m.real();
     h5_dump_matrix(file, path, mr);
   }
#endif
}

#endif
