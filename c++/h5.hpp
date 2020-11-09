#ifndef _h5_hpp_
#define _h5_hpp_

#include <iostream>
#include <string>
#include <complex>
#include <boost/numeric/ublas/matrix.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#define H5_USE_BOOST
#include <highfive/H5Easy.hpp>

//#define H5_DEBUG

namespace NRG {
   using namespace boost::numeric;

   template<typename T>
   void h5_dump_scalar(H5Easy::File &file, const std::string &path, const T x) {
     std::vector<T> vec = {x};
     H5Easy::dump(file, path, vec);
   }
   
   inline void h5_dump_matrix(H5Easy::File &file, const std::string &path, const ublas::matrix<double> &m) {
#ifdef H5_DEBUG
     fmt::print("h5_dump_matrix path={} size1={} size2={}\n", path, m.size1(), m.size2());
#endif
     H5Easy::detail::createGroupsToDataSet(file, path);
     HighFive::DataSet dataset = file.createDataSet<double>(path, HighFive::DataSpace::From(m));
     dataset.write(m);
   }

   inline void h5_dump_matrix(H5Easy::File &file, const std::string &path, const ublas::matrix<std::complex<double>> &m) {
     ublas::matrix<double> mr = ublas::real(m);
     h5_dump_matrix(file, path, mr);
   }
}

#endif
