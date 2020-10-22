// Numerical renormalization group
// (c) Rok Zitko, rok.zitko@ijs.si, 2005-2020

#ifndef _nrg_general_h_
#define _nrg_general_h_

#include <utility>
#include <functional>
#include <iterator>
#include <algorithm>
#include <complex>
#include <numeric>
#include <limits>
#include <memory>
#include <string>
using namespace std::string_literals;
#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>
#include <unordered_map>
#include <list>
#include <deque>
#include <set>
#include <stdexcept>

// C headers
#include <cassert>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstring>
#include <unistd.h>
#include <cstdlib> // mkdtemp

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/io/ios_state.hpp>

// ublas matrix & vector containers
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>
using namespace boost::numeric;
using namespace boost::numeric::ublas;

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>
namespace atlas = boost::numeric::bindings::atlas;

// Serialization support (used for storing to files and for MPI)
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/serialization/complex.hpp>

// This must be set manually through CXXFLAG=-DNRG_MPI
#ifdef NRG_MPI
#include <boost/optional.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
namespace mpi = boost::mpi;
#endif

#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include <fmt/color.h>
#include <fmt/ranges.h>
using namespace fmt::literals;

#include <range/v3/all.hpp>

// Support for compiler dependant optimizations

#ifdef __GNUC__
#define PUREFNC __attribute__((pure))
#define CONSTFNC __attribute__((const))
#else
#define PUREFNC
#define CONSTFNC
#endif

int myrank(); // XXX

#endif
