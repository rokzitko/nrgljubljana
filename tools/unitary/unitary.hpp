// Unitary transformation tool
// Rotates matrices by performing unitary transformations
// Rok Zitko, rok.zitko@ijs.si, May 2009, June 2010

#ifndef _unitary_unitary_hpp_
#define _unitary_unitary_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include <map>

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cstdint>

#include <unistd.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <basicio.hpp>

namespace NRG::Unitary {

using namespace std;
using namespace boost::numeric;

inline auto count_words_in_string(const std::string &s) {
  std::stringstream stream(s);
  return std::distance(std::istream_iterator<std::string>(stream), std::istream_iterator<std::string>());
}

// Determine the matrix dimensions from a stream of rows of white-space-separated tabulated values
inline auto get_dims(std::ifstream &F) {
  auto dim1 = 0; // number of rows
  auto dim2 = 0; // number of columns
  while (F.good()) {
    std::string s;
    std::getline(F, s);
    if (!F.fail()) {
      auto n = count_words_in_string(s);
      if (dim2 > 0 && dim2 != n) throw std::runtime_error("All matrix rows must be equally long");
      dim2 = n;
      dim1++;
    }
  }
  return std::make_pair(dim1, dim2);
}

using MAT = ublas::matrix<double>;

// Read dim1 x dim2 matrix from stream. Use function next_value to extract consecutive values.
template<typename FNC>
auto read_matrix_data(std::ifstream &F, FNC next_value, const size_t dim1, const size_t dim2, const bool check_is_finite = true) {
  MAT M(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto x = next_value(F);
      if (check_is_finite && !finite(x)) throw std::runtime_error("Non-finite number detected.");      
      M(i, j) = x;
    }
  }
  if (F.fail()) throw std::runtime_error("read_matrix_text() failed. Input corrupted?");
  return M;
}

// Read a matrix from stream (text)
inline auto read_matrix_text(const std::string &filename, const bool verbose) {
  auto F = safe_open_for_reading(filename, false);
  const auto [dim1, dim2] = get_dims(F);
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  F.clear();
  F.seekg (0, ios::beg);
  return read_matrix_data(F, [](auto &F) { double x; F >> x; return x; }, dim1, dim2);
}

// Read a matrix from stream (binary). Format: two unit32_t for matrix size, followed by
// dim1 x dim2 double values;
inline auto read_matrix_bin(const std::string &filename, const bool verbose) {
  auto F = safe_open_for_reading(filename, true);
  uint32_t dim1, dim2;
  F.read((char *)&dim1, sizeof(uint32_t));
  F.read((char *)&dim2, sizeof(uint32_t));
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  return read_matrix_data(F, [](auto &F) { double x; F.read((char *)&x, sizeof(double)); return x; }, dim1, dim2);
}

inline auto read_matrix(const std::string &filename, const bool bin = false, const bool verbose = false, const bool veryverbose = false) {
  auto M = bin ? read_matrix_bin(filename, verbose) : read_matrix_text(filename, verbose);
  if (veryverbose) std::cout << M << std::endl;
  return M;
}

inline void save_matrix(const std::string &filename, const MAT &M, const bool verbose = false,
                        const double chop_tol = 1e-14, const int output_prec = 18)
{
  if (verbose) std::cout << "Saving result to " << filename << std::endl;
  auto F = safe_open(filename);
  F << std::setprecision(output_prec);
  const auto dim1 = M.size1();
  const auto dim2 = M.size2();
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      const auto val = M(i, j);
      F << (std::abs(val) > chop_tol ? val : 0.0) << (j != dim2 - 1 ? " " : "");
    }
    F << std::endl;
  }
  F.close();
}

class Unitary {
 private:
   bool quiet            = false;
   bool verbose          = false;
   bool veryverbose      = false;
   bool transpose_first  = false;
   bool transpose_last   = false;
   double scale_factor   = 1.0;
   double chop_tol       = 1e-14;
   bool input_ac_bin     = false;
   std::optional<std::string> output_filename;

   void usage(std::ostream &F = std::cout) {
     F << "Usage: unitary [-h] [-b | -B] [-qvV] [-tl] [-s scale] [-o output_fn] [-c chop_tol] <A> <B> <C>" << std::endl;
   }
   
   void parse_param(int argc, char *argv[]) {
     char c;
     while (c = getopt(argc, argv, "hbBqvVtls:o:c:"), c != -1) {
       switch (c) {
       case 'h':
         usage();
         exit(EXIT_SUCCESS);
       case 'b': input_ac_bin = false; break;
       case 'B': input_ac_bin = true; break;
       case 'q': quiet = true; break;
       case 'v': verbose = true; break;
       case 'V': veryverbose = true; break;
       case 't': transpose_first = true; break;
       case 'l': transpose_last = true; break;
       case 's': scale_factor = atof(optarg); break;
       case 'o': output_filename = std::string(optarg); break;
       case 'c': chop_tol = atof(optarg); break;
       default: throw std::runtime_error("Unknown argument " + c);
       }
     }
   }
   
   void about() {
     if (!quiet) {
       std::cout << "# unitary -- command line unitary transformation tool" << std::endl;
       std::cout << "# Rok Zitko, rok.zitko@ijs.si, 2009-2020" << std::endl;
     }
   }

 public:
   void run(MAT &A, const MAT &B, MAT &C) { // A and B may be transposed
     if (transpose_first) A = ublas::trans(A);
     if (transpose_last) C = ublas::trans(C);
     assert(A.size2() == B.size1());
     assert(B.size2() == C.size1());
     MAT N = ublas::prod(B, C);
     MAT M = ublas::prod(A, N);
     if (scale_factor != 1.0) M = scale_factor * M;
     if (veryverbose) std::cout << "M=" << M << std::endl;
     if (output_filename) save_matrix(output_filename.value(), M, verbose, chop_tol);
   }

   void run(const std::string &fnA, const std::string &fnB, const std::string &fnC) {
     auto A = read_matrix(fnA, input_ac_bin, verbose, veryverbose);
     const auto B = read_matrix(fnB, false);
     auto C = read_matrix(fnC, input_ac_bin, verbose, veryverbose);
     run(A, B, C);
   }

   void run(int argc, char *argv[]) {
     const auto remaining = argc - optind; // arguments left after switch parsing
     if (remaining == 3) {
       const std::string fnA = argv[optind];
       const std::string fnB = argv[optind+1];
       const std::string fnC = argv[optind+2];
       run(fnA, fnB, fnC);
     } else {
       usage();
     }
   }
   
   Unitary(int argc, char *argv[]) {
     parse_param(argc, argv);
     about();
   }
};
   
} // namespace

#endif
