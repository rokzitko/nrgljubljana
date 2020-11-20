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

bool quiet            = false;
bool verbose          = false;
bool veryverbose      = false;
bool transpose_first  = false;
bool transpose_last   = false;
double scale_factor   = 1.0;
char *output_filename = nullptr;

bool input_ac_bin = false;

const int OUTPUT_PREC = 18;
double CHOP_TOL       = 1e-14;

using MAT = ublas::matrix<double>;

auto count_breaks(const std::string &s) {
  auto found = 0;
  for (auto i = 0; i < s.length(); i++)
    if (isspace(s[i]))
      // Only non-consecutive spaces are counted
      if (i == 0 || (i > 0 && !isspace(s[i-1]))) found++;
  return found;
}

auto get_dims(std::ifstream &F) {
  auto dim1 = 0; // number of rows
  auto dim2 = 0; // number of columns
  while (F.good()) {
    std::string s;
    std::getline(F, s);
    if (!F.fail()) {
      auto n = count_breaks(s);
      if (dim2>0 && dim2 != n+1) throw std::runtime_error("All matrix rows must be equally long");
      dim2 = n+1;
      dim1++;
    }
  }
  return std::make_pair(dim1, dim2);
}

auto read_matrix_text_data(std::ifstream &F, const size_t dim1, const size_t dim2) {
  MAT M(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      double x;
      F >> x;
      assert(finite(x));
      M(i, j) = x;
    }
  }
  if (veryverbose) std::cout << M << std::endl;
  if (F.fail()) throw std::runtime_error("read_matrix_text() failed. Input corrupted?");
  return M;
}

auto read_matrix_bin_data(std::ifstream &F, const size_t dim1, const size_t dim2) {
  MAT M(dim1, dim2);
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      double x;
      F.read((char *)&x, sizeof(double));
      if (!F) throw std::runtime_error("Failed reading binary data.");
      assert(finite(x));
      M(i, j) = x;
    }
  }
  if (veryverbose) std::cout << M << std::endl;
  if (F.fail()) throw std::runtime_error("read_matrix_bin() failed. Input corrupted?");
  return M;
}

auto read_text_matrix(const std::string &filename) {
  auto F = safe_open_for_reading(filename);
  const auto [dim1, dim2] = get_dims(F);
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  F.clear();
  F.seekg (0, ios::beg);
  return read_matrix_text_data(F, dim1, dim2);
}

auto read_bin_matrix(const std::string &filename) {
  auto F = safe_open_for_reading(filename, true);
  uint32_t dim1, dim2;
  F.read((char *)&dim1, sizeof(uint32_t));
  F.read((char *)&dim2, sizeof(uint32_t));
  if (verbose) std::cout << filename << " [" << dim1 << " x " << dim2 << "]" << std::endl;
  return read_matrix_bin_data(F, dim1, dim2);
}

auto read_matrix(const std::string &filename, const bool bin = false) {
  return bin ? read_bin_matrix(filename) : read_text_matrix(filename);
}

void save_matrix(const std::string &filename, const MAT &M) {
  if (verbose) std::cout << "Saving result to " << filename << std::endl;
  auto F = safe_open(filename);
  F << std::setprecision(OUTPUT_PREC);
  const auto dim1 = M.size1();
  const auto dim2 = M.size2();
  for (auto i = 0; i < dim1; i++) {
    for (auto j = 0; j < dim2; j++) {
      double val = M(i, j);
      if (abs(val) < CHOP_TOL) val = 0.0;
      F << val << (j != dim2 - 1 ? " " : "");
    }
    F << std::endl;
  }
  F.close();
}

void usage(std::ostream &F = std::cout)
{
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
      case 'o': output_filename = optarg; break;
      case 'c': CHOP_TOL = atof(optarg); break;
      default: throw std::runtime_error("Unknown argument " + c);
    }
  }
}

void run(MAT &A, const MAT &B, MAT &C) // A and B may be transposed
{
  if (transpose_first) A = ublas::trans(A);
  if (transpose_last) C = ublas::trans(C);
  assert(A.size2() == B.size1());
  assert(B.size2() == C.size1());
  //    MAT N(B.size1(), C.size2());
  MAT N = ublas::prod(B, C);
  //    MAT M(A.size1(), C.size2());
  MAT M = ublas::prod(A, N);
  if (scale_factor != 1.0) M = scale_factor * M;
  if (veryverbose) std::cout << "M=" << M << std::endl;
  if (output_filename != nullptr) save_matrix(output_filename, M);
}

void run(const std::string &fnA, const std::string &fnB, const std::string &fnC) {
  auto A = read_matrix(fnA, input_ac_bin);
  const auto B = read_matrix(fnB);
  auto C = read_matrix(fnC, input_ac_bin);
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

void about() {
  if (!quiet) {
    std::cout << "# unitary -- command line unitary transformation tool" << std::endl;
    std::cout << "# Rok Zitko, rok.zitko@ijs.si, 2009-2020" << std::endl;
  }
}

} // namespace

#endif
