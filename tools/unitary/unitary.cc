// Unitary transformation tool
// Rotates matrices by performing unitary transformations
// Rok Zitko, rok.zitko@ijs.si, May 2009, June 2010

// CHANGE LOG
// 9. 6. 2010 - chop small values (parameter CHOP_TOL)
// 25.6. 2010 - support for binary input files

#define NDEBUG

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <unistd.h>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>

using namespace std;
using namespace boost::numeric;
using namespace boost::numeric::ublas;

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

void about() {
  if (!quiet) {
    cout << "# unitary -- command line unitary transformation tool" << endl;
    cout << "# Rok Zitko, rok.zitko@ijs.si, May 2009" << endl;
  }
}

void usage() { cout << "Usage: unitary <A> <B> <C>" << endl; }

void safeopen(ifstream &F, char *filename) {
  F.open(filename);
  if (!F) {
    cerr << "Error opening file " << filename << endl;
    abort();
  }
}

unsigned int countspaces(const string &s) {
  unsigned int found = 0;
  for (unsigned int i = 0; i < s.length(); i++)
    if (isspace(s[i]))
      // Only non-consecutive spaces are counted
      if (i == 0 || (i > 0 && !isspace(s[i - 1]))) found++;
  return found;
}

void get_dims(ifstream &F, unsigned int &dim1, unsigned int &dim2) {
  dim1 = dim2 = 0;
  while (F.good()) {
    string s;
    getline(F, s);
    if (!F.fail()) {
      unsigned int n = countspaces(s);
      if (dim2 == 0) {
        dim2 = n + 1; // number of columns
      } else {
        assert(dim2 == n + 1); // all rows must be equally long
      }
      dim1++; // count rows
    }
  }
}

void read_matrix_text(ifstream &F, MAT &M, unsigned int dim1, unsigned int dim2) {
  M.resize(dim1, dim2);
  for (unsigned int i = 0; i < dim1; i++) {
    for (unsigned int j = 0; j < dim2; j++) {
      double x;
      F >> x;
      assert(finite(x));
      M(i, j) = x;
    }
  }
  if (veryverbose) cout << M << endl;
  if (F.fail()) {
    cerr << "read_matrix_text() failed. Input corrupted?" << endl;
    exit(1);
  }
}

void read_matrix_bin(ifstream &F, MAT &M, unsigned int dim1, unsigned int dim2) {
  M.resize(dim1, dim2);
  for (unsigned int i = 0; i < dim1; i++) {
    for (unsigned int j = 0; j < dim2; j++) {
      double x;
      F.read((char *)&x, sizeof(double));
      if (!F) {
        cerr << "Failed reading binary data." << endl;
        exit(1);
      }
      assert(finite(x));
      M(i, j) = x;
    }
  }
  if (veryverbose) cout << M << endl;
  if (F.fail()) {
    cerr << "read_matrix_bin() failed. Input corrupted?" << endl;
    exit(1);
  }
}

void read_text(char *filename, MAT &M) {
  ifstream F;
  safeopen(F, filename);
  unsigned int dim1, dim2;
  get_dims(F, dim1, dim2);
  if (verbose) cout << filename << " [" << dim1 << " x " << dim2 << "]" << endl;
  F.close();
  safeopen(F, filename);
  read_matrix_text(F, M, dim1, dim2);
}

void read_bin(char *filename, MAT &M) {
  ifstream F;
  F.open(filename, ios_base::binary);
  if (!F) {
    cerr << "Can't open " << filename << " for reading." << endl;
    exit(1);
  }
  unsigned int dim1, dim2;
  F.read((char *)&dim1, sizeof(unsigned int));
  F.read((char *)&dim2, sizeof(unsigned int));
  if (verbose) cout << filename << " [" << dim1 << " x " << dim2 << "]" << endl;
  read_matrix_bin(F, M, dim1, dim2);
}

void read(char *filename, MAT &M, bool bin = false) {
  if (bin)
    read_bin(filename, M);
  else
    read_text(filename, M);
}

void save(char *filename, MAT &M) {
  if (verbose) cout << "Saving result to " << filename << endl;
  ofstream F(filename);
  if (!F) {
    cerr << "Can't open file for writing." << endl;
    abort();
  }
  F << setprecision(OUTPUT_PREC);
  unsigned int dim1 = M.size1();
  unsigned int dim2 = M.size2();
  for (unsigned int i = 0; i < dim1; i++) {
    for (unsigned int j = 0; j < dim2; j++) {
      double val = M(i, j);
      if (abs(val) < CHOP_TOL) val = 0.0;
      F << val << (j != dim2 - 1 ? " " : "");
    }
    F << endl;
  }
  F.close();
}

void parse_param(int argc, char *argv[]) {
  char c;
  while ((c = getopt(argc, argv, "bBqvVtls:o:c:")) != -1) {
    switch (c) {
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
      default: cerr << "Unknown argument " << c << endl; abort();
    }
  }
}

void run(int argc, char *argv[]) {
  int remaining = argc - optind; // arguments left
  if (remaining == 3) {
    char *fnA = argv[optind];
    char *fnB = argv[optind + 1];
    char *fnC = argv[optind + 2];
    MAT A, B, C;
    read(fnA, A, input_ac_bin);
    read(fnB, B);
    read(fnC, C, input_ac_bin);
    if (transpose_first) A = trans(A);
    if (transpose_last) C = trans(C);
    assert(A.size2() == B.size1());
    assert(B.size2() == C.size1());
    MAT N(B.size1(), C.size2());
    N = prod(B, C);
    MAT M(A.size1(), C.size2());
    M = prod(A, N);
    if (scale_factor != 1.0) M = scale_factor * M;
    if (veryverbose) cout << "M=" << M << endl;
    if (output_filename != nullptr) save(output_filename, M);
    return;
  }
  usage();
}

int main(int argc, char *argv[]) {
  parse_param(argc, argv);
  about();
  run(argc, argv);
}
