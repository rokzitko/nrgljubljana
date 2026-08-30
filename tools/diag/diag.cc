// Diagonalisation tool
// Computes eigenvalues and eigenvectors of a matrix

#include <cstddef>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <ios>
#include <istream>
#include <ostream>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <utility>
#include <algorithm>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <limits>
#include <unistd.h>

#include <common/version.hpp>

#include "../common/diagnostics.hpp"

using namespace std;

using DVEC = std::vector<double>;

bool output_bin       = false;
bool output_text      = true;
bool quiet            = false;
bool verbose          = false;
bool veryverbose      = false;
const int OUTPUT_PREC = 18;

const char *fn_val = nullptr;
const char *fn_vec = nullptr;
const char *fn_input = nullptr;

double scale_factor = 1.0;

void about(ostream &OUT = cout) {
  if (!quiet) {
    OUT << "# diag -- command line eigensolver tool" << endl;
  }
}

#include "linalg.hpp"

#define IJ(i, j) (dim * (i) + (j))

void report_configuration(const unsigned int dim,
                          const std::size_t input_elements,
                          const lapack_int lwork,
                          const bool lwork_raised) {
  NRG::Tools::ConfigurationReport report("diag");
  report.value("verbosity", veryverbose ? 2 : 1);
  report.value("quiet", quiet);
  report.value("input_file", fn_input);
  report.value("input_mode", "text_matrix");
  report.value("input_elements", input_elements);
  report.resolved("dimension", dim, "square root of input element count");
  report.value("eigensolver", "LAPACK_dsyev");
  report.value("eigenvectors", true);
  report.value("stored_triangle", "lower");
  report.value("scale_factor", scale_factor);
  report.resolved("LWORK", lwork,
                  lwork_raised ? "LAPACK workspace query raised to documented minimum" : "LAPACK workspace query");
  report.value("eigenvalue_file", fn_val != nullptr ? fn_val : "disabled");
  report.value("eigenvector_file", fn_vec != nullptr ? fn_vec : "disabled");
  report.value("eigenvector_text", output_text);
  report.value("eigenvector_binary", output_bin);
  report.value("output_precision", OUTPUT_PREC);
  report.write(cerr);
}

void diagonalize(unsigned int dim, DVEC &d) {
  if (dim > static_cast<unsigned int>(std::numeric_limits<lapack_int>::max())) {
    cerr << "matrix dimension exceeds LAPACK integer range" << endl;
    abort();
  }

  double *ham = &d[0];     // contiguous storage guaranteed
  std::vector<double> eigenvalues(dim); // eigenvalues on exit

  char jobz = 'V'; // eigenvalues and eigenvectors
  char UPLO = 'L'; // lower triangle of a is stored
  lapack_int NN   = dim; // the order of the matrix
  lapack_int LDA  = dim; // the leading dimension of the array a
  lapack_int INFO = 0;   // 0 on successful exit

  lapack_int LWORK0 = -1; // length of the WORK array
  std::vector<double> WORK0(1);

  // Step 1: determine optimal LWORK
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK0.data(), &LWORK0, &INFO);

  assert(INFO == 0);

  lapack_int LWORK = lapack_int(WORK0[0]);
  assert(LWORK > 0);

  const lapack_int minLWORK = std::max(lapack_int(1), 3 * NN - 1); // cf. LAPACK 3.1 dsyev.f
  bool lwork_raised = false;
  if (LWORK < minLWORK) {
    cerr << "Buggy dsyev. Fixing LWORK." << endl;
    LWORK = minLWORK;
    lwork_raised = true;
  }
  if (verbose) { cerr << "LWORK=" << LWORK << endl; }
  if (verbose) report_configuration(dim, d.size(), LWORK, lwork_raised);

  std::vector<double> WORK(static_cast<std::size_t>(LWORK));

  // Step 2: perform the diagonalisation
  LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, eigenvalues.data(), WORK.data(), &LWORK, &INFO);

  if (INFO != 0) {
    cerr << "eigensolver failed. INFO=" << INFO;
    abort();
  }

  // Perform the optional scaling
  for (unsigned int r = 0; r < dim; r++) { eigenvalues[r] *= scale_factor; }

  if (verbose) {
    cerr << "Eigenvalues: " << endl;
    for (unsigned int r = 0; r < dim; r++) { cerr << r + 1 << "  " << eigenvalues[r] << endl; }
    cerr << endl;
  }

  if (veryverbose) {
    cerr << "Eigenvectors: " << endl;
    for (unsigned int r = 0; r < dim; r++) {
      cerr << r + 1 << "  ";
      for (unsigned int j = 0; j < dim; j++) { cerr << ham[IJ(r, j)] << (j != dim - 1 ? " " : ""); }
      cerr << endl;
    }
    cerr << endl;
  }

  if (fn_val != nullptr) {
    if (verbose) { cerr << "Saving eigenvalues to " << fn_val << " [text]" << endl; }

    ofstream F(fn_val);
    if (!F) {
      cerr << "Can't open file for writing." << endl;
      abort();
    }
    F << setprecision(OUTPUT_PREC);
    for (unsigned int r = 0; r < dim; r++) { F << eigenvalues[r] << endl; }
    F.close();
  }

  if (fn_vec != nullptr && output_text) {
    if (verbose) { cerr << "Saving eigenvectors to " << fn_vec << " [text]" << endl; }

    ofstream F(fn_vec);
    if (!F) {
      cerr << "Can't open file for writing." << endl;
      abort();
    }
    F << setprecision(OUTPUT_PREC);
    for (unsigned int r = 0; r < dim; r++) {
      for (unsigned int j = 0; j < dim; j++) { F << ham[IJ(r, j)] << (j != dim - 1 ? " " : ""); }
      F << endl;
    }
    F.close();
  }

  if (fn_vec != nullptr && output_bin) {
    if (verbose) { cerr << "Saving eigenvectors to " << fn_vec << " [bin]" << endl; }

    ofstream F(fn_vec, ios_base::binary);
    if (!F) {
      cerr << "Can't open file for writing." << endl;
      abort();
    }

    // Save matrix dimensions
    F.write((char *)&dim, sizeof(unsigned int));
    F.write((char *)&dim, sizeof(unsigned int));

    for (unsigned int r = 0; r < dim; r++) {
      for (unsigned int j = 0; j < dim; j++) {
        const double el = ham[IJ(r, j)];
        F.write((char *)&el, sizeof(double));
      }
    }
    F.close();
  }
}

// Diagonalise a matrix obtained by reading data from stream F
void diag_stream(istream &F) {
  DVEC data;

  while (F.good()) {
    double x;
    F >> x;
    if (!F.fail()) { data.push_back(x); }
  }

  const auto size = data.size();
  if (size == 0) {
    cerr << "ERROR: matrix input is empty!" << endl;
    exit(1);
  }
  const auto N = static_cast<int>(sqrt(size));
  if (verbose) { cerr << "size=" << size << " N=" << N << endl; }

  if (N == 0 || size != static_cast<size_t>(N) * static_cast<size_t>(N)) {
    cerr << "ERROR: matrix must be square!" << endl;
    exit(1);
  }

  diagonalize(static_cast<unsigned int>(N), data);
}

void usage(ostream &OUT = cout) { 
  OUT << "Usage: diag [-h] [-t | -T] [-b | -B] [-v|-vv] [-q] [-o fn_val] [-O fn_vec] [-s scale] <input file>" << endl;
  OUT << "  -v              print resolved configuration and enable verbose diagnostics" << endl;
  OUT << "  -vv             enable very verbose diagnostics" << endl;
  OUT << "  -V, --version   print version and exit" << endl;
}

void parse_param(int argc, char *argv[]) {
  int c;
  while ((c = getopt(argc, argv, "htTbBvqo:O:s:")) != -1) {
    switch (c) {
      case 'h':
        usage();
        exit(EXIT_SUCCESS);

      case 't': output_text = false; break;

      case 'T': output_text = true; break;

      case 'b': output_bin = false; break;

      case 'B': output_bin = true; break;

      case 'v':
        if (verbose) veryverbose = true;
        verbose = true;
        break;

      case 'q': quiet = true; break;

      case 'o': fn_val = optarg; break;

      case 'O': fn_vec = optarg; break;

      case 's': scale_factor = atof(optarg); break;

      default: cerr << "Unknown argument " << c << endl; abort();
    }
  }
}

void run(int argc, char *argv[]) {
  int remaining = argc - optind; // arguments left

  if (remaining == 1) {
    char *filename = argv[optind];
    fn_input = filename;

    ifstream F(filename);
    if (!F) {
      cerr << "Error opening file " << filename << endl;
      abort();
    }
    diag_stream(F);
    F.close();
    return;
  }

  usage();
}

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "diag")) return EXIT_SUCCESS;
  clock_t start_clock = clock();
  cout << setprecision(OUTPUT_PREC);
  cerr << setprecision(OUTPUT_PREC);
  parse_param(argc, argv);
  about();
  run(argc, argv);
  clock_t end_clock = clock();
  if (!quiet) { cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << endl; }
}
