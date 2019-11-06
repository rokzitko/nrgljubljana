// Diagonalisation tool
// Computes eigenvalues and eigenvectors of a matrix
// Rok Zitko, rok.zitko@ijs.si, May 2009, June 2010

// CHANGE LOG
// 25.6.2010 - binary file output

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

using namespace std;

typedef std::vector<double> DVEC;

bool output_bin = false;
bool output_text = true;
bool quiet = false;
bool verbose = false;
bool veryverbose = false;
const int OUTPUT_PREC = 18;

const char *fn_val = 0;
const char *fn_vec = 0;

double scale_factor = 1.0;

void about(ostream &OUT = cout)
{
   if (!quiet) {
      OUT << "# diag -- command line eigensolver tool" << endl;
      OUT << "# Rok Zitko, rok.zitko@ijs.si, May 2009" << endl;
   }   
}

void usage(ostream &OUT = cout)
{
   OUT << "Usage: diag <input file>" << endl;
}

inline int MAX(int a, int b) { return (a > b ? a : b); }

#include "lapack.h"

#define IJ(i, j) (dim*(i)+(j))

void diagonalize(unsigned int dim, DVEC &d)
{
   double *ham = &d[0]; // contiguous storage guaranteed
   double eigenvalues[dim]; // eigenvalues on exit
   
   char jobz = 'V'; // eigenvalues and eigenvectors
   char UPLO = 'L'; // lower triangle of a is stored
   int NN = dim; // the order of the matrix
   int LDA = dim; // the leading dimension of the array a
   int INFO = 0; // 0 on successful exit
   
   int LWORK0 = -1; // length of the WORK array
   double WORK0[1];
   
   // Step 1: determine optimal LWORK
   LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, (double*)eigenvalues,
	 WORK0, &LWORK0, &INFO);
      
   assert(INFO == 0);
      
   int LWORK = int(WORK0[0]);
   if (verbose) {
      cout << "LWORK=" << LWORK << endl;
   }
   assert(LWORK > 0);
      
   const int minLWORK = MAX(1, 3*dim-1); // cf. LAPACK 3.1 dsyev.f
   if (LWORK < minLWORK) {
      cerr << "Buggy dsyev. Fixing LWORK." << endl;
      LWORK = minLWORK;
   }
      
   double WORK[LWORK];
   
   // Step 2: perform the diagonalisation
   LAPACK_dsyev(&jobz, &UPLO, &NN, ham, &LDA, (double*)eigenvalues,
	 WORK, &LWORK, &INFO);
   
   if (INFO != 0) {
      cerr << "eigensolver failed. INFO=" << INFO;
      abort();
   }
   
   // Perform the optional scaling
   for (unsigned int r = 0; r < dim; r++) {
      eigenvalues[r] *= scale_factor;
   }   

   if (verbose) {
      cout << "Eigenvalues: " << endl;
      for (unsigned int r = 0; r < dim; r++) {
	 cout << r+1 << "  " << eigenvalues[r] << endl;
      }
      cout << endl;
   }
   
   if (veryverbose) {
      cout << "Eigenvectors: " << endl;
      for (unsigned int r = 0; r < dim; r++) {
	 cout << r+1 << "  ";
	 for (unsigned int j = 0; j < dim; j++) {
	    cout << ham[ IJ(r, j) ] << (j != dim-1 ? " " : "");
	 }
	 cout << endl;
      }
      cout << endl;
   }
   
   if (fn_val != 0) {
      if (verbose) {
	 cout << "Saving eigenvalues to " << fn_val << " [text]" << endl;
      }
      
      ofstream F(fn_val);
      if (!F) {
	 cerr << "Can't open file for writing." << endl;
	 abort();
      }
      F << setprecision(OUTPUT_PREC);
      for (unsigned int r = 0; r < dim; r++) {
	 F << eigenvalues[r] << endl;
      }
      F.close();
   }

   if (fn_vec != 0 && output_text) {
      if (verbose) {
	 cout << "Saving eigenvectors to " << fn_vec << " [text]" << endl;
      }
      
      ofstream F(fn_vec);
      if (!F) {
	 cerr << "Can't open file for writing." << endl;
	 abort();
      }
      F << setprecision(OUTPUT_PREC);
      for (unsigned int r = 0; r < dim; r++) {
	 for (unsigned int j = 0; j < dim; j++) {
	    F << ham[ IJ(r, j) ] <<  (j != dim-1 ? " " : "");
	 }
	 F << endl;
      }
      F.close();
   }

   if (fn_vec != 0 && output_bin) {
      if (verbose) {
	 cout << "Saving eigenvectors to " << fn_vec << " [bin]" << endl;
      }
      
      ofstream F(fn_vec, ios_base::binary);
      if (!F) {
	 cerr << "Can't open file for writing." << endl;
	 abort();
      }

      // Save matrix dimensions
      F.write((char*)&dim, sizeof(unsigned int));
      F.write((char*)&dim, sizeof(unsigned int));
      
      for (unsigned int r = 0; r < dim; r++) {
	 for (unsigned int j = 0; j < dim; j++) {
	    const double el = ham[ IJ(r, j) ];
	    F.write((char*)&el, sizeof(double));
	 }
      }
      F.close();
   }
}

// Diagonalise a matrix obtained by reading data from stream F
void diag_stream(istream &F)
{
   DVEC data;
   
   while (F.good()) {
      double x;
      F >> x;
      if (!F.fail()) {
	 data.push_back(x);
      }
   }
   
   int size = data.size();
   int N = (int)sqrt(size);
   if (verbose) {
     cout << "size=" << size << " N=" << N << endl;
   }

   if (size % (N*N) != 0) {
      cerr << "ERROR: matrix must be square!" << endl;
   }
   
   diagonalize(N, data);
}

void parse_param(int argc, char *argv[])
{
   char c;
   
   while ((c = getopt(argc, argv, "tTbBvVqo:O:s:")) != -1) {
      switch (c) {
       case 't':
	 output_text = false;
	 break;
	 
       case 'T':
	 output_text = true;
	 break;
	 
       case 'b':
	 output_bin = false;
	 break;
	 
       case 'B':
	 output_bin = true;
	 break;
	 
       case 'v':
	 verbose = true;
	 break;
	
       case 'V':
	 veryverbose = true;
	 break;
	 
       case 'q':
	 quiet = true;
	 break;
	 
       case 'o':
	 fn_val = optarg;
	 break;
	 
       case 'O':
	 fn_vec = optarg;
	 break;
	 
       case 's':
	 scale_factor = atof(optarg);
	 break;
	 
       default:
	 cerr << "Unknown argument " << c << endl;
	 abort();
      }
   }
}

void run(int argc, char *argv[])
{
   int remaining = argc-optind; // arguments left
   
   if (remaining == 1) {
      char *filename = argv[optind];
      
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

int main(int argc, char *argv[])
{
   clock_t start_clock = clock();
   cout << setprecision(OUTPUT_PREC);
   parse_param(argc, argv);
   about();
   run(argc, argv);
   clock_t end_clock = clock();
   if (!quiet) {
      cout << "# Elapsed " << double(end_clock-start_clock)/CLOCKS_PER_SEC << " s" << endl;
   }   
}

