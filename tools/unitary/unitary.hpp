// Unitary transformation tool
// Rotates matrices by performing unitary transformations
// Rok Zitko, rok.zitko@ijs.si, 2009-2020

#ifndef _unitary_unitary_hpp_
#define _unitary_unitary_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
using namespace std::string_literals;
#include <map>

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cstdint>

#include <unistd.h>

#include <basicio.hpp>
#include <io.hpp>
#include <numerics.hpp>

namespace NRG::Unitary {

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
       default: throw std::runtime_error("Unknown argument "s + c);
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
   template<real_matrix RM>
   void run(const RM &A_, const RM &B, const RM &C_) {
     auto A = transpose_first ? NRG::trans(A_) : A_;
     auto C = transpose_last ? NRG::trans(C_) : C_;
     assert(size2(A) == size1(B));
     assert(size2(B) == size1(C));
     auto N = matrix_prod<double>(B, C); // XXX: use rotate() instead?
     auto M = matrix_prod<double>(A, N);
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
