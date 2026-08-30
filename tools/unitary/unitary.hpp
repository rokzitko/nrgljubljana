// Unitary transformation tool
// Rotates matrices by performing unitary transformations
// Rok Zitko, rok.zitko@ijs.si

#ifndef _unitary_unitary_hpp_
#define _unitary_unitary_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <utility>
#include <vector>
#include <string>
#include <optional>
#include <ostream>
#include <stdexcept>
using namespace std::string_literals;
#include <map>

#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <cstdint>

#include <unistd.h>
#include <getopt.h>

#include <basicio.hpp>
#include <io.hpp>
#include <numerics.hpp>

#include "../common/diagnostics.hpp"

namespace NRG::Unitary {

class Unitary {
 private:
   bool quiet            = false;
   int verbosity         = 0;
   bool transpose_first  = false;
   bool transpose_last   = false;
   double scale_factor   = 1.0;
   double chop_tol       = 1e-14;
   bool input_ac_bin     = false;
   std::optional<std::string> output_filename;

   void usage(std::ostream &F = std::cout) {
     F << "Usage: unitary [-h] [-b | -B] [-q] [-v|-vv] [-tl] [-s scale] [-o output_fn] [-c chop_tol] <A> <B> <C>" << std::endl;
     F << " -h, --help -- show help" << std::endl;
     F << " -v -- show resolved configuration and verbose diagnostics on stderr" << std::endl;
     F << " -vv -- include matrix dumps" << std::endl;
     F << " -V, --version -- show project version and exit" << std::endl;
     F << " -b, -B -- read A and C as text or binary matrices" << std::endl;
     F << " -q -- suppress the banner" << std::endl;
     F << " -t, -l -- transpose A or C" << std::endl;
     F << " -s scale -- scale the result" << std::endl;
     F << " -o output_fn -- write the result matrix" << std::endl;
     F << " -c chop_tol -- zero output values below this magnitude" << std::endl;
   }
   
   void parse_param(int argc, char *argv[]) {
     const option long_options[] = {
       {"help", no_argument, nullptr, 'h'},
       {nullptr, 0, nullptr, 0},
     };
     int c;
     while ((c = getopt_long(argc, argv, "hbBqvtls:o:c:", long_options, nullptr)) != -1) {
       switch (c) {
       case 'h':
         usage();
         std::exit(EXIT_SUCCESS);
       case 'b': input_ac_bin = false; break;
       case 'B': input_ac_bin = true; break;
       case 'q': quiet = true; break;
       case 'v': ++verbosity; break;
       case 't': transpose_first = true; break;
       case 'l': transpose_last = true; break;
       case 's': scale_factor = atof(optarg); break;
       case 'o': output_filename = std::string(optarg); break;
       case 'c': chop_tol = atof(optarg); break;
       default: throw std::runtime_error("Unknown argument "s + std::string(1, static_cast<char>(c)));
       }
     }
   }
   
   void about() {
     if (!quiet) {
       std::cout << "# unitary -- command line unitary transformation tool" << std::endl;
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
      if (verbosity >= 2) std::cerr << "M=" << M << std::endl;
      if (output_filename) {
        if (verbosity >= 1) std::cerr << "Saving result to " << output_filename.value() << std::endl;
        save_matrix(output_filename.value(), M, false, chop_tol);
      }
    }

    void run(const std::string &fnA, const std::string &fnB, const std::string &fnC) {
      auto A = read_matrix(fnA, input_ac_bin);
      const auto B = read_matrix(fnB, false);
      auto C = read_matrix(fnC, input_ac_bin);
      if (verbosity > 0) {
        NRG::Tools::ConfigurationReport report("unitary");
        report.value("verbosity", verbosity);
        report.value("input_a", fnA);
        report.resolved("input_a_rows", size1(A), "loaded matrix");
        report.resolved("input_a_columns", size2(A), "loaded matrix");
        report.value("input_b", fnB);
        report.resolved("input_b_rows", size1(B), "loaded matrix");
        report.resolved("input_b_columns", size2(B), "loaded matrix");
        report.value("input_c", fnC);
        report.resolved("input_c_rows", size1(C), "loaded matrix");
        report.resolved("input_c_columns", size2(C), "loaded matrix");
        report.value("input_a_c_format", input_ac_bin ? "binary" : "text");
        report.value("input_b_format", "text");
        report.value("operation", "A*B*C");
        report.value("output", output_filename.value_or("none"));
        report.value("transpose_a", transpose_first);
        report.value("transpose_c", transpose_last);
        report.value("scale", scale_factor);
        report.value("chop_tolerance", chop_tol);
        report.value("output_precision", 18);
        report.value("output_format", "text-matrix");
        report.value("quiet", quiet);
        report.write(std::cerr);
      }
      if (verbosity >= 1) std::cerr << fnA << " [" << size1(A) << " x " << size2(A) << "]" << std::endl;
      if (verbosity >= 2) {
        std::cerr << A << std::endl;
        std::cerr << B << std::endl;
      }
      if (verbosity >= 1) std::cerr << fnC << " [" << size1(C) << " x " << size2(C) << "]" << std::endl;
      if (verbosity >= 2) std::cerr << C << std::endl;
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
