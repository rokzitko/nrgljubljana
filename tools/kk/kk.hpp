// Kramers-Kronig transformation tool
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

// The input file must consist of a table of space-separated (energy, value) pairs. The energy grid must be symmetric
// with respect to zero and the file must contain an even number of lines. Gauss-Kronrod quadrature rules are used.
// At singularity points the derivative is computed using GSL interpolation routines (cubic splines).

// NOTE about Gauss-Kronrod: The higher-order rules give better accuracy for smooth functions, while lower-order
// rules save time when the function contains local difficulties, such as discontinuities. [GSL manual] On each
// iteration the adaptive integration strategy bisects the interval with the largest error estimate. The subintervals
// and their results are stored in the memory provided by workspace. The maximum number of subintervals is given by
// limit, which may not exceed the allocated size of the workspace. [GSL manual]

#ifndef _kk_kk_hpp_
#define _kk_kk_hpp_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>

namespace NRG::KK {

using XYPOINT = std::pair<double, double>;
using XYFUNC = std::vector<XYPOINT>;
using DVEC = std::vector<double>;

// number of digits of precision in the generated output file
constexpr auto OUTPUT_PRECISION = 16;

// Output width for verbose mode
constexpr auto WIDTH = 24;

// Read data from stream F.
inline auto read(std::istream &F) {
  XYFUNC v;
  while (F) {
    if (F.peek() == '#') { // skip comment lines
      std::string line;
      std::getline(F, line);
    } else {
      double x, y;
      F >> x >> y;
      if (F.fail()) break;
      assert(std::isfinite(x) && std::isfinite(y));
      v.push_back({x, y});
    }
  }
  return v;
}

inline void write(const XYFUNC &re, std::ostream &F, const int prec = OUTPUT_PRECISION) {
  F << std::setprecision(prec);
  for (const auto & [x,y] : re) F << x << " " << y << std::endl;
}

inline auto x_range(const XYFUNC &l) 
{
  return std::make_pair(l.front().first, l.back().first);
}

// Transform a vector of pairs into a pair of vectors
template<typename S, typename T>
auto split_vector_of_pairs(const std::vector<std::pair<S,T>> &v)
{
  std::vector<S> x;
  x.reserve(v.size());
  std::vector<T> y;
  y.reserve(v.size());
  for (const auto &[i,j]: v) {
    x.push_back(i);
    y.push_back(j);
  }
  return std::make_pair(x,y);
}

inline double unwrap(const double x, void *p) {
  auto fp = static_cast<std::function<double(double)> *>(p);
  return (*fp)(x);
}

// mode==FILES: we read from a file and write to a file
// mode==STD: we read from stdin and write to stdout. No other output is sent to stdout (errors are reported to stderr).
enum class MODE { LIBRARY, FILES, STD };

class KK {
 private:
   bool verbose     = false;
   bool veryverbose = false;

   MODE mode = MODE::LIBRARY;
   
   int len;           // number of data points
   DVEC Xpts, Ypts;
   DVEC Xpos;         // Only positive X points [grid]
   double Xmin, Xmax; // Interval boundaries for the frequency grid
   gsl_interp_accel *acc;
   gsl_spline *spline;
   gsl_integration_workspace *w;

   std::ifstream Fin;
   std::ofstream Fout;

   inline static const size_t workspace_limit = 1000; // workspace size for integration routine
   
   // Initialize the KK transformer
   void init(XYFUNC im) {  // pass by value
     std::sort(im.begin(), im.end());
     len = im.size();
     assert(len % 2 == 0);
     std::tie (Xmin, Xmax) = x_range(im);
     if (mode == MODE::FILES) std::cout << "Range: [" << Xmin << " ; " << Xmax << "]" << std::endl;
     if (gsl_fcmp(-Xmin, Xmax, 1.e-8) != 0) throw std::runtime_error("Only symmetric intervals are supported!");
     tie(Xpts, Ypts) = split_vector_of_pairs(im);
     acc = gsl_interp_accel_alloc();
     // NOTE: With akime splines the might be problems with the loss of the floating point precision in the numeric
     // integration step. In cubic splines instead no such difficulties seem to appear.
     // gsl_interp_linear;
     // gsl_interp_cspline;
     const auto Interp_type = gsl_interp_akima;
     spline                 = gsl_spline_alloc(Interp_type, len);
     gsl_spline_init(spline, Xpts.data(), Ypts.data(), len);
     const auto sum = gsl_spline_eval_integ(spline, Xmin, Xmax, acc);
     if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
     if (mode == MODE::FILES) std::cout << "Sum=" << sum << std::endl;
     const auto nr = Xpts.size()/2;
     for (auto i = nr; i < len; i++)
       assert(gsl_fcmp(Xpts[len - i - 1], -Xpts[i], 1e-8) == 0);     // Check for the symmetry of the grid
     Xpos = DVEC(nr);   // Xpos are positive and increasing!
     std::copy(Xpts.begin() + nr, Xpts.end(), Xpos.begin());
     w = gsl_integration_workspace_alloc(workspace_limit);
     gsl_set_error_handler_off();
   }
   
   // Integrand. Method: we take a sum of the contributions for positive and negative x and return their sum. This is
   // helpful for even integrands, since it leads to possible cancellations and better accuracy of the final result, in
   // particular for small Z.
   auto f(const double X, const double Z) const {
     // [ f(x) - f(z) ] / (x-z)
     const auto a = X != Z ? (gsl_spline_eval(spline, X, acc) - gsl_spline_eval(spline, Z, acc)) / (X - Z)
                           : gsl_spline_eval_deriv(spline, X, acc);
     // [ f(-x) - f(z) ] / (-x-z)
     const auto b = -X != Z ? (gsl_spline_eval(spline, -X, acc) - gsl_spline_eval(spline, Z, acc)) / (-X - Z) 
                            : gsl_spline_eval_deriv(spline, -X, acc);
     return a + b;
   }

   void handle_qag(const int status) {
     if (status && veryverbose) std::cerr << "WARNING - qag error: " << status << " -- " << gsl_strerror(status) << std::endl;
   }

   class Wrap {
    private:
      gsl_function F;
      std::function<double(double)> fnc;
    public:
      Wrap(std::function<double(double)> fnc_) : fnc(fnc_) {
        F.function = &unwrap;
        F.params   = &fnc;
      }
      auto get() { return &F; }
   };
   
   void gsl_done() {
     gsl_spline_free(spline);
     gsl_interp_accel_free(acc);
     gsl_integration_workspace_free(w);
   }

   void about() {
     std::cout << "Kramers-Kronig transformation tool, RZ 2007-2020" << std::endl;
   }
   
   void usage() {
     std::cout << "\nUsage: kk [-h] <input> <output>\n";
     std::cout << "\nAlternative usage: kk -\n";
     std::cout << "\nIn this mode, kk reads from STDIN and outputs to STDOUT." << std::endl;
   }

   void parse_cmd_line(int argc, char *argv[]) {
     if (argc == 2 && strcmp(argv[1], "-h") == 0) {
       usage();
       exit(EXIT_SUCCESS);
     }
     if (argc == 3) mode = MODE::FILES;
     if (argc == 2 && strcmp(argv[1], "-") == 0) mode = MODE::STD;
     if (mode != MODE::STD) about();
     if (mode == MODE::LIBRARY) {
       usage();
       exit(1);
     }
     if (mode == MODE::FILES) {
       const std::string inputfn  = argv[1];
       const std::string outputfn = argv[2];
       std::cout << inputfn << " --> " << outputfn << std::endl;
       Fin.open(inputfn);
       if (!Fin) {
         std::cerr << "Can't open " << inputfn << " for reading." << std::endl;
         exit(2);
       }
       Fout.open(outputfn);
       if (!Fout) {
         std::cerr << "Can't open " << outputfn << " for writing." << std::endl;
         exit(2);
       }
     }
   }

 public:
   // Perform the calculation for one point. Note: this is the critical part of the code, both for computational
   // requirements as for the accuracy of the results. Optimise wisely!
   auto calc(const double Z, 
             const double EPSABS = 1e-12, // numeric integration epsilon (absolute)
             const double EPSREL = 1e-8)  // numeric integration epsilon (relative)
   {
     auto F = Wrap([Z,this](double X) -> double { return f(X,Z); }); // wrap a C++ lambda for the C interface of GSL
     double integral;
     double integration_error;
     int status = gsl_integration_qag(F.get(),             // integrand
                                      0,                   // lower integration boundary
                                      Xmax,                // upper integration boundary
                                      EPSABS, EPSREL,      // convergence criteria
                                      workspace_limit,     // size of workspace w
                                      GSL_INTEG_GAUSS15,   // Gauss-Kronrod rule
                                      w,                   // integration workspace
                                      &integral,           // final approximation
                                      &integration_error); // estimate of absolute error
     handle_qag(status);
     // Add an approximation of the (-inf,-Xmax] and [Xmax,+inf) intervals.
     const auto correction = std::abs(Z) != Xmax ? -gsl_spline_eval(spline, Z, acc) * 2. * gsl_atanh(Z / Xmax) : 0.0;
     const auto sum = integral + correction;
     if (mode != MODE::STD && verbose) {
       std::cout << std::scientific;
       std::cout << std::setw(WIDTH) << Z << " t=" << std::setw(WIDTH) << integral / M_PI << " c=" << std::setw(WIDTH) << correction / M_PI << " ratio=" << std::setw(WIDTH)
         << correction / integral << std::endl;
       std::cout.unsetf(std::ios_base::scientific);
     }
     return sum/M_PI;  // Divide by pi in the definition of the KK relation!
   }
   
   // Perform the calculations for all points on a grid
   auto calc(const DVEC &grid) {
     XYFUNC result;
     result.reserve(grid.size());
     for (const auto x : grid) result.push_back({x, calc(x)});
     return result;
   }

   // Legacy interface when kk is used as a command-line tool
   KK(int argv, char *argc[]) {
     parse_cmd_line(argv, argc);
     const auto im = read(mode == MODE::FILES ? Fin : std::cin);
     init(im);
     const auto re = calc(Xpts);
     write(re, mode == MODE::FILES ? Fout : std::cout);
     std::cout << "KK done!" << std::endl;
   }
   
   // Modern interface when kk is used as a library
   KK(XYFUNC im) {
     init(im);
   }
   
   ~KK() { gsl_done(); }
};

} // namespace

#endif
