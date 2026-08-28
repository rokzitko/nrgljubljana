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
#include <cstddef>
#include <cstdlib>
#include <ios>
#include <istream>
#include <ostream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <tuple>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <getopt.h>

#include "../common/gsl_config.hpp"

namespace NRG::KK {

using std::size_t;

struct gsl_accel_deleter {
  void operator()(gsl_interp_accel *acc) const {
    if (acc) gsl_interp_accel_free(acc);
  }
};

struct gsl_spline_deleter {
  void operator()(gsl_spline *spline) const {
    if (spline) gsl_spline_free(spline);
  }
};

struct gsl_workspace_deleter {
  void operator()(gsl_integration_workspace *workspace) const {
    if (workspace) gsl_integration_workspace_free(workspace);
  }
};

using XYPOINT = std::pair<double, double>;
using XYFUNC = std::vector<XYPOINT>;
using DVEC = std::vector<double>;

struct NumericalOptions {
  NRG::Tools::InterpolationMethod interpolation = NRG::Tools::InterpolationMethod::akima;
  double epsabs = 1e-12;
  double epsrel = 1e-8;
  size_t workspace_limit = 1000;
  NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;
};

// number of digits of precision in the generated output file
constexpr auto OUTPUT_PRECISION = 16;

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
   MODE mode = MODE::LIBRARY;
   
   int len;           // number of data points
   DVEC Xpts, Ypts;
   DVEC Xpos;         // Only positive X points [grid]
   double Xmin, Xmax; // Interval boundaries for the frequency grid
   std::unique_ptr<gsl_interp_accel, gsl_accel_deleter> acc;
   std::unique_ptr<gsl_spline, gsl_spline_deleter> spline;
   std::unique_ptr<gsl_integration_workspace, gsl_workspace_deleter> w;

   std::ifstream Fin;
   std::ofstream Fout;

   NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::akima;
   double epsabs = 1e-12;
   double epsrel = 1e-8;
   size_t workspace_limit = 1000; // workspace size for integration routine
   NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
   NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;

   void configure(const NumericalOptions &options) {
     NRG::Tools::validate_tolerances(options.epsabs, options.epsrel);
     NRG::Tools::validate_qag_workspace_limit(options.workspace_limit);
     interpolation_method = options.interpolation;
     epsabs = options.epsabs;
     epsrel = options.epsrel;
     workspace_limit = options.workspace_limit;
     quadrature_rule = options.quadrature_rule;
     gsl_error_policy = options.gsl_error_policy;
   }
   
    // Initialize the KK transformer
    void init(XYFUNC im) {  // pass by value
      if (im.empty()) throw std::runtime_error("No input data points provided.");
      NRG::Tools::validate_qag_workspace_limit(workspace_limit);
      std::sort(im.begin(), im.end());
      len = im.size();
      if (len % 2 != 0) throw std::runtime_error("Input grid must contain an even number of points.");
      const auto minimum_size = NRG::Tools::interpolation_minimum_size(interpolation_method);
      if (im.size() < minimum_size)
        throw std::runtime_error("Interpolation method " + std::string(NRG::Tools::interpolation_method_name(interpolation_method))
                                 + " requires at least " + std::to_string(minimum_size) + " input points.");
      std::tie (Xmin, Xmax) = x_range(im);
      if (mode == MODE::FILES) std::cout << "Range: [" << Xmin << " ; " << Xmax << "]" << std::endl;
      if (gsl_fcmp(-Xmin, Xmax, 1.e-8) != 0) throw std::runtime_error("Only symmetric intervals are supported!");
      tie(Xpts, Ypts) = split_vector_of_pairs(im);
      gsl_set_error_handler_off();
      acc.reset(gsl_interp_accel_alloc());
      if (!acc) throw std::runtime_error("Failed to allocate GSL interpolation accelerator.");
      // NOTE: With Akima splines there might be problems with the loss of floating point precision in the numeric
      // integration step. In cubic splines instead no such difficulties seem to appear.
      const auto Interp_type = NRG::Tools::gsl_interpolation_type(interpolation_method);
      spline.reset(gsl_spline_alloc(Interp_type, len));
      if (!spline) throw std::runtime_error("Failed to allocate GSL spline.");
      if (const auto status = gsl_spline_init(spline.get(), Xpts.data(), Ypts.data(), len); status != 0)
        throw std::runtime_error(std::string("Failed to initialize GSL spline: ") + gsl_strerror(status));
      const auto sum = gsl_spline_eval_integ(spline.get(), Xmin, Xmax, acc.get());
      if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
      if (mode == MODE::FILES) std::cout << "Sum=" << sum << std::endl;
      const auto nr = Xpts.size()/2;
      for (auto i = nr; i < len; i++)
        if (gsl_fcmp(Xpts[len - i - 1], -Xpts[i], 1e-8) != 0) {
          throw std::runtime_error("Input grid is not symmetric around zero.");
        }
      Xpos = DVEC(nr);   // Xpos are positive and increasing!
      std::copy(Xpts.begin() + nr, Xpts.end(), Xpos.begin());
      w.reset(gsl_integration_workspace_alloc(workspace_limit));
      if (!w) throw std::runtime_error("Failed to allocate GSL integration workspace.");
    }
   
   // Integrand. Method: we take a sum of the contributions for positive and negative x and return their sum. This is
   // helpful for even integrands, since it leads to possible cancellations and better accuracy of the final result, in
   // particular for small Z.
    auto f(const double X, const double Z) const {
      // [ f(x) - f(z) ] / (x-z)
      const auto a = X != Z ? (gsl_spline_eval(spline.get(), X, acc.get()) - gsl_spline_eval(spline.get(), Z, acc.get())) / (X - Z)
                            : gsl_spline_eval_deriv(spline.get(), X, acc.get());
      // [ f(-x) - f(z) ] / (-x-z)
      const auto b = -X != Z ? (gsl_spline_eval(spline.get(), -X, acc.get()) - gsl_spline_eval(spline.get(), Z, acc.get())) / (-X - Z) 
                            : gsl_spline_eval_deriv(spline.get(), -X, acc.get());
      return a + b;
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
   
    void about() {
     std::cout << "Kramers-Kronig transformation tool, RZ 2007-2020" << std::endl;
   }
   
   void usage() {
     std::cout << "\nUsage: kk [options] <input> <output>\n";
     std::cout << "\nAlternative usage: kk [options] -\n";
     std::cout << "\nIn this mode, kk reads from STDIN and outputs to STDOUT.\n\n";
     std::cout << "Options:\n"
               << "  -h, --help                     show this help\n"
               << "  -i, --interpolation METHOD     linear, cspline, or akima (default: akima)\n"
               << "      --epsabs VALUE             absolute tolerance (default: 1e-12)\n"
               << "      --epsrel VALUE             relative tolerance (default: 1e-8)\n"
               << "      --workspace-limit N        integration workspace size (default: 1000)\n"
               << "      --quadrature-rule RULE     15, 21, 31, 41, 51, or 61 (default: 15)\n"
               << "      --gsl-error-policy POLICY  ignore, warn, or fail (default: ignore)" << std::endl;
   }

   void parse_cmd_line(int argc, char *argv[]) {
     enum {
       OPT_EPSABS = 1000,
       OPT_EPSREL,
       OPT_WORKSPACE_LIMIT,
       OPT_QUADRATURE_RULE,
       OPT_GSL_ERROR_POLICY
     };
     const option long_options[] = {
       {"help", no_argument, nullptr, 'h'},
       {"interpolation", required_argument, nullptr, 'i'},
       {"epsabs", required_argument, nullptr, OPT_EPSABS},
       {"epsrel", required_argument, nullptr, OPT_EPSREL},
       {"workspace-limit", required_argument, nullptr, OPT_WORKSPACE_LIMIT},
       {"quadrature-rule", required_argument, nullptr, OPT_QUADRATURE_RULE},
       {"gsl-error-policy", required_argument, nullptr, OPT_GSL_ERROR_POLICY},
       {nullptr, 0, nullptr, 0}
     };
     int option;
     while ((option = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
       switch (option) {
         case 'h':
           usage();
           std::exit(EXIT_SUCCESS);
         case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
         case OPT_EPSABS: epsabs = NRG::Tools::parse_finite_double(optarg, "Absolute integration tolerance"); break;
         case OPT_EPSREL: epsrel = NRG::Tools::parse_finite_double(optarg, "Relative integration tolerance"); break;
         case OPT_WORKSPACE_LIMIT:
           workspace_limit = NRG::Tools::parse_positive_size(optarg, "Integration workspace limit");
           break;
         case OPT_QUADRATURE_RULE: quadrature_rule = NRG::Tools::parse_qag_rule(optarg); break;
         case OPT_GSL_ERROR_POLICY: gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg); break;
         default: throw std::invalid_argument("Unknown command-line option.");
       }
     }
     NRG::Tools::validate_tolerances(epsabs, epsrel);
     NRG::Tools::validate_qag_workspace_limit(workspace_limit);
     const auto remaining = argc - optind;
     if (remaining == 2) mode = MODE::FILES;
     if (remaining == 1 && std::strcmp(argv[optind], "-") == 0) mode = MODE::STD;
     if (mode != MODE::STD) about();
     if (mode == MODE::LIBRARY) {
       usage();
       std::exit(1);
     }
     if (mode == MODE::FILES) {
       const std::string inputfn  = argv[optind];
       const std::string outputfn = argv[optind + 1];
       std::cout << inputfn << " --> " << outputfn << std::endl;
       Fin.open(inputfn);
       if (!Fin) {
         std::cerr << "Can't open " << inputfn << " for reading." << std::endl;
         std::exit(2);
       }
       Fout.open(outputfn);
       if (!Fout) {
         std::cerr << "Can't open " << outputfn << " for writing." << std::endl;
         std::exit(2);
       }
     }
   }

 public:
   // Perform the calculation for one point. Note: this is the critical part of the code, both for computational
   // requirements as for the accuracy of the results. Optimise wisely!
   auto calc(const double Z, 
             const double EPSABS,          // numeric integration epsilon (absolute)
             const double EPSREL = 1e-8)  // numeric integration epsilon (relative)
   {
      NRG::Tools::validate_tolerances(EPSABS, EPSREL);
      auto F = Wrap([Z,this](double X) -> double { return f(X,Z); }); // wrap a C++ lambda for the C interface of GSL
      double integral = std::numeric_limits<double>::quiet_NaN();
      double integration_error = std::numeric_limits<double>::quiet_NaN();
      const int status = gsl_integration_qag(F.get(),             // integrand
                                        0,                   // lower integration boundary
                                        Xmax,                // upper integration boundary
                                        EPSABS, EPSREL,      // convergence criteria
                                        workspace_limit,     // size of workspace w
                                        NRG::Tools::gsl_qag_rule(quadrature_rule), // Gauss-Kronrod rule
                                        w.get(),             // integration workspace
                                        &integral,           // final approximation
                                        &integration_error); // estimate of absolute error
      if (NRG::Tools::gsl_integration_failed(status, integral, integration_error)
          && gsl_error_policy != NRG::Tools::GslErrorPolicy::ignore) {
        std::ostringstream message;
        message << std::setprecision(17) << "GSL QAG failed for z=" << Z << ": status=" << status << " ("
                << gsl_strerror(status) << "), result=" << integral << ", estimated_error=" << integration_error
                << ", epsabs=" << EPSABS << ", epsrel=" << EPSREL << ", workspace_limit=" << workspace_limit
                << ", quadrature_rule=" << static_cast<int>(quadrature_rule);
        if (gsl_error_policy == NRG::Tools::GslErrorPolicy::warn)
          std::cerr << "kk: warning: " << message.str() << std::endl;
        else
          throw std::runtime_error(message.str());
      }
      // Add an approximation of the (-inf,-Xmax] and [Xmax,+inf) intervals.
      const auto correction = std::abs(Z) != Xmax ? -gsl_spline_eval(spline.get(), Z, acc.get()) * 2. * gsl_atanh(Z / Xmax) : 0.0;
      const auto sum = integral + correction;
      return sum/M_PI;  // Divide by pi in the definition of the KK relation!
    }

   auto calc(const double Z) { return calc(Z, epsabs, epsrel); }
   
   // Perform the calculations for all points on a grid
   auto calc(const DVEC &grid) {
     XYFUNC result;
     result.reserve(grid.size());
     for (const auto x : grid) result.push_back({x, calc(x, epsabs, epsrel)});
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

    KK(XYFUNC im, const NumericalOptions &options) {
      configure(options);
      init(im);
    }
};

} // namespace

#endif
