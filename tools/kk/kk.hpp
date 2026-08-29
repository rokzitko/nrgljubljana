// Kramers-Kronig transformation tool
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

// The input file must consist of a table of space-separated (energy, value) pairs. The energy grid must be symmetric
// with respect to zero and the file must contain an even number of lines. The selected GSL interpolant is materialized
// as intervalwise polynomials whose principal-value Cauchy transform is integrated analytically.

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
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>

#include <gsl/gsl_math.h>

#include <unistd.h>
#include <getopt.h>

#include "../common/gsl_piecewise_polynomial.hpp"

namespace NRG::KK {

using std::size_t;

using XYPOINT = std::pair<double, double>;
using XYFUNC = std::vector<XYPOINT>;
using DVEC = std::vector<double>;

struct NumericalOptions {
  NRG::Tools::InterpolationMethod interpolation = NRG::Tools::InterpolationMethod::akima;
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

// mode==FILES: we read from a file and write to a file
// mode==STD: we read from stdin and write to stdout. No other output is sent to stdout (errors are reported to stderr).
enum class MODE { LIBRARY, FILES, STD };

class KK {
 private:
   MODE mode = MODE::LIBRARY;
   
   size_t len;        // number of data points
   DVEC Xpts, Ypts;
   double Xmin, Xmax; // Interval boundaries for the frequency grid
   std::optional<NRG::Tools::PiecewisePolynomial<double>> polynomial;

   std::ifstream Fin;
   std::ofstream Fout;

   NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::akima;

   void configure(const NumericalOptions &options) {
     interpolation_method = options.interpolation;
   }
   
     // Initialize the KK transformer
     void init(XYFUNC im) {  // pass by value
       if (im.empty()) throw std::runtime_error("No input data points provided.");
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
       polynomial.emplace(NRG::Tools::make_gsl_piecewise_polynomial(Xpts, Ypts, interpolation_method));
       const auto sum = polynomial->integral();
       if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
       if (mode == MODE::FILES) std::cout << "Sum=" << sum << std::endl;
       const auto nr = Xpts.size()/2;
      for (auto i = nr; i < len; i++)
        if (gsl_fcmp(Xpts[len - i - 1], -Xpts[i], 1e-8) != 0) {
           throw std::runtime_error("Input grid is not symmetric around zero.");
         }
     }
   
    void about() {
     std::cout << "Kramers-Kronig transformation tool, RZ 2007-2020" << std::endl;
   }
   
   void usage() {
     std::cout << "\nUsage: kk [options] <input> <output>\n";
     std::cout << "\nAlternative usage: kk [options] -\n";
     std::cout << "\nIn this mode, kk reads from STDIN and outputs to STDOUT.\n\n";
      std::cout << "Options:\n"
                << "  -h, --help                     show this help\n"
                << "  -i, --interpolation METHOD     linear, cspline, akima, or steffen (default: akima)" << std::endl;
    }

    void parse_cmd_line(int argc, char *argv[]) {
      const option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"interpolation", required_argument, nullptr, 'i'},
        {nullptr, 0, nullptr, 0}
     };
     int option;
     while ((option = getopt_long(argc, argv, "hi:", long_options, nullptr)) != -1) {
       switch (option) {
         case 'h':
            usage();
            std::exit(EXIT_SUCCESS);
          case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
          default: throw std::invalid_argument("Unknown command-line option.");
        }
      }
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
    auto calc(const double Z) const {
      if (!polynomial) throw std::logic_error("Piecewise polynomial is not initialized.");
      const auto endpoint = Z == Xmin || Z == Xmax;
      const auto policy = endpoint ? NRG::Tools::CauchyEndpointPolicy::subtracted
                                   : NRG::Tools::CauchyEndpointPolicy::reject;
      const auto canonical = NRG::Tools::cauchy_principal_value(*polynomial, Z, policy);
      return -canonical.real() / M_PI;
    }
   
   // Perform the calculations for all points on a grid
    auto calc(const DVEC &grid) const {
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

    KK(XYFUNC im, const NumericalOptions &options) {
      configure(options);
      init(im);
    }
};

} // namespace

#endif
