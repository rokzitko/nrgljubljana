// Kramers-Kronig transformation tool
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si

// The input file must consist of one space-separated (energy, value) pair per data line. Blank lines and full-line
// comments beginning with '#' are ignored. The energy grid must be symmetric with respect to zero and contain an even
// number of points. The selected GSL interpolant is materialized as intervalwise polynomials whose principal-value
// Cauchy transform is integrated analytically.

#ifndef _kk_kk_hpp_
#define _kk_kk_hpp_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cerrno>
#include <cstddef>
#include <cstdlib>
#include <filesystem>
#include <ios>
#include <istream>
#include <ostream>
#include <sstream>
#include <vector>
#include <utility>
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

#include "../common/diagnostics.hpp"
#include "../common/gsl_piecewise_polynomial.hpp"
#include "../common/output_file.hpp"

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

inline auto parse_finite_field(const std::string &text, const std::string &context) {
  char *end = nullptr;
  errno = 0;
  const auto value = std::strtod(text.c_str(), &end);
  const auto underflowed_to_zero = errno == ERANGE && value == 0.0;
  if (underflowed_to_zero || end == text.c_str() || end != text.c_str() + text.size() || !std::isfinite(value))
    throw std::runtime_error(context + ": expected a finite representable number; got '" + text + "'.");
  return value;
}

// Read one physical record per line from stream F.
inline auto read(std::istream &F, const std::string &source = "<input>") {
  XYFUNC v;
  size_t line_number = 0;
  auto process_line = [&](const std::string &line) {
    ++line_number;
    const auto first = line.find_first_not_of(" \t\r\f\v");
    if (first == std::string::npos || line[first] == '#') return;

    std::istringstream fields(line);
    std::vector<std::string> tokens;
    std::string token;
    while (fields >> token) tokens.push_back(token);
    if (tokens.size() != 2)
      throw std::runtime_error(source + ":" + std::to_string(line_number)
                               + ": expected exactly 2 numeric fields; found " + std::to_string(tokens.size()) + ".");

    const auto context = source + ":" + std::to_string(line_number) + ": field ";
    v.push_back({parse_finite_field(tokens[0], context + "1"), parse_finite_field(tokens[1], context + "2")});
  };
  std::string line;
  try {
    while (std::getline(F, line)) {
      process_line(line);
    }
  } catch (const std::ios_base::failure &error) {
    if (F.eof() && !F.bad()) {
      if (!line.empty()) process_line(line);
    } else {
      throw std::runtime_error(source + ": I/O error after line " + std::to_string(line_number) + ": " + error.what());
    }
  }
  if (F.bad() || (F.fail() && !F.eof()))
    throw std::runtime_error(source + ": I/O error after line " + std::to_string(line_number) + ".");
  return v;
}

inline void write(const XYFUNC &re, std::ostream &F, const int prec = OUTPUT_PRECISION,
                  const std::string &source = "<output>") {
  try {
    F << std::setprecision(prec);
    for (const auto & [x,y] : re) F << x << " " << y << '\n';
    F.flush();
  } catch (const std::ios_base::failure &error) {
    throw std::runtime_error(source + ": output write or flush failed: " + error.what());
  }
  if (!F) throw std::runtime_error(source + ": output write or flush failed.");
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
   std::string input_source = "<stdin>";
   std::string output_source = "<stdout>";

   NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::akima;
   int verbosity = 0;

   void configure(const NumericalOptions &options) {
     interpolation_method = options.interpolation;
    }

    void report_configuration() const {
      if (verbosity == 0) return;
      NRG::Tools::ConfigurationReport report("kk");
      report.value("verbosity", verbosity);
      report.resolved("mode", mode == MODE::FILES ? "files" : "stdio", "positional arguments");
      report.value("input", input_source);
      report.value("output", output_source);
      report.value("interpolation", NRG::Tools::interpolation_method_name(interpolation_method));
      report.value("input.points", len);
      report.resolved("input.lower_bound", Xmin, "sorted input grid");
      report.resolved("input.upper_bound", Xmax, "sorted input grid");
      report.value("output.precision", OUTPUT_PRECISION);
      report.value("endpoint_policy", "subtracted");
      report.write(std::cerr);
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
        const auto nr = Xpts.size()/2;
       for (auto i = nr; i < len; i++)
         if (gsl_fcmp(Xpts[len - i - 1], -Xpts[i], 1e-8) != 0) {
            throw std::runtime_error("Input grid is not symmetric around zero.");
          }
        report_configuration();
        polynomial.emplace(NRG::Tools::make_gsl_piecewise_polynomial(Xpts, Ypts, interpolation_method));
        const auto sum = polynomial->integral();
        if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
        if (mode == MODE::FILES) std::cout << "Sum=" << sum << std::endl;
      }
   
    void about() {
     std::cout << "Kramers-Kronig transformation tool" << std::endl;
   }
   
   void usage() {
     std::cout << "\nUsage: kk [options] <input> <output>\n";
     std::cout << "\nAlternative usage: kk [options] -\n";
     std::cout << "\nIn this mode, kk reads from STDIN and outputs to STDOUT.\n\n";
       std::cout << "Options:\n"
                 << "  -h, --help                     show this help\n"
                 << "  -i, --interpolation METHOD     linear, cspline, akima, or steffen (default: akima)\n"
                 << "  -v                             show resolved configuration on standard error\n"
                 << "  -vv                            increase verbosity further\n"
                 << "  -V, --version                  show project version" << std::endl;
    }

    void parse_cmd_line(int argc, char *argv[]) {
      const option long_options[] = {
        {"help", no_argument, nullptr, 'h'},
        {"interpolation", required_argument, nullptr, 'i'},
        {nullptr, 0, nullptr, 0}
     };
     int option;
     while ((option = getopt_long(argc, argv, "hi:v", long_options, nullptr)) != -1) {
       switch (option) {
         case 'h':
            usage();
           std::exit(EXIT_SUCCESS);
          case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
          case 'v': ++verbosity; break;
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
        input_source = inputfn;
        output_source = outputfn;
        std::cout << inputfn << " --> " << outputfn << std::endl;
        Fin.open(inputfn);
        if (!Fin) {
          std::cerr << "Can't open " << inputfn << " for reading." << std::endl;
          std::exit(2);
        }
        std::error_code status_error;
        const auto output_status = std::filesystem::status(outputfn, status_error);
        if (status_error && status_error != std::errc::no_such_file_or_directory)
          throw std::runtime_error("Can't inspect " + outputfn + " before writing: " + status_error.message());
        if (!status_error && std::filesystem::exists(output_status)) {
          std::error_code equivalent_error;
          const auto equivalent = std::filesystem::equivalent(inputfn, outputfn, equivalent_error);
          if (equivalent_error)
            throw std::runtime_error("Can't compare input and output files: " + equivalent_error.message());
          if (equivalent)
            throw std::runtime_error("Input and output files must be different; '" + inputfn + "' and '" + outputfn
                                     + "' refer to the same filesystem object.");
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
      const auto im = read(mode == MODE::FILES ? Fin : std::cin, input_source);
      init(im);
      const auto re = calc(Xpts);
      if (mode == MODE::FILES) {
        std::ostringstream output;
        write(re, output, OUTPUT_PRECISION, output_source);
        NRG::Tools::write_output_file(output_source, output.str());
        std::cout << "KK done!" << std::endl;
        if (!std::cout) throw std::runtime_error("<stdout>: output write or flush failed.");
      } else {
        write(re, std::cout, OUTPUT_PRECISION, output_source);
      }
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
