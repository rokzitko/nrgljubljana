// Kramers-Kronig transformation tool
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si

// The input file must consist of one space-separated (energy, value) pair per data line. Blank lines and full-line
// comments beginning with '#' are ignored. The energy grid must be symmetric with respect to zero and contain an even
// number of points. By default, the principal-value transform is evaluated with adaptive QAG integration. An analytic
// interval-polynomial transform with guarded wider and exact-arithmetic fallbacks is also selectable.

#ifndef _kk_kk_hpp_
#define _kk_kk_hpp_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <exception>
#include <filesystem>
#include <ios>
#include <istream>
#include <memory>
#include <mutex>
#include <ostream>
#include <sstream>
#include <string_view>
#include <thread>
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

#include <unistd.h>
#include <getopt.h>

#include "../common/diagnostics.hpp"
#include "../common/gsl_piecewise_polynomial.hpp"
#include "../common/output_file.hpp"
#include "../common/parallel.hpp"
#include "../common/tabulated.hpp"

namespace NRG::KK {

using std::size_t;

using XYPOINT = std::pair<double, double>;
using XYFUNC = std::vector<XYPOINT>;
using DVEC = std::vector<double>;

enum class Algorithm { analytic, qag };

inline auto algorithm_name(const Algorithm algorithm) -> std::string_view {
  switch (algorithm) {
    case Algorithm::analytic: return "analytic";
    case Algorithm::qag: return "qag";
  }
  throw std::logic_error("Unknown KK algorithm.");
}

inline auto parse_algorithm(const std::string_view value) {
  if (value == "analytic") return Algorithm::analytic;
  if (value == "qag") return Algorithm::qag;
  throw std::invalid_argument("Algorithm must be one of: analytic, qag.");
}

inline auto parse_worker_count(const std::string_view value, const std::string_view source) {
  return NRG::Tools::parse_worker_count(value, source);
}

struct NumericalOptions {
  NRG::Tools::InterpolationMethod interpolation = NRG::Tools::InterpolationMethod::akima;
  Algorithm algorithm = Algorithm::qag;
  double epsabs = 1e-12;
  double epsrel = 1e-8;
  size_t workspace_limit = 1000;
  NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
  NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;
  std::optional<size_t> jobs;
};

struct gsl_accel_deleter {
  void operator()(gsl_interp_accel *accelerator) const {
    if (accelerator) gsl_interp_accel_free(accelerator);
  }
};

struct gsl_workspace_deleter {
  void operator()(gsl_integration_workspace *workspace) const {
    if (workspace) gsl_integration_workspace_free(workspace);
  }
};

// number of digits of precision in the generated output file
constexpr auto OUTPUT_PRECISION = 16;

// Read one physical record per line from stream F.
inline auto read(std::istream &F, const std::string &source = "<input>") {
  return NRG::Tools::read_strict_pairs(F, source);
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
   std::unique_ptr<gsl_spline, NRG::Tools::GslSplineDeleter> qag_spline;

   std::ifstream Fin;
   std::string input_source = "<stdin>";
   std::string output_source = "<stdout>";

   NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::akima;
   Algorithm algorithm = Algorithm::qag;
   double epsabs = 1e-12;
   double epsrel = 1e-8;
   size_t workspace_limit = 1000;
   NRG::Tools::QagRule quadrature_rule = NRG::Tools::QagRule::gauss15;
   NRG::Tools::GslErrorPolicy gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;
   size_t jobs = 1;
   std::string jobs_source = "default";
   bool jobs_resolved = false;
   int verbosity = 0;

   struct AnalyticWorker {};

   struct QagWorker {
     std::unique_ptr<gsl_interp_accel, gsl_accel_deleter> accelerator;
     std::unique_ptr<gsl_integration_workspace, gsl_workspace_deleter> workspace;

     explicit QagWorker(const size_t limit)
         : accelerator{gsl_interp_accel_alloc()}, workspace{gsl_integration_workspace_alloc(limit)} {
       if (!accelerator) throw std::runtime_error("Failed to allocate GSL interpolation accelerator.");
       if (!workspace) throw std::runtime_error("Failed to allocate GSL integration workspace.");
     }
   };

    struct QagContext {
      const gsl_spline *spline = nullptr;
      gsl_interp_accel *accelerator = nullptr;
      double argument = 0.0;
      double value_at_argument = 0.0;
      bool paired_symmetric_domain = true;
      int interpolation_status = GSL_SUCCESS;
    };

   static auto qag_integrand(const double x, void *raw_context) noexcept {
     auto &context = *static_cast<QagContext *>(raw_context);
     auto evaluate = [&context](const double point) {
       double value = std::numeric_limits<double>::quiet_NaN();
       if (context.interpolation_status == GSL_SUCCESS)
         context.interpolation_status = gsl_spline_eval_e(context.spline, point, context.accelerator, &value);
       return value;
     };
     auto derivative = [&context](const double point) {
       double value = std::numeric_limits<double>::quiet_NaN();
       if (context.interpolation_status == GSL_SUCCESS)
         context.interpolation_status = gsl_spline_eval_deriv_e(context.spline, point, context.accelerator, &value);
       return value;
     };

      const auto regularized = [&](const double point) {
        return point != context.argument
                 ? (evaluate(point) - context.value_at_argument) / (point - context.argument)
                 : derivative(point);
      };
      const auto positive = regularized(x);
      if (!context.paired_symmetric_domain) return positive;
      const auto negative_x = -x;
      return positive + regularized(negative_x);
   }

   void resolve_jobs(const std::optional<size_t> requested, const std::string_view requested_source) {
     if (requested) {
       if (*requested == 0) throw std::invalid_argument("Worker count must be positive.");
       jobs = *requested;
       jobs_source = requested_source;
     } else if (const auto *environment = std::getenv("OMP_NUM_THREADS")) {
       jobs = parse_worker_count(environment, "OMP_NUM_THREADS");
       jobs_source = "OMP_NUM_THREADS";
     } else {
       jobs = 1;
       jobs_source = "default";
     }
     jobs_resolved = true;
   }

   void configure(const NumericalOptions &options) {
     interpolation_method = options.interpolation;
     algorithm = options.algorithm;
     epsabs = options.epsabs;
     epsrel = options.epsrel;
     workspace_limit = options.workspace_limit;
     quadrature_rule = options.quadrature_rule;
     gsl_error_policy = options.gsl_error_policy;
     if (algorithm == Algorithm::qag) {
       NRG::Tools::validate_tolerances(epsabs, epsrel);
       NRG::Tools::validate_qag_workspace_limit(workspace_limit);
     }
     if (options.jobs) {
       resolve_jobs(options.jobs, "NumericalOptions.jobs");
     } else {
       jobs = 1;
       jobs_source = "library default";
       jobs_resolved = true;
     }
   }

    void report_configuration() const {
      if (verbosity == 0) return;
      NRG::Tools::ConfigurationReport report("kk");
      report.value("verbosity", verbosity);
      report.resolved("mode", mode == MODE::FILES ? "files" : "stdio", "positional arguments");
       report.value("input", input_source);
       report.value("output", output_source);
       report.value("interpolation", NRG::Tools::interpolation_method_name(interpolation_method));
       report.value("algorithm", algorithm_name(algorithm));
       if (jobs_source == "--jobs")
         report.value("jobs", jobs);
       else
         report.resolved("jobs", jobs, jobs_source);
       report.value("jobs.source", jobs_source);
       if (algorithm == Algorithm::qag) {
         report.value("qag.epsabs", epsabs);
         report.value("qag.epsrel", epsrel);
         report.value("qag.workspace_limit", workspace_limit);
         report.value("qag.quadrature_rule", static_cast<int>(quadrature_rule));
         report.value("qag.error_policy", NRG::Tools::gsl_error_policy_name(gsl_error_policy));
       }
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
       if (!jobs_resolved) {
         jobs = 1;
         jobs_source = "library default";
         jobs_resolved = true;
       }
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
         auto sum = 0.0;
         if (algorithm == Algorithm::analytic) {
           polynomial.emplace(NRG::Tools::make_gsl_piecewise_polynomial(Xpts, Ypts, interpolation_method));
           sum = polynomial->integral();
         } else {
           const NRG::Tools::GslErrorHandlerGuard error_handler;
           qag_spline.reset(gsl_spline_alloc(NRG::Tools::gsl_interpolation_type(interpolation_method), len));
           if (!qag_spline) throw std::runtime_error("Failed to allocate GSL interpolation spline.");
           if (const auto status = gsl_spline_init(qag_spline.get(), Xpts.data(), Ypts.data(), len);
               status != GSL_SUCCESS)
             throw std::runtime_error("Failed to initialize GSL interpolation spline: "
                                      + std::string(gsl_strerror(status)));
           std::unique_ptr<gsl_interp_accel, gsl_accel_deleter> accelerator{gsl_interp_accel_alloc()};
           if (!accelerator) throw std::runtime_error("Failed to allocate GSL interpolation accelerator.");
           if (const auto status = gsl_spline_eval_integ_e(qag_spline.get(), Xmin, Xmax, accelerator.get(), &sum);
               status != GSL_SUCCESS)
             throw std::runtime_error("Failed to integrate GSL interpolation spline: "
                                      + std::string(gsl_strerror(status)));
         }
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
                  << "  -a, --algorithm ALGORITHM      qag or analytic (default: qag)\n"
                  << "  -i, --interpolation METHOD     linear, cspline, akima, or steffen (default: akima)\n"
                  << "  -j, --jobs N                   worker count (default: first OMP_NUM_THREADS value, or 1)\n"
                  << "      --epsabs VALUE             QAG absolute tolerance (default: 1e-12)\n"
                  << "      --epsrel VALUE             QAG relative tolerance (default: 1e-8)\n"
                  << "      --workspace-limit N        QAG integration workspace size (default: 1000)\n"
                  << "      --quadrature-rule RULE     QAG rule: 15, 21, 31, 41, 51, or 61 (default: 15)\n"
                  << "      --gsl-error-policy POLICY  QAG policy: ignore, warn, or fail (default: ignore)\n"
                  << "  -v                             show resolved configuration on standard error\n"
                 << "  -vv                            increase verbosity further\n"
                 << "  -V, --version                  show project version" << std::endl;
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
        {"algorithm", required_argument, nullptr, 'a'},
        {"interpolation", required_argument, nullptr, 'i'},
        {"jobs", required_argument, nullptr, 'j'},
        {"epsabs", required_argument, nullptr, OPT_EPSABS},
        {"epsrel", required_argument, nullptr, OPT_EPSREL},
        {"workspace-limit", required_argument, nullptr, OPT_WORKSPACE_LIMIT},
        {"quadrature-rule", required_argument, nullptr, OPT_QUADRATURE_RULE},
        {"gsl-error-policy", required_argument, nullptr, OPT_GSL_ERROR_POLICY},
        {nullptr, 0, nullptr, 0}
      };
      std::optional<size_t> requested_jobs;
      auto qag_control_requested = false;
      int option;
      while ((option = getopt_long(argc, argv, "ha:i:j:v", long_options, nullptr)) != -1) {
        switch (option) {
          case 'h':
             usage();
            std::exit(EXIT_SUCCESS);
          case 'a': algorithm = parse_algorithm(optarg); break;
          case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
          case 'j': requested_jobs = NRG::Tools::parse_positive_size(optarg, "Worker count"); break;
          case OPT_EPSABS:
            epsabs = NRG::Tools::parse_finite_double(optarg, "Absolute integration tolerance");
            qag_control_requested = true;
            break;
          case OPT_EPSREL:
            epsrel = NRG::Tools::parse_finite_double(optarg, "Relative integration tolerance");
            qag_control_requested = true;
            break;
          case OPT_WORKSPACE_LIMIT:
            workspace_limit = NRG::Tools::parse_positive_size(optarg, "Integration workspace limit");
            qag_control_requested = true;
            break;
          case OPT_QUADRATURE_RULE:
            quadrature_rule = NRG::Tools::parse_qag_rule(optarg);
            qag_control_requested = true;
            break;
          case OPT_GSL_ERROR_POLICY:
            gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg);
            qag_control_requested = true;
            break;
          case 'v': ++verbosity; break;
          default: throw std::invalid_argument("Unknown command-line option.");
        }
      }
      if (qag_control_requested && algorithm != Algorithm::qag)
        throw std::invalid_argument("QAG numerical controls require --algorithm qag.");
      if (algorithm == Algorithm::qag) {
        NRG::Tools::validate_tolerances(epsabs, epsrel);
        NRG::Tools::validate_qag_workspace_limit(workspace_limit);
      }
      resolve_jobs(requested_jobs, "--jobs");
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

    auto calc_analytic(const double argument) const {
      if (!polynomial) throw std::logic_error("Piecewise polynomial is not initialized.");
      const auto endpoint = argument == Xmin || argument == Xmax;
      const auto policy = endpoint ? NRG::Tools::CauchyEndpointPolicy::subtracted
                                   : NRG::Tools::CauchyEndpointPolicy::reject;
      const auto canonical = NRG::Tools::cauchy_principal_value(*polynomial, argument, policy);
      return -canonical.real() / M_PI;
    }

    auto calc_qag(const double argument, QagWorker &worker) const {
      if (!qag_spline) throw std::logic_error("QAG interpolation spline is not initialized.");
      if (!std::isfinite(argument)) throw std::invalid_argument("QAG argument must be finite.");
      if (argument < Xmin || argument > Xmax)
        throw std::domain_error("QAG arguments must lie inside the interpolation domain.");

      gsl_interp_accel_reset(worker.accelerator.get());
      double value_at_argument = std::numeric_limits<double>::quiet_NaN();
      if (const auto status = gsl_spline_eval_e(qag_spline.get(), argument, worker.accelerator.get(),
                                                &value_at_argument);
          status != GSL_SUCCESS)
        throw std::runtime_error("Failed to evaluate GSL interpolation spline: "
                                 + std::string(gsl_strerror(status)));

      const auto paired_symmetric_domain = Xmin == -Xmax;
      QagContext context{qag_spline.get(), worker.accelerator.get(), argument, value_at_argument,
                         paired_symmetric_domain};
      gsl_function integrand;
      integrand.function = &qag_integrand;
      integrand.params = &context;
      auto integral = std::numeric_limits<double>::quiet_NaN();
      auto integration_error = std::numeric_limits<double>::quiet_NaN();
      const auto integration_lower = paired_symmetric_domain ? 0.0 : Xmin;
      const auto status = gsl_integration_qag(&integrand, integration_lower, Xmax, epsabs, epsrel, workspace_limit,
                                               NRG::Tools::gsl_qag_rule(quadrature_rule), worker.workspace.get(),
                                               &integral, &integration_error);
      if (context.interpolation_status != GSL_SUCCESS)
        throw std::runtime_error("Failed to evaluate GSL interpolation spline during QAG: "
                                 + std::string(gsl_strerror(context.interpolation_status)));

      const auto failed = NRG::Tools::gsl_integration_failed(status, integral, integration_error);
      if (failed && gsl_error_policy != NRG::Tools::GslErrorPolicy::ignore) {
        std::ostringstream message;
        message << std::setprecision(17) << "GSL QAG failed for z=" << argument << ": status=" << status << " ("
                << gsl_strerror(status) << "), result=" << integral << ", estimated_error=" << integration_error
                << ", epsabs=" << epsabs << ", epsrel=" << epsrel << ", workspace_limit=" << workspace_limit
                << ", quadrature_rule=" << static_cast<int>(quadrature_rule);
        if (gsl_error_policy == NRG::Tools::GslErrorPolicy::warn) {
          static std::mutex warning_mutex;
          const std::lock_guard lock(warning_mutex);
          std::cerr << "kk: warning: " << message.str() << std::endl;
        } else {
          throw std::runtime_error(message.str());
        }
      }
      if (!std::isfinite(integral) || !std::isfinite(integration_error))
        throw std::overflow_error("QAG Cauchy transform is not finite.");

      auto correction = 0.0;
      if (argument != Xmin && argument != Xmax) {
        const auto extended_argument = static_cast<long double>(argument);
        const auto right_distance = static_cast<long double>(Xmax) - extended_argument;
        const auto left_distance = extended_argument - static_cast<long double>(Xmin);
        const auto smaller_distance = std::min(left_distance, right_distance);
        const auto larger_distance = std::max(left_distance, right_distance);
        const auto extended_logarithm = smaller_distance >= larger_distance / 2.0L
                                          ? std::log1p((static_cast<long double>(Xmin)
                                                        + static_cast<long double>(Xmax)
                                                        - 2.0L * extended_argument)
                                                       / left_distance)
                                          : std::log(right_distance) - std::log(left_distance);
        const auto logarithm = static_cast<double>(extended_logarithm);
        correction = value_at_argument * logarithm;
      }
      const auto result = (integral + correction) / M_PI;
      if (!std::isfinite(result)) throw std::overflow_error("QAG Cauchy transform is not finite.");
      return result;
    }

    void report_timing(const double wall_seconds, const std::optional<double> cpu_seconds) const {
      std::ostringstream report;
      report << std::setprecision(6) << "Time elapsed: " << wall_seconds << " s\n";
      if (verbosity != 0) {
        const auto workers = std::min(jobs, Xpts.size());
        const auto throughput = wall_seconds > 0.0 ? static_cast<double>(Xpts.size()) / wall_seconds : 0.0;
        report << "Performance: wall=" << wall_seconds << "s";
        if (cpu_seconds && wall_seconds > 0.0)
          report << " cpu=" << *cpu_seconds << "s effective_parallelism=" << *cpu_seconds / wall_seconds << "x";
        else
          report << " cpu=n/a effective_parallelism=n/a";
        report << " throughput=" << throughput << " points/s workers=" << workers
               << " algorithm=" << algorithm_name(algorithm) << '\n';
      }
      std::cout << report.str();
      std::cout.flush();
      if (!std::cout) throw std::runtime_error("<stdout>: timing report write failed.");
    }

    template<typename WorkerFactory, typename Calculate>
    auto calculate_grid(const DVEC &grid, WorkerFactory worker_factory, Calculate calculate) const {
      XYFUNC result(grid.size());
      if (grid.empty()) return result;

      const auto worker_count = std::min(jobs, grid.size());
      if (worker_count == 1) {
        auto worker = worker_factory();
        for (std::size_t index = 0; index < grid.size(); ++index)
          result[index] = {grid[index], calculate(grid[index], worker)};
        return result;
      }

      std::atomic<std::size_t> next_index{0};
      std::atomic<bool> stopped{false};
      std::exception_ptr failure;
      auto failure_index = std::numeric_limits<std::size_t>::max();
      std::mutex failure_mutex;
      auto run_worker = [&] {
        auto active_index = std::size_t{0};
        try {
          auto worker = worker_factory();
          while (!stopped.load(std::memory_order_relaxed)) {
            const auto index = next_index.fetch_add(1, std::memory_order_relaxed);
            if (index >= grid.size()) break;
            active_index = index;
            result[index] = {grid[index], calculate(grid[index], worker)};
          }
        } catch (...) {
          {
            const std::lock_guard lock(failure_mutex);
            if (!failure || active_index < failure_index) {
              failure = std::current_exception();
              failure_index = active_index;
            }
          }
          stopped.store(true, std::memory_order_relaxed);
        }
      };

      std::vector<std::thread> workers;
      workers.reserve(worker_count);
      try {
        for (std::size_t worker = 0; worker < worker_count; ++worker) workers.emplace_back(run_worker);
      } catch (...) {
        stopped.store(true, std::memory_order_relaxed);
        for (auto &worker : workers) worker.join();
        throw;
      }
      for (auto &worker : workers) worker.join();
      if (failure) std::rethrow_exception(failure);
      return result;
    }

  public:
    auto calc(const double argument) const {
      if (algorithm == Algorithm::analytic) return calc_analytic(argument);
      const NRG::Tools::GslErrorHandlerGuard error_handler;
      QagWorker worker(workspace_limit);
      return calc_qag(argument, worker);
    }
   
    // Perform the calculations for all points on a grid
    auto calc(const DVEC &grid) const {
      if (algorithm == Algorithm::analytic)
        return calculate_grid(grid, [] { return AnalyticWorker{}; },
                              [this](const double argument, AnalyticWorker &) {
                                return calc_analytic(argument);
                              });
      const NRG::Tools::GslErrorHandlerGuard error_handler;
      return calculate_grid(grid, [this] { return QagWorker{workspace_limit}; },
                            [this](const double argument, QagWorker &worker) {
                              return calc_qag(argument, worker);
                            });
    }

    // Legacy interface when kk is used as a command-line tool
    KK(int argv, char *argc[]) {
      const auto wall_start = std::chrono::steady_clock::now();
      const auto cpu_start = std::clock();
      parse_cmd_line(argv, argc);
      const auto im = read(mode == MODE::FILES ? Fin : std::cin, input_source);
      init(im);
      const auto re = calc(Xpts);
      if (mode == MODE::FILES) {
        std::ostringstream output;
        write(re, output, OUTPUT_PRECISION, output_source);
        NRG::Tools::write_output_file(output_source, output.str());
        const auto wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
        const auto cpu_end = std::clock();
        std::optional<double> cpu_seconds;
        if (cpu_start != static_cast<std::clock_t>(-1) && cpu_end != static_cast<std::clock_t>(-1)
            && cpu_end >= cpu_start)
          cpu_seconds = static_cast<double>(cpu_end - cpu_start) / static_cast<double>(CLOCKS_PER_SEC);
        report_timing(wall_seconds, cpu_seconds);
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
