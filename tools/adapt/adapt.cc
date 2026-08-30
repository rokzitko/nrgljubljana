// Discretization ODE solver for NRG
// Adaptable discretization mesh code
//
// Selectable shape-preserving interpolation of rho(omega) with an exact interval primitive. 4th-order Runge-Kutta
// ODE solver. Secant method for refinement of parameter A.

#include <cstddef>
#include <exception>
#include <iostream>
#include <ostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cctype>
#include <cmath>
#include <optional>
#include <stdexcept>

#include <common/version.hpp>

using namespace std::string_literals;

#include "adapt.hpp"
#include "../common/diagnostics.hpp"
using namespace NRG::Adapt;

void about(std::ostream &F = std::cout) {
  F << "# Discretization ODE solver" << std::endl;
}

void help(int argc, char **argv, const std::string &help_message)
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() >= 1 && args[0] == "-h") {
    std::cout << help_message << std::endl;
    std::exit(EXIT_SUCCESS);
  }
}

const auto usage = "Usage: adapt [-h|--help] [-v|-vv] [--flat GG] [--integral] [--epsabs VALUE] [--epsrel VALUE] "
                   "[--workspace-limit N] [--gsl-error-policy ignore|warn|fail] [P|N] [param_filename]\n"
                   "  -v              print resolved configuration to stderr\n"
                   "  -vv             enable very verbose diagnostics\n"
                   "  -V, --version   print version and exit"s;

struct CommandLineOptions {
  Sign sign = Sign::POS;
  std::string param_fn = "param";
  std::optional<double> flat_gamma;
  bool integral = false;
  bool verbose = false;
  bool veryverbose = false;
  CquadOptions cquad;
};

auto uppercase(std::string text) {
  for (auto &ch : text) { ch = static_cast<char>(std::toupper(static_cast<unsigned char>(ch))); }
  return text;
}

bool parse_sign_arg(const std::string &arg, Sign &sign) {
  const auto token = uppercase(arg);
  if (token == "P" || token == "POS" || token == "POSITIVE") {
    sign = Sign::POS;
    return true;
  }
  if (token == "N" || token == "NEG" || token == "NEGATIVE") {
    sign = Sign::NEG;
    return true;
  }
  return false;
}

double parse_flat_gamma(const std::string &value) {
  std::size_t parsed = 0;
  double gamma       = 0.0;
  try {
    gamma = std::stod(value, &parsed);
  } catch (const std::exception &) {
    throw std::invalid_argument("--flat expects a positive finite number.");
  }
  if (parsed != value.size() || !(std::isfinite(gamma) && gamma > 0.0)) {
    throw std::invalid_argument("--flat expects a positive finite number.");
  }
  return gamma;
}

bool matches_value_option(const std::string &arg, const std::string &option) {
  return arg == option || arg.starts_with(option + "=");
}

std::string value_for_option(const std::string &arg,
                             const std::string &option,
                             int &index,
                             const int argc,
                             char *argv[]) {
  if (arg == option) {
    if (index + 1 >= argc) { throw std::invalid_argument("Missing value for " + option + ".\n" + usage); }
    return argv[++index];
  }
  return arg.substr(option.size() + 1);
}

CommandLineOptions cmd_line(int argc, char *argv[]) {
  CommandLineOptions options;
  bool sign_set  = false;
  bool param_set = false;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      std::cout << usage << std::endl;
      std::exit(EXIT_SUCCESS);
    }
    if (arg == "-v") {
      if (options.verbose) options.veryverbose = true;
      options.verbose = true;
      continue;
    }
    if (arg == "-vv") {
      options.verbose     = true;
      options.veryverbose = true;
      continue;
    }
    if (arg == "--flat") {
      if (options.flat_gamma) { throw std::invalid_argument("--flat specified more than once.\n" + usage); }
      if (i + 1 >= argc) { throw std::invalid_argument("Missing value for --flat.\n" + usage); }
      options.flat_gamma = parse_flat_gamma(argv[++i]);
      continue;
    }
    if (arg.starts_with("--flat=")) {
      if (options.flat_gamma) { throw std::invalid_argument("--flat specified more than once.\n" + usage); }
      options.flat_gamma = parse_flat_gamma(arg.substr(7));
      continue;
    }
    if (arg == "--integral") {
      if (options.integral) { throw std::invalid_argument("--integral specified more than once.\n" + usage); }
      options.integral = true;
      continue;
    }
    if (matches_value_option(arg, "--epsabs")) {
      if (options.cquad.epsabs) { throw std::invalid_argument("--epsabs specified more than once.\n" + usage); }
      options.cquad.epsabs = NRG::Tools::parse_finite_double(
        value_for_option(arg, "--epsabs", i, argc, argv), "Absolute integration tolerance");
      continue;
    }
    if (matches_value_option(arg, "--epsrel")) {
      if (options.cquad.epsrel) { throw std::invalid_argument("--epsrel specified more than once.\n" + usage); }
      options.cquad.epsrel = NRG::Tools::parse_finite_double(
        value_for_option(arg, "--epsrel", i, argc, argv), "Relative integration tolerance");
      continue;
    }
    if (matches_value_option(arg, "--workspace-limit")) {
      if (options.cquad.workspace_limit) {
        throw std::invalid_argument("--workspace-limit specified more than once.\n" + usage);
      }
      options.cquad.workspace_limit = NRG::Tools::parse_positive_size(
        value_for_option(arg, "--workspace-limit", i, argc, argv), "Integration workspace limit");
      continue;
    }
    if (matches_value_option(arg, "--gsl-error-policy")) {
      if (options.cquad.gsl_error_policy) {
        throw std::invalid_argument("--gsl-error-policy specified more than once.\n" + usage);
      }
      options.cquad.gsl_error_policy = NRG::Tools::parse_gsl_error_policy(
        value_for_option(arg, "--gsl-error-policy", i, argc, argv));
      continue;
    }
    if (!arg.empty() && arg[0] == '-') { throw std::invalid_argument("Unknown option: " + arg + "\n" + usage); }

    Sign parsed_sign;
    if (parse_sign_arg(arg, parsed_sign)) {
      if (sign_set) { throw std::invalid_argument("Sign specified more than once.\n" + usage); }
      options.sign = parsed_sign;
      sign_set     = true;
      continue;
    }

    if (param_set) { throw std::invalid_argument("Unexpected argument: " + arg + "\n" + usage); }
    options.param_fn = arg;
    param_set        = true;
  }

  if (options.cquad.epsrel) {
    NRG::Tools::validate_cquad_tolerances(options.cquad.epsabs.value_or(0.0), *options.cquad.epsrel);
  } else if (options.cquad.epsabs) {
    NRG::Tools::validate_cquad_tolerances(*options.cquad.epsabs, 1.0);
  }
  if (options.cquad.workspace_limit) NRG::Tools::validate_cquad_workspace_limit(*options.cquad.workspace_limit);

  std::cout << "# ++ " << (options.sign == Sign::POS ? "POSITIVE" : "NEGATIVE") << std::endl;
  return options;
}

void report_configuration(const CommandLineOptions &options, const Adapt &calc) {
  NRG::Tools::ConfigurationReport report("adapt");
  report.value("verbosity", options.veryverbose ? 2 : 1);
  report.value("frequency_branch", calc.sign == Sign::POS ? "positive" : "negative");
  report.value("parameter_file", options.param_fn);
  if (calc.flat_gamma) {
    report.value("density_mode", "flat");
    report.value("flat_gamma", *calc.flat_gamma);
  } else {
    report.value("density_mode", "file");
    const auto density_file = calc.P.Pstr("dos", "Delta.dat");
    if (calc.P.contains("dos")) {
      report.value("density_file", density_file);
    } else {
      report.resolved("density_file", density_file, "parameter default");
    }
  }
  const auto density_method = NRG::Tools::interpolation_method_name(calc.density_interpolation);
  if (calc.P.contains("density_interpolation"))
    report.value("density_interpolation", density_method);
  else
    report.resolved("density_interpolation", density_method, "parameter default");
  report.value("density_integration", "exact interpolant primitive");
  report.value("frequency_min", 0.0);
  report.value("frequency_max", 1.0);
  report.value("Lambda", static_cast<double>(calc.Lambda));
  report.value("adaptive_mesh", calc.adapt);
  report.value("hardgap", calc.hardgap);
  report.value("boundary", calc.boundary);
  report.value("bandrescale", calc.bandrescale);
  report.value("x_min", 1.0);
  report.value("xmax", calc.xmax);
  report.value("xfine", calc.xfine);
  report.value("output_step", calc.output_step);
  report.value("dx_fine", calc.dx_fine);
  report.value("dx_fast", calc.dx_fast);
  report.value("allowed_error", calc.allowed_error);
  report.value("max_subdiv", calc.max_subdiv);
  report.value("max_abs", calc.max_abs);
  report.value("secant_eps", calc.convergence_eps);
  report.value("secant_factor", calc.P.P("secant_factor", 1e-7));
  report.value("secant_max_iter", calc.max_iter);
  if (!calc.adapt) {
    report.resolved("load_g", "inactive", "adaptive_mesh=false");
    report.resolved("g_file", "inactive", "adaptive_mesh=false");
  } else {
    const auto load_g = calc.P.Pbool("loadg", false);
    if (calc.P.contains("loadg")) {
      report.value("load_g", load_g);
    } else {
      report.resolved("load_g", load_g, "parameter default");
    }
    report.value("g_file", calc.sign == Sign::POS ? "GSOL.dat" : "GSOLNEG.dat");
  }
  report.value("f_file", calc.sign == Sign::POS ? "FSOL.dat" : "FSOLNEG.dat");

  const auto method = calc.f_method == FMethod::ODE ? "ode" : "integral";
  if (options.integral) {
    report.value("f_method", method);
    report.value("f_method_source", "--integral override");
  } else if (calc.P.contains("f_method")) {
    report.value("f_method", method);
  } else {
    report.resolved("f_method", method, "parameter default");
  }
  report.value("cquad_active", calc.f_method == FMethod::INTEGRAL);
  if (options.cquad.epsabs) {
    report.value("cquad_epsabs", *options.cquad.epsabs);
  } else {
    report.resolved("cquad_epsabs", 0.0, "CQUAD default");
  }
  if (options.cquad.epsrel) {
    report.value("cquad_epsrel", *options.cquad.epsrel);
  } else {
    report.resolved("cquad_epsrel", calc.allowed_error, "allowed_error");
  }
  if (options.cquad.workspace_limit) {
    report.value("cquad_workspace_limit", *options.cquad.workspace_limit);
  } else {
    report.resolved("cquad_workspace_limit", std::size_t{1000}, "CQUAD default");
  }
  const auto error_policy = options.cquad.gsl_error_policy.value_or(NRG::Tools::GslErrorPolicy::fail);
  if (options.cquad.gsl_error_policy) {
    report.value("cquad_error_policy", NRG::Tools::gsl_error_policy_name(error_policy));
  } else {
    report.resolved("cquad_error_policy", NRG::Tools::gsl_error_policy_name(error_policy), "default policy");
  }
  report.write(std::cerr);
}

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "adapt")) return EXIT_SUCCESS;
  try {
    const std::clock_t start_clock = std::clock();
    about();
    help(argc, argv, usage);
    const auto options = cmd_line(argc, argv);
    Params P(options.param_fn);
    Adapt calc(P, options.sign, options.flat_gamma, options.integral, options.cquad);
    if (options.verbose) report_configuration(options, calc);
    calc.run();
    const std::clock_t end_clock = std::clock();
    std::cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "adapt: error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "adapt: error: unknown exception" << std::endl;
    return EXIT_FAILURE;
  }
}
