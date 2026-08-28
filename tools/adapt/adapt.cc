// Discretization ODE solver for NRG
// Adaptable discretization mesh code
//
// Linear interpolation of rho(omega), integration using trapezoidal method. 4th-order Runge-Kutta ODE solver. Secant
// method for refinement of parameter A.

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

using namespace std::string_literals;

#include "adapt.hpp"
using namespace NRG::Adapt;

void about(std::ostream &F = std::cout) {
  F << "# Discretization ODE solver" << std::endl;
  F << "# Rok Zitko, rok.zitko@ijs.si, 2008-2020" << std::endl;
}

void help(int argc, char **argv, const std::string &help_message)
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() >= 1 && args[0] == "-h") {
    std::cout << help_message << std::endl;
    std::exit(EXIT_SUCCESS);
  }
}

const auto usage = "Usage: adapt [-h|--help] [--flat GG] [--integral] [--epsabs VALUE] [--epsrel VALUE] "
                   "[--workspace-limit N] [--gsl-error-policy ignore|warn|fail] [P|N] [param_filename]"s;

struct CommandLineOptions {
  Sign sign = Sign::POS;
  std::string param_fn = "param";
  std::optional<double> flat_gamma;
  bool integral = false;
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

int main(int argc, char *argv[]) {
  try {
    const std::clock_t start_clock = std::clock();
    about();
    help(argc, argv, usage);
    const auto options = cmd_line(argc, argv);
    Params P(options.param_fn);
    Adapt calc(P, options.sign, options.flat_gamma, options.integral, options.cquad);
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
