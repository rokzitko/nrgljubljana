// Channel-mixing discretization for NRG
//
// Discretizes a matrix hybridisation function Gamma(omega) into a star Hamiltonian, following
// J.-G. Liu, D. Wang and Q.-H. Wang, PRB 93, 035102 (2016). See README.md.

#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <limits>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <common/version.hpp>

#include "../common/diagnostics.hpp"
#include "../common/gsl_config.hpp"
#include "../common/tabulated_density.hpp"
#include "branches.hpp"
#include "chain.hpp"
#include "chain_io.hpp"
#include "load.hpp"
#include "mesh.hpp"
#include "parser.hpp"
#include "precision.hpp"
#include "star.hpp"
#include "star_io.hpp"
#include "types.hpp"

using namespace NRG::MixChain;

namespace {

enum class Mode {
  Full, // star and chain
  Star, // 's': compute and save the star
  Chain // 'l': load the star and tridiagonalize
};

constexpr auto usage_text = "Usage: mixchain [options] [s|l] [parameter_file]";

void about(std::ostream &F = std::cout) { F << "# Discretization of a matrix hybridisation function" << std::endl; }

void usage(std::ostream &F = std::cout) {
  F << usage_text << std::endl;
  F << " -h, --help -- show help" << std::endl;
  F << " -v -- show resolved configuration and verbose diagnostics on stderr" << std::endl;
  F << " -vv -- increase verbosity further" << std::endl;
  F << " -V, --version -- show project version and exit" << std::endl;
  F << " --epsabs VALUE -- absolute tolerance of the integral method" << std::endl;
  F << " --epsrel VALUE -- relative tolerance of the integral method" << std::endl;
  F << " --workspace-limit VALUE -- size of the integration workspace" << std::endl;
  F << " --gsl-error-policy fail|warn|ignore -- how to report integration failures" << std::endl;
  F << " s -- discretize Gamma and write the star to star.dat" << std::endl;
  F << " l -- read the star from star.dat, tridiagonalize, and write the chain to chain.dat" << std::endl;
  F << " (no mode) -- both, one after the other" << std::endl;
}

struct CommandLineOptions {
  std::string param_filename = "param";
  Mode mode                  = Mode::Full;
  int verbosity              = 0;
  NRG::Tools::CquadOptions cquad;
};

bool matches_value_option(const std::string &arg, const std::string &option) {
  return arg == option || arg.starts_with(option + "=");
}

std::string value_for_option(const std::string &arg, const std::string &option, int &index, const int argc,
                             char *argv[]) {
  if (arg == option) {
    if (index + 1 >= argc) throw std::invalid_argument("Missing value for " + option + ".\n" + usage_text);
    return argv[++index];
  }
  return arg.substr(option.size() + 1);
}

CommandLineOptions cmd_line(int argc, char *argv[]) {
  CommandLineOptions options;
  bool mode_set  = false;
  bool param_set = false;

  for (int i = 1; i < argc; i++) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      exit(EXIT_SUCCESS);
    }
    if (arg.size() >= 2 && arg[0] == '-'
        && std::all_of(arg.begin() + 1, arg.end(), [](const char ch) { return ch == 'v'; })) {
      options.verbosity += static_cast<int>(arg.size() - 1);
      continue;
    }
    if (matches_value_option(arg, "--epsabs")) {
      if (options.cquad.epsabs) throw std::invalid_argument("--epsabs specified more than once.\n" + std::string(usage_text));
      options.cquad.epsabs = NRG::Tools::parse_finite_double(value_for_option(arg, "--epsabs", i, argc, argv),
                                                             "Absolute integration tolerance");
      continue;
    }
    if (matches_value_option(arg, "--epsrel")) {
      if (options.cquad.epsrel) throw std::invalid_argument("--epsrel specified more than once.\n" + std::string(usage_text));
      options.cquad.epsrel = NRG::Tools::parse_finite_double(value_for_option(arg, "--epsrel", i, argc, argv),
                                                             "Relative integration tolerance");
      continue;
    }
    if (matches_value_option(arg, "--workspace-limit")) {
      if (options.cquad.workspace_limit)
        throw std::invalid_argument("--workspace-limit specified more than once.\n" + std::string(usage_text));
      options.cquad.workspace_limit = NRG::Tools::parse_positive_size(
        value_for_option(arg, "--workspace-limit", i, argc, argv), "Integration workspace limit");
      continue;
    }
    if (matches_value_option(arg, "--gsl-error-policy")) {
      if (options.cquad.gsl_error_policy)
        throw std::invalid_argument("--gsl-error-policy specified more than once.\n" + std::string(usage_text));
      options.cquad.gsl_error_policy =
        NRG::Tools::parse_gsl_error_policy(value_for_option(arg, "--gsl-error-policy", i, argc, argv));
      continue;
    }
    if (arg == "s" || arg == "l") {
      if (mode_set) throw std::invalid_argument("Mode specified more than once.\n" + std::string(usage_text));
      options.mode = arg == "s" ? Mode::Star : Mode::Chain;
      mode_set     = true;
      continue;
    }
    if (!arg.empty() && arg[0] == '-')
      throw std::invalid_argument("Unknown option: " + arg + "\n" + usage_text);
    if (param_set) throw std::invalid_argument("Unexpected argument: " + arg + "\n" + usage_text);
    options.param_filename = arg;
    param_set              = true;
  }

  if (options.cquad.epsrel)
    NRG::Tools::validate_cquad_tolerances(options.cquad.epsabs.value_or(0.0), *options.cquad.epsrel);
  else if (options.cquad.epsabs)
    NRG::Tools::validate_cquad_tolerances(*options.cquad.epsabs, 1.0);
  if (options.cquad.workspace_limit) NRG::Tools::validate_cquad_workspace_limit(*options.cquad.workspace_limit);

  return options;
}

struct Configuration {
  GammaOptions gamma;
  StarOptions star;
  bool mmax_from_nmax{};
  ChainOptions chain;
  unsigned int preccpp{};
};

auto builds_star(const Mode mode) { return mode != Mode::Chain; }
auto builds_chain(const Mode mode) { return mode != Mode::Star; }

void read_star_configuration(const Params &P, const CommandLineOptions &command_line, Configuration &configuration) {
  configuration.gamma.prefix                = P.Pstr("dos_prefix", "Gamma");
  configuration.gamma.channels              = P.Pint("channels", 1);
  configuration.gamma.bandrescale           = P.P("bandrescale", 1.0);
  configuration.gamma.hermiticity_tolerance = P.P("hermiticity_tolerance", 1e-8);

  auto &star       = configuration.star;
  star.Lambda      = NRG::Tools::LambdaCache(P.P("Lambda", 2.0));
  star.z           = P.P("z", 1.0);
  star.bandrescale = configuration.gamma.bandrescale;
  star.adapt       = P.Pbool("adapt", false);
  star.hardgap     = P.Pbool("hardgap", false);
  star.boundary    = P.P("boundary", 0.0);
  star.mesh_weight = mesh_weight_from_string(P.Pstr("mesh_weight", "frobenius"));
  star.interpolation =
    NRG::Tools::parse_density_interpolation_method(P.Pstr("density_interpolation", "linear"));
  star.branches.ordering = branch_ordering_from_string(P.Pstr("branch_ordering", "tracked"));
  star.allowed_error     = P.P("allowed_error", 1e-10);
  star.cquad             = command_line.cquad;

  // The star stage needs mMAX, the number of intervals. Nmax is the length of the Wilson chain, which only the chain
  // stage uses; it serves here as the customary default for mMAX.
  if (P.contains("mMAX")) {
    const auto mmax = P.Pint("mMAX", 0);
    if (mmax <= 0) throw std::invalid_argument("mMAX must be greater than 0.");
    star.mMAX = static_cast<unsigned int>(mmax);
  } else if (P.contains("Nmax")) {
    const auto nmax = P.Pint("Nmax", 0);
    if (nmax <= 0) throw std::invalid_argument("Nmax must be greater than 0.");
    star.mMAX                  = static_cast<unsigned int>(2 * nmax);
    configuration.mmax_from_nmax = true;
  } else {
    throw std::invalid_argument("Either mMAX or Nmax must be given.");
  }
}

void read_chain_configuration(const Params &P, Configuration &configuration) {
  if (!P.contains("Nmax")) throw std::invalid_argument("Nmax must be given to build the chain.");
  const auto nmax = P.Pint("Nmax", 0);
  if (nmax <= 0) throw std::invalid_argument("Nmax must be greater than 0.");
  configuration.chain.Nmax = static_cast<unsigned int>(nmax);

  configuration.chain.breakdown_tolerance = P.P("breakdown_tolerance", 1e-20);
  if (!(std::isfinite(configuration.chain.breakdown_tolerance) && configuration.chain.breakdown_tolerance > 0.0))
    throw std::invalid_argument("breakdown_tolerance must be a positive finite number.");

  // As in nrgchain, in bits. It is rounded up to the precision ladder of precision.hpp.
  const auto preccpp = P.Pint("preccpp", 2000);
  if (preccpp <= 10) throw std::invalid_argument("preccpp must be greater than 10.");
  configuration.preccpp = static_cast<unsigned int>(preccpp);
  resolve_precision(configuration.preccpp); // fail before the star stage runs, not after it
}

Configuration read_configuration(const Params &P, const CommandLineOptions &command_line) {
  Configuration configuration;
  if (builds_star(command_line.mode)) read_star_configuration(P, command_line, configuration);
  if (builds_chain(command_line.mode)) read_chain_configuration(P, configuration);
  return configuration;
}

const char *mode_name(const Mode mode) {
  switch (mode) {
    case Mode::Full: return "star-and-chain";
    case Mode::Star: return "star";
    case Mode::Chain: return "chain";
  }
  return "unknown";
}

void report_star_configuration(const Configuration &configuration, NRG::Tools::ConfigurationReport &report);

void report_configuration(const Configuration &configuration, const CommandLineOptions &command_line) {
  if (command_line.verbosity == 0) return;
  NRG::Tools::ConfigurationReport report("mixchain");
  report.value("verbosity", command_line.verbosity);
  report.value("parameter_file", command_line.param_filename);
  report.value("mode", mode_name(command_line.mode));
  if (builds_star(command_line.mode)) report_star_configuration(configuration, report);
  if (builds_chain(command_line.mode)) {
    report.value("Nmax", configuration.chain.Nmax);
    report.value("preccpp", configuration.preccpp);
    report.resolved("digits", resolve_precision(configuration.preccpp), "smallest precision rung covering preccpp");
    report.value("breakdown_tolerance", configuration.chain.breakdown_tolerance);
  }
  report.write(std::cerr);
}

void report_star_configuration(const Configuration &configuration, NRG::Tools::ConfigurationReport &report) {
  const auto &star = configuration.star;
  report.value("channels", configuration.gamma.channels);
  report.value("dos_prefix", configuration.gamma.prefix);
  report.value("Lambda", static_cast<double>(star.Lambda));
  report.value("z", star.z);
  if (configuration.mmax_from_nmax)
    report.resolved("mMAX", star.mMAX, "2*Nmax");
  else
    report.value("mMAX", star.mMAX);
  report.value("bandrescale", star.bandrescale);
  report.value("adapt", star.adapt);
  report.value("hardgap", star.hardgap);
  // boundary is a fraction of the rescaled band edge, as in adapt; its value in the units of the input file is
  // boundary*bandrescale.
  report.value("boundary", star.boundary);
  report.value("boundary_in_input_units", star.boundary * star.bandrescale);
  report.value("mesh_weight", star.adapt ? mesh_weight_name(star.mesh_weight) : std::string("inactive"));
  report.value("density_interpolation", NRG::Tools::interpolation_method_name(star.interpolation));
  report.value("branch_ordering", branch_ordering_name(star.branches.ordering));
  report.value("allowed_error", star.allowed_error);
  report.value("hermiticity_tolerance", configuration.gamma.hermiticity_tolerance);
}

template<typename S> void report_star(const Star<S> &star, std::ostream &out) {
  const auto &diagnostics = star.diagnostics;
  out << "# levels=" << star.levels.size() << " complex=" << (is_complex_v<S> ? 1 : 0) << std::endl;
  out << "# max_interval_deviation=" << diagnostics.max_interval_deviation
      << " at omega=" << diagnostics.max_interval_omega << std::endl;
  out << "# max_cquad_error=" << diagnostics.max_cquad_error << std::endl;
  out << "# crossings: " << diagnostics.crossings_pos.size() << " positive, " << diagnostics.crossings_neg.size()
      << " negative" << std::endl;

  // Where the mesh reaches below the innermost tabulated frequency, the density is the constant continuation of the
  // input: exact for a flat band, an approximation for anything with structure at low frequency.
  // These frequencies are compared with the input grid and can sit very close to an accumulation point, so they are
  // printed with every digit that distinguishes them.
  const auto precision = out.precision(std::numeric_limits<double>::max_digits10);
  for (const auto &[name, coverage] : {std::pair{"POS", &diagnostics.coverage_pos},
                                       std::pair{"NEG", &diagnostics.coverage_neg}}) {
    if (coverage->collapsed_levels > 0)
      out << "# " << name << ": " << coverage->collapsed_levels
          << " representative energy levels are indistinguishable from the accumulation point "
          << coverage->accumulation_point << " in double precision" << std::endl;
    if (coverage->unresolved_intervals == 0) continue;
    out << "# " << name << ": " << coverage->unresolved_intervals << " of " << star.mMAX + 1
        << " intervals contain no tabulated point of the input, the outermost being [" << coverage->unresolved_to
        << ", " << coverage->unresolved_from << "]; ";
    if (coverage->continued())
      out << "the input ends at omega=" << coverage->innermost_input
          << " and below that the density is its constant continuation" << std::endl;
    else
      out << "there the star follows the interpolant between neighbouring points" << std::endl;
  }
  out.precision(precision);

  // The trace of Theta agrees with that of the integral of Gamma by construction; the off-diagonal elements only
  // approximately, because a single level per interval and branch cannot follow a rotating eigenvector.
  const auto trace_star  = std::real(star.theta.trace());
  const auto trace_exact = std::real(star.theta_exact.trace());
  out << "# tr(theta)=" << trace_star << " tr(int Gamma)=" << trace_exact << " difference="
      << trace_star - trace_exact << std::endl;
  const auto scale = star.theta_exact.norm();
  if (scale > 0.0)
    out << "# ||theta - int Gamma||/||int Gamma||=" << (star.theta - star.theta_exact).norm() / scale << std::endl;
}

template<typename S> void run_star(const Configuration &configuration) {
  const auto input = load_gamma<S>(configuration.gamma);
  const auto star  = build_star(input, configuration.star);
  report_star(star, std::cout);
  save_star(star, star_default_filename);
  std::cout << "# star written to " << star_default_filename << std::endl;
}

// The star is self-contained: the chain is built with the Lambda, z and bandrescale it records. A parameter file that
// sets any of them to something else was meant for a different star, so that is an error rather than a silent choice.
template<typename S> void check_star_against_parameters(const Star<S> &star, const Params &P) {
  const auto check = [&P](const std::string &name, const double recorded) {
    if (!P.contains(name)) return;
    const auto requested = P.P(name, recorded);
    if (std::abs(requested - recorded) <= 1e-12 * std::max(1.0, std::abs(recorded))) return;
    std::ostringstream message;
    message << std::setprecision(17) << "The star in " << star_default_filename << " was built with " << name << "="
            << recorded << ", but the parameter file gives " << name << "=" << requested << ".";
    throw std::invalid_argument(message.str());
  };
  check("Lambda", star.Lambda);
  check("z", star.z);
  check("bandrescale", star.bandrescale);
}

template<typename S> void report_chain(const Chain<S> &chain, const unsigned digits, std::ostream &out) {
  const auto &d = chain.diagnostics;
  out << "# chain: sites=" << chain.Nmax + 1 << " channels=" << chain.channels << " digits=" << digits << std::endl;
  out << "# theta_condition=" << d.theta_condition << " min_residual_condition=" << d.min_residual_condition
      << std::endl;
  out << "# max_antihermitian=" << d.max_antihermitian << " max_reorthogonalization=" << d.max_reorthogonalization
      << std::endl;
}

// The chain is built from star.dat also in the default mode, right after the star stage has written it, so that the
// default mode and 's' followed by 'l' produce the same chain by construction.
template<typename S0> void run_chain(const Configuration &configuration, const Params &P) {
  const auto star = load_star<S0>(star_default_filename);
  check_star_against_parameters(star, P);
  const auto digits = resolve_precision(configuration.preccpp);
  with_precision_like<S0>(configuration.preccpp, [&]<typename S>() {
    const auto chain = build_chain<S>(star, configuration.chain);
    report_chain(chain, digits, std::cout);
    save_chain(chain, ChainFileHeader{star.z, star.Lambda, star.bandrescale, digits}, chain_default_filename);
  });
  std::cout << "# chain written to " << chain_default_filename << std::endl;
}

// Wall-clock time of one stage, for comparing the cost of the two.
template<typename F> void timed(const char *stage, F &&f) {
  const auto start = std::chrono::steady_clock::now();
  f();
  const auto seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
  std::cout << "# " << stage << " stage: " << seconds << " s" << std::endl;
}

void run(const CommandLineOptions &command_line) {
  const Params P(command_line.param_filename);
  const auto configuration = read_configuration(P, command_line);
  report_configuration(configuration, command_line);

  if (builds_star(command_line.mode))
    timed("star", [&] {
      if (gamma_is_complex(configuration.gamma))
        run_star<std::complex<double>>(configuration);
      else
        run_star<double>(configuration);
    });
  if (builds_chain(command_line.mode))
    timed("chain", [&] {
      if (star_is_complex(star_default_filename))
        run_chain<std::complex<double>>(configuration, P);
      else
        run_chain<double>(configuration, P);
    });
}

} // namespace

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "mixchain")) return EXIT_SUCCESS;
  try {
    const clock_t start_clock = clock();
    about();
    run(cmd_line(argc, argv));
    const clock_t end_clock = clock();
    std::cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "mixchain: error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
