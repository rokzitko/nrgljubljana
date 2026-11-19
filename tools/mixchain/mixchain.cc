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
#include <filesystem>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <common/version.hpp>

#include "../common/diagnostics.hpp"
#include "../common/gsl_config.hpp"
#include "../common/tabulated_density.hpp"
#include "blocks.hpp"
#include "branches.hpp"
#include "chain.hpp"
#include "chain_io.hpp"
#include "configuration.hpp"
#include "discretize.hpp"
#include "load.hpp"
#include "mesh.hpp"
#include "parser.hpp"
#include "report.hpp"
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
  F << " --Nz N -- discretize for z = i/N, i = 1..N, into the directories 1/ .. N/; z in the parameter file is ignored"
    << std::endl;
  F << " s -- discretize Gamma and write the star to star.dat" << std::endl;
  F << " l -- read the star from star.dat, tridiagonalize, and write the chain to chain.dat" << std::endl;
  F << " (no mode) -- both, one after the other" << std::endl;
}

struct CommandLineOptions {
  std::string param_filename = "param";
  Mode mode                  = Mode::Full;
  int verbosity              = 0;
  NRG::Tools::CquadOptions cquad;
  std::optional<unsigned int> Nz; // several values of z, each in its own directory
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
    if (matches_value_option(arg, "--Nz")) {
      if (options.Nz) throw std::invalid_argument("--Nz specified more than once.\n" + std::string(usage_text));
      const auto value = value_for_option(arg, "--Nz", i, argc, argv);
      const auto count = NRG::Tools::parse_parameter_int(value, "--Nz");
      if (count <= 0) throw std::invalid_argument("--Nz must be a positive integer: " + value);
      options.Nz = static_cast<unsigned int>(count);
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

auto builds_star(const Mode mode) { return mode != Mode::Chain; }
auto builds_chain(const Mode mode) { return mode != Mode::Star; }

const char *mode_name(const Mode mode) {
  switch (mode) {
    case Mode::Full: return "star-and-chain";
    case Mode::Star: return "star";
    case Mode::Chain: return "chain";
  }
  return "unknown";
}

void report_star_configuration(const Configuration &configuration, const CommandLineOptions &command_line,
                               NRG::Tools::ConfigurationReport &report);

void report_configuration(const Configuration &configuration, const CommandLineOptions &command_line) {
  if (command_line.verbosity == 0) return;
  NRG::Tools::ConfigurationReport report("mixchain");
  report.value("verbosity", command_line.verbosity);
  report.value("parameter_file", command_line.param_filename);
  report.value("mode", mode_name(command_line.mode));
  if (command_line.Nz) report.value("Nz", *command_line.Nz);
  if (builds_star(command_line.mode)) report_star_configuration(configuration, command_line, report);
  if (builds_chain(command_line.mode)) {
    report.value("Nmax", configuration.chain.Nmax);
    report.value("rank_tolerance", configuration.chain.rank_tolerance);
    report.value("sensitivity_tolerance", configuration.chain.sensitivity_tolerance);
    report.value("discretization_files", configuration.discretization_files);
    report.value("chain_gauge", chain_gauge_name(configuration.chain.gauge));
    if (configuration.chain.gauge == ChainGauge::nambu)
      report.value("nambu_tolerance", configuration.chain.nambu_tolerance);
  }
  report.write(std::cerr);
}

void report_star_configuration(const Configuration &configuration, const CommandLineOptions &command_line,
                               NRG::Tools::ConfigurationReport &report) {
  const auto &star = configuration.star;
  report.value("channels", configuration.gamma.channels);
  report.value("dos_prefix", configuration.gamma.prefix);
  report.value("Lambda", static_cast<double>(star.Lambda));
  if (command_line.Nz)
    report.resolved("z", "i/" + std::to_string(*command_line.Nz) + " for i=1.." + std::to_string(*command_line.Nz),
                    "--Nz; the z of the parameter file is ignored");
  else
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
  report.value("split_blocks", star.split_blocks ? std::string("true") : std::string("false"));
  report.value("allowed_error", star.allowed_error);
  report.value("hermiticity_tolerance", configuration.gamma.hermiticity_tolerance);
}

// One value of z and the directory its files go to: the working directory for a single z, and i/ for z = i/Nz.
struct Target {
  // The z to discretize for. In the chain stage of a single-z run it is empty, since the star records its own z.
  std::optional<double> z;
  std::filesystem::path directory;

  [[nodiscard]] auto file(const char *name) const { return (directory / name).string(); }
};

std::vector<Target> targets(const CommandLineOptions &command_line, const std::optional<double> single_z) {
  if (!command_line.Nz) return {Target{single_z, {}}};
  std::vector<Target> list;
  for (unsigned int i = 1; i <= *command_line.Nz; i++)
    list.push_back(Target{static_cast<double>(i) / *command_line.Nz, std::to_string(i)});
  return list;
}

auto seconds_since(const std::chrono::steady_clock::time_point start) {
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}

// The z-independent part of the star stage is done once; the stars for the individual values of z share it.
template<typename S> void run_star(const Configuration &configuration, const std::vector<Target> &list) {
  const auto setup_start = std::chrono::steady_clock::now();
  const auto input       = load_gamma<S>(configuration.gamma);
  StarDiscretizer<S> discretizer(input, configuration.star);
  std::cout << "# star setup: " << seconds_since(setup_start) << " s" << std::endl;
  report_blocks(discretizer, configuration.star.split_blocks, std::cout);

  for (const auto &target : list) {
    const auto start = std::chrono::steady_clock::now();
    if (!target.directory.empty()) std::cout << "# --- z=" << *target.z << " in " << target.directory.string() << "/" << std::endl;
    const auto star = discretizer.star(*target.z);
    report_star(star, std::cout);
    if (!target.directory.empty()) std::filesystem::create_directories(target.directory);
    const auto filename = target.file(star_default_filename);
    save_star(star, filename);
    std::cout << "# star written to " << filename << std::endl;
    std::cout << "# star z=" << star.z << ": " << seconds_since(start) << " s" << std::endl;
  }
}

// The star is self-contained: the chain is built with the Lambda, z and bandrescale it records. A parameter file that
// sets Lambda or bandrescale to something else was meant for a different star, so that is an error rather than a
// silent choice. The z must be the one asked for: i/Nz with --Nz, otherwise the z of the parameter file if it gives one.
template<typename S>
void check_star_against_parameters(const Star<S> &star, const Params &P, const Target &target,
                                   const std::string &filename) {
  const auto check = [&](const std::string &name, const double recorded, const double requested, const char *source) {
    if (std::abs(requested - recorded) <= 1e-12 * std::max(1.0, std::abs(recorded))) return;
    std::ostringstream message;
    message << std::setprecision(17) << "The star in " << filename << " was built with " << name << "=" << recorded
            << ", but " << source << " gives " << name << "=" << requested << ".";
    throw std::invalid_argument(message.str());
  };
  if (P.contains("Lambda")) check("Lambda", star.Lambda, P.P("Lambda", star.Lambda), "the parameter file");
  if (P.contains("bandrescale"))
    check("bandrescale", star.bandrescale, P.P("bandrescale", star.bandrescale), "the parameter file");
  if (target.z)
    check("z", star.z, *target.z, "--Nz");
  else if (P.contains("z"))
    check("z", star.z, P.P("z", star.z), "the parameter file");
}

// The chain is built from star.dat also in the default mode, right after the star stage has written it, so that the
// default mode and 's' followed by 'l' produce the same chain by construction.
template<typename S0> void run_chain(const Configuration &configuration, const Params &P, const Target &target) {
  const auto start     = std::chrono::steady_clock::now();
  const auto star_file = target.file(star_default_filename);
  if (!target.directory.empty()) std::cout << "# --- z=" << *target.z << " in " << target.directory.string() << "/" << std::endl;
  const auto star = load_star<S0>(star_file);
  check_star_against_parameters(star, P, target, star_file);
  const auto chain_file = target.file(chain_default_filename);
  const auto header = chain_header(star);
  const auto chain  = make_chain(star, configuration.chain);
  report_chain(chain, header, configuration.chain.sensitivity_tolerance, std::cout);
  save_chain(chain, header, chain_file);
  if (configuration.discretization_files) {
    save_chain_matrix_files(chain, header, target.directory);
    std::cout << "# matrix files written to " << (target.directory.empty() ? "." : target.directory.string())
              << std::endl;
  }
  std::cout << "# chain written to " << chain_file << std::endl;
  std::cout << "# chain z=" << star.z << ": " << seconds_since(start) << " s" << std::endl;
}

// Wall-clock time of one stage, for comparing the cost of the two.
template<typename F> void timed(const char *stage, F &&f) {
  const auto start = std::chrono::steady_clock::now();
  f();
  std::cout << "# " << stage << " stage: " << seconds_since(start) << " s" << std::endl;
}

void run(const CommandLineOptions &command_line) {
  const Params P(command_line.param_filename);
  const auto configuration =
    read_configuration(P, builds_star(command_line.mode), builds_chain(command_line.mode), command_line.cquad);
  report_configuration(configuration, command_line);

  if (builds_star(command_line.mode))
    timed("star", [&] {
      const auto list = targets(command_line, configuration.star.z);
      if (gamma_is_complex(configuration.gamma))
        run_star<std::complex<double>>(configuration, list);
      else
        run_star<double>(configuration, list);
    });
  if (builds_chain(command_line.mode))
    timed("chain", [&] {
      for (const auto &target : targets(command_line, std::nullopt)) {
        if (star_is_complex(target.file(star_default_filename)))
          run_chain<std::complex<double>>(configuration, P, target);
        else
          run_chain<double>(configuration, P, target);
      }
    });
}

} // namespace

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "mixchain")) return EXIT_SUCCESS;
  try {
    // Wall-clock time, like the stage timings, which then add up to it; the CPU time is shown alongside, since the
    // difference is time spent waiting, typically on the filesystem.
    const auto wall_start     = std::chrono::steady_clock::now();
    const clock_t cpu_start   = clock();
    about();
    run(cmd_line(argc, argv));
    const auto cpu_seconds = double(clock() - cpu_start) / CLOCKS_PER_SEC;
    std::cout << "# Elapsed " << seconds_since(wall_start) << " s (CPU " << cpu_seconds << " s)" << std::endl;
  } catch (const std::exception &e) {
    std::cerr << "mixchain: error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
}
