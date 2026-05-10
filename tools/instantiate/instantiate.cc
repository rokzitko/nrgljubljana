// Initial data instantiation tool
// First production slice: in-process Wilson-chain generation.

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "matrix_evaluator.hpp"
#include "nrgchain.hpp"

namespace {

struct Options {
  std::string param_filename = "param";
  bool wilson_only = false;
};

void usage(std::ostream &out = std::cout) {
  out << "Usage: instantiate [--wilson-only] [--param FILE]\n";
}

Options parse_options(const int argc, char *argv[]) {
  Options options;
  for (int i = 1; i < argc; ++i) {
    const std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      std::exit(EXIT_SUCCESS);
    }
    if (arg == "--wilson-only") {
      options.wilson_only = true;
      continue;
    }
    if (arg == "--param") {
      if (i + 1 >= argc) throw std::invalid_argument("--param requires a filename.");
      options.param_filename = argv[++i];
      continue;
    }
    throw std::invalid_argument("Unknown argument: " + arg);
  }
  return options;
}

void write_values(const std::string &filename, const std::vector<double> &values) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(16);
  for (const auto value : values) out << value << '\n';
}

void write_scalar(const std::string &filename, const double value) {
  std::ofstream out(filename);
  if (!out) throw std::runtime_error("Can't open " + filename + " for writing.");
  out << std::setprecision(16) << value << '\n';
}

void write_wilson_channel(const NRG::Tools::NrgChain::WilsonChannel &channel, const size_t index) {
  const auto suffix = std::to_string(index) + ".dat";
  write_values("xi" + suffix, channel.xi);
  write_values("zeta" + suffix, channel.zeta);
  write_scalar("theta" + suffix, channel.theta);
}

void run_wilson_only(const Options &options) {
  const auto wilson = NRG::Tools::NrgChain::calculate_from_file(options.param_filename);
  if (wilson.channels.size() != 1)
    throw std::runtime_error("Only single-channel Wilson generation is supported in this instantiate slice.");
  write_wilson_channel(wilson.channels.front(), 1);
}

} // namespace

int main(int argc, char *argv[]) {
  try {
    const auto options = parse_options(argc, argv);
    if (!options.wilson_only)
      throw std::runtime_error("Full data instantiation is not implemented yet; use --wilson-only for this slice.");
    run_wilson_only(options);
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << "instantiate: error: " << e.what() << '\n';
    return EXIT_FAILURE;
  }
}
