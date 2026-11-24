// instantiate with supplied Wilson coefficients (solver_chain.hpp), against the instantiate that computes them.
//
// Usage: workflow_instantiate_chain <work directory>
// The work directory holds A, B, C and D, each with template/ and param (Nmax=6); A and B have Delta.dat, C and D
// the same parameters with band=flat. The test:
//   - A: discretize and instantiate (nrgchain inside). B: the overload with nrgchain_coefficients of A. data of A and
//     B must be byte-identical, and nrgchain_coefficients leaves no files in A.
//   - C: instantiate with band=flat. D: discretize_hybridization of the same flat band, diagonal_coefficients and the
//     overload. The numbers of data agree to rounding (see solver_workflow_vs_nrgchain).
//   - Errors: z, Lambda or bandrescale other than the parameter file's, too few sites, two channels.

#include <nrgljubljana/solver_chain.hpp>
#include <nrgljubljana/solver_preparation.hpp>
#include <nrgljubljana/solver_workflow.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace NRG::Workflow;
unsigned assertions = 0;

void check(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &file) {
  std::ifstream in(file, std::ios::binary);
  check(bool(in), "Cannot read " + file.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

template <typename Function> void rejects(Function function, const std::string &detail) {
  try {
    function();
  } catch (const std::exception &error) {
    check(std::string(error.what()).find(detail) != std::string::npos,
          "Expected an error containing '" + detail + "', got '" + error.what() + "'");
    return;
  }
  throw std::runtime_error("Expected an error containing: " + detail);
}

std::map<std::string, std::string> param_section(const fs::path &directory) {
  return read_parameter_sections(read(directory / "param")).at("param");
}

double number(const std::map<std::string, std::string> &param, const std::string &key, const double default_value) {
  const auto it = param.find(key);
  return it == param.end() ? default_value : std::stod(it->second);
}

// The names in a directory
std::vector<std::string> listing(const fs::path &directory) {
  std::vector<std::string> names;
  for (const auto &entry : fs::directory_iterator(directory)) names.push_back(entry.path().filename().string());
  std::sort(names.begin(), names.end());
  return names;
}

// The numbers of two data files agree to tolerance (relative, or absolute below 1); everything else is equal
void same_numbers(const std::string &a, const std::string &b, const double tolerance) {
  std::istringstream sa(a), sb(b);
  std::string ta, tb;
  std::size_t tokens = 0;
  double worst       = 0;
  while (sa >> ta) {
    check(bool(sb >> tb), "data: fewer tokens in the second file");
    ++tokens;
    char *ea = nullptr, *eb = nullptr;
    const double x = std::strtod(ta.c_str(), &ea), y = std::strtod(tb.c_str(), &eb);
    if (*ea == '\0' && *eb == '\0' && ea != ta.c_str()) {
      worst = std::max(worst, std::abs(x - y) / std::max(1.0, std::abs(x)));
      check(std::abs(x - y) <= tolerance * std::max(1.0, std::abs(x)), "data: " + ta + " differs from " + tb);
    } else
      check(ta == tb, "data: '" + ta + "' differs from '" + tb + "'");
  }
  check(!(sb >> tb), "data: more tokens in the second file");
  std::cout << "  flat band: " << tokens << " tokens, largest difference " << worst << '\n';
}

void run(const fs::path &work) {
  const unsigned nmax = 6;
  const auto a = work / "A", b = work / "B", c = work / "C", d = work / "D";

  // A and B: byte-identical data
  discretize({a.string(), "param", Branch::positive});
  discretize({a.string(), "param", Branch::negative});
  const auto ra           = instantiate({a.string(), "param", "template", nmax});
  const auto before       = listing(a);
  const auto coefficients = nrgchain_coefficients({a.string(), "param", "template", nmax});
  check(listing(a) == before, "nrgchain_coefficients left files behind");
  check(coefficients.z == 0.7 && coefficients.Lambda == 2 && coefficients.bandrescale == 1, "z, Lambda, bandrescale");
  check(coefficients.channels.size() == 1 && coefficients.channels[0].xi.size() == nmax + 1, "Coefficient channels and sites");
  const auto rb           = instantiate({b.string(), "param", "template", nmax}, coefficients, "instantiate.log");
  check(fs::path(rb.data_file).is_absolute(), "data_file of B is not absolute: " + rb.data_file);
  check(fs::equivalent(rb.data_file, b / "data"),
        "data_file of B: expected " + (b / "data").string() + ", got " + rb.data_file);
  check(rb.ground_energy == ra.ground_energy, "Ground energy of A and B");
  check(read(a / "data") == read(b / "data"), "data of A and B differ");
  check(read(b / "instantiate.log").find("# Wilson coefficients supplied: channels=1") != std::string::npos, "Log of B");
  std::cout << "  nrgchain_coefficients and the overload: data byte-identical\n";

  // A longer chain is cut to nmax
  auto longer = coefficients;
  for (auto &channel : longer.channels)
    for (int k = 0; k < 4; k++) {
      channel.xi.push_back(1.0);
      channel.zeta.push_back(0.5);
    }
  fs::remove(b / "data");
  instantiate({b.string(), "param", "template", nmax}, longer);
  check(read(a / "data") == read(b / "data"), "data of a longer chain differs");

  // C and D: the mixchain route on a flat band
  instantiate({c.string(), "param", "template", nmax});
  const auto param = param_section(d);
  HybridizationFunction gamma;
  gamma.channels = 1;
  for (int k = -100; k <= 100; k++) {
    gamma.omega.push_back(k / 100.0);
    gamma.real.push_back(0.5);
  }
  const auto chain = discretize_hybridization(gamma, param, {number(param, "z", 1.0)}, nmax).at(0);
  instantiate({d.string(), "param", "template", nmax}, diagonal_coefficients(chain));
  same_numbers(read(c / "data"), read(d / "data"), 1e-12);

  // Errors
  const auto refused = [&](ChainCoefficients bad, const std::string &detail) {
    rejects([&] { instantiate({b.string(), "param", "template", nmax}, bad); }, detail);
  };
  auto bad = coefficients;
  bad.z    = 0.5;
  refused(bad, "computed for z=0.5");
  bad        = coefficients;
  bad.Lambda = 3;
  refused(bad, "computed for Lambda=3");
  bad             = coefficients;
  bad.bandrescale = 2;
  refused(bad, "computed for bandrescale=2");
  bad = coefficients;
  bad.channels[0].xi.resize(nmax);
  refused(bad, "do not cover the sites 0..6 required");
  bad = coefficients;
  bad.channels.push_back(bad.channels[0]);
  refused(bad, "channel count");
  rejects([&] { instantiate({(work / "missing").string(), "param", "template", nmax}, coefficients); }, "calculation_directory");
  std::cout << "  errors: ok\n";
}

} // namespace

int main(int argc, char **argv) {
  try {
    check(argc == 2, "Usage: workflow_instantiate_chain <work directory>");
    run(fs::absolute(argv[1]));
    std::cout << "instantiate with supplied coefficients: tests passed (" << assertions << " checks).\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "instantiate with supplied coefficients: tests failed: " << error.what() << '\n';
    return 1;
  }
}
