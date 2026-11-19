// The SIAM through the SIAM/QS template, from the same Gamma by both routes, through the workflow library only:
//   nrgchain: discretize (both branches), nrgchain_coefficients, instantiate, solve, in N1/ and N2/
//   mixchain: discretize_hybridization for both twists, then diagonal_coefficients, instantiate, solve, in M1/ and M2/
// For z = 1/2 and 1, which lie on the FSOL grid, the chains agree to rounding (solver_workflow_vs_nrgchain), so the
// coefficients must agree to 1e-12, measured as there, and the solver's outputs to 1e-8 relative. The floor of 1e-10 is the magnitude
// below which two values count as equal: quantities that vanish by symmetry are rounding on one side and exact zeros
// on the other. No numbers are stored.
//
// Usage: workflow_routes <work directory>, prepared by routes_test: N1, N2, M1, M2 with template/ and param (with its
// z), N1 and N2 with Delta.dat, and Gamma_11-re.dat (the same Gamma) in the work directory.

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

void check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &file) {
  std::ifstream in(file);
  check(bool(in), "Cannot read " + file.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

// The largest relative difference |a-b|/max(|a|,|b|) between the numbers of two texts, values below floor on both
// sides counting as equal; the other tokens must be equal
double difference(const std::string &a, const std::string &b, const double floor, const std::string &what) {
  std::istringstream sa(a), sb(b);
  std::string ta, tb;
  double worst = 0;
  while (sa >> ta) {
    check(bool(sb >> tb), what + ": fewer tokens in the mixchain route");
    char *ea = nullptr, *eb = nullptr;
    const double x = std::strtod(ta.c_str(), &ea), y = std::strtod(tb.c_str(), &eb);
    if (ea != ta.c_str() && *ea == '\0' && eb != tb.c_str() && *eb == '\0') {
      if (std::abs(x) < floor && std::abs(y) < floor) continue;
      worst = std::max(worst, std::abs(x - y) / std::max(std::abs(x), std::abs(y)));
    } else
      check(ta == tb, what + ": '" + ta + "' differs from '" + tb + "'");
  }
  check(!(sb >> tb), what + ": more tokens in the mixchain route");
  return worst;
}

// The largest difference of the coefficients, as solver_workflow_vs_nrgchain measures it: theta and xi relatively,
// zeta relative to the scale of its site, max(|zeta|, xi), since zeta vanishes at particle-hole symmetry
double coefficient_difference(const CoefficientChannel &a, const CoefficientChannel &b) {
  check(a.xi.size() == b.xi.size() && a.zeta.size() == b.zeta.size(), "coefficients: the number of sites differs");
  double worst = std::abs(a.theta - b.theta) / a.theta;
  for (std::size_t n = 0; n < a.xi.size(); n++) {
    worst = std::max(worst, std::abs(a.xi[n] - b.xi[n]) / a.xi[n]);
    worst = std::max(worst, std::abs(a.zeta[n] - b.zeta[n]) / std::max(std::abs(a.zeta[n]), a.xi[n]));
  }
  return worst;
}

void run(const fs::path &work) {
  const auto param = read_parameter_sections(read(work / "N1" / "param")).at("param");
  const auto nmax  = resolve_chain_length(param);
  validate_preparation_parameters(param, ChainSource::nrgchain);
  validate_preparation_parameters(param, ChainSource::mixchain);

  HybridizationFunction gamma;
  gamma.channels = 1;
  std::istringstream table(read(work / "Gamma_11-re.dat"));
  for (double w, g; table >> w >> g;) {
    gamma.omega.push_back(w);
    gamma.real.push_back(g);
  }
  const std::vector<double> twists{0.5, 1.0};
  const auto chains = discretize_hybridization(gamma, param, twists, nmax, (work / "discretize_hybridization.log").string());
  fs::create_directories(work / "scratch");
  const auto scratch = (work / "scratch").string();

  bool passed = true;
  for (std::size_t i = 0; i < twists.size(); i++) {
    const auto n = work / ("N" + std::to_string(i + 1)), m = work / ("M" + std::to_string(i + 1));

    discretize({n.string(), "param", Branch::positive}, "discretize-P.log");
    discretize({n.string(), "param", Branch::negative}, "discretize-N.log");
    const auto from_nrgchain = nrgchain_coefficients({n.string(), "param", "template", nmax}, "nrgchain.log");
    instantiate({n.string(), "param", "template", nmax}, from_nrgchain, "instantiate.log");
    solve(n.string(), scratch, "solve.log");

    const auto from_mixchain = diagonal_coefficients(chains[i]);
    instantiate({m.string(), "param", "template", nmax}, from_mixchain, "instantiate.log");
    solve(m.string(), scratch, "solve.log");

    const auto where = "z=" + std::to_string(twists[i]);
    const auto chain = coefficient_difference(from_nrgchain.channels.at(0), from_mixchain.channels.at(0));
    std::cout << where << ": coefficients " << chain << (chain <= 1e-12 ? "" : " FAILED") << '\n';
    passed = passed && chain <= 1e-12;
    for (const auto *file : {"td", "custom", "tdfdm", "customfdm", "spec_FDM_dens_A_d-A_d.dat"}) {
      check(fs::exists(n / "DONE") && fs::exists(m / "DONE"), where + ": no DONE");
      const auto d = difference(read(n / file), read(m / file), 1e-10, where + " " + file);
      std::cout << where << ": " << file << ' ' << d << (d <= 1e-8 ? "" : " FAILED") << '\n';
      passed = passed && d <= 1e-8;
    }
  }
  check(passed, "the routes differ beyond tolerance");
}

} // namespace

int main(int argc, char **argv) {
  try {
    check(argc == 2, "Usage: workflow_routes <work directory>");
    run(fs::absolute(argv[1]));
    std::cout << "Both routes agree.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Routes: " << error.what() << '\n';
    return 1;
  }
}
