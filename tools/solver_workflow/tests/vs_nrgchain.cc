// discretize_hybridization compared with adapt --integral followed by nrgchain, for one channel; see vs_nrgchain_test.
//
// Usage: workflow_vs_nrgchain <case>=<tolerance on the FSOL grid>:<tolerance off it>...
// Each case directory holds param, Gamma_11-re.dat and nz<Nz>/<i>/{theta,xi,zeta}.dat from nrgchain at z = i/Nz.
// A z is on the FSOL grid when z/outputstep is an integer: then every x = j + z that nrgchain reads from FSOL is a
// node of the table. The comparison is relative, as in test/tools/mixchain/compare_nrgchain.sh:
//   theta:  |V^2 - theta| <= tol * theta
//   xi_n:   |T_n - xi_n|  <= tol * xi_n
//   zeta_n: |E_n - zeta_n| <= tol * max(|zeta_n|, xi_n)

#include <nrgljubljana/solver_chain.hpp>
#include <nrgljubljana/solver_preparation.hpp>

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

std::string read(const fs::path &file) {
  std::ifstream in(file);
  if (!in) throw std::runtime_error("Cannot read " + file.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

// The columns of a numeric file, comment lines skipped
std::vector<std::vector<double>> rows(const fs::path &file) {
  std::vector<std::vector<double>> result;
  std::istringstream lines(read(file));
  for (std::string line; std::getline(lines, line);) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream fields(line);
    std::vector<double> row;
    for (double x; fields >> x;) row.push_back(x);
    if (!row.empty()) result.push_back(row);
  }
  return result;
}

std::vector<double> column(const fs::path &file) {
  std::vector<double> result;
  for (const auto &row : rows(file)) result.push_back(row.at(0));
  return result;
}

struct Worst {
  double value = 0;
  std::string where;
  void update(const double difference, const std::string &at) {
    if (difference > value) value = difference, where = at;
  }
};

bool compare_case(const fs::path &directory, const double on_grid, const double off_grid) {
  const auto param = NRG::Workflow::read_parameter_sections(read(directory / "param")).at("param");
  NRG::Workflow::HybridizationFunction gamma;
  gamma.channels = 1;
  for (const auto &row : rows(directory / "Gamma_11-re.dat")) {
    gamma.omega.push_back(row.at(0));
    gamma.real.push_back(row.at(1));
  }
  const auto nmax = static_cast<unsigned>(std::stoul(param.at("Nmax")));

  const auto step = param.contains("outputstep") ? std::stod(param.at("outputstep")) : 0.0;
  const auto on   = [&](const double z) { return step == 0.0 || std::floor(z / step) == z / step; };
  Worst theta[2], xi[2], zeta[2]; // on and off the FSOL grid
  std::vector<unsigned> nzs;
  for (const auto &entry : fs::directory_iterator(directory))
    if (const auto name = entry.path().filename().string(); entry.is_directory() && name.rfind("nz", 0) == 0)
      nzs.push_back(static_cast<unsigned>(std::stoul(name.substr(2))));
  std::sort(nzs.begin(), nzs.end());
  if (nzs.empty()) throw std::runtime_error(directory.string() + ": no nz<Nz> directories");
  for (const unsigned nz : nzs) {
    std::vector<double> twists;
    for (unsigned i = 1; i <= nz; i++) twists.push_back(static_cast<double>(i) / nz);
    const auto log    = directory / ("discretize_hybridization-nz" + std::to_string(nz) + ".log");
    fs::remove(log);
    const auto chains = NRG::Workflow::discretize_hybridization(gamma, param, twists, nmax, log.string());
    for (unsigned i = 1; i <= nz; i++) {
      const auto &chain = chains.at(i - 1);
      const auto dir    = directory / ("nz" + std::to_string(nz)) / std::to_string(i);
      const auto where  = "Nz=" + std::to_string(nz) + " z=" + std::to_string(chain.z);
      const auto t      = column(dir / "theta.dat").at(0);
      const auto x      = column(dir / "xi.dat");
      const auto e      = column(dir / "zeta.dat");
      if (x.size() != chain.T.size() || e.size() != chain.E.size())
        throw std::runtime_error(directory.filename().string() + " " + where + ": the number of sites differs");
      const auto g = on(chain.z) ? 0 : 1;
      const auto v = chain.V.real.at(0);
      theta[g].update(std::abs(v * v - t) / t, where);
      for (std::size_t n = 0; n < x.size(); n++) {
        const auto tn = chain.T[n].real.at(0), en = chain.E[n].real.at(0);
        xi[g].update(std::abs(tn - x[n]) / x[n], where + " site " + std::to_string(n));
        zeta[g].update(std::abs(en - e[n]) / std::max(std::abs(e[n]), x[n]), where + " site " + std::to_string(n));
      }
    }
  }
  bool passed = true;
  for (const int g : {0, 1}) {
    const auto tolerance = g == 0 ? on_grid : off_grid;
    const bool ok        = theta[g].value <= tolerance && xi[g].value <= tolerance && zeta[g].value <= tolerance;
    passed               = passed && ok;
    std::cout << directory.filename().string() << (g == 0 ? " on" : " off") << " the FSOL grid (tolerance " << tolerance
              << "): theta " << theta[g].value << " (" << theta[g].where << "), xi " << xi[g].value << " (" << xi[g].where
              << "), zeta " << zeta[g].value << " (" << zeta[g].where << ") " << (ok ? "passed" : "FAILED") << '\n';
  }
  return passed;
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2) throw std::invalid_argument("Usage: workflow_vs_nrgchain <case>=<tolerance>:<tolerance>...");
    bool passed = true;
    for (int i = 1; i < argc; i++) {
      const std::string argument = argv[i];
      const auto eq              = argument.find('=');
      const auto colon           = argument.find(':', eq);
      if (eq == std::string::npos || colon == std::string::npos)
        throw std::invalid_argument("Expected <case>=<tolerance>:<tolerance>: " + argument);
      passed = compare_case(fs::absolute(argument.substr(0, eq)), std::stod(argument.substr(eq + 1, colon - eq - 1)),
                            std::stod(argument.substr(colon + 1)))
            && passed;
    }
    if (!passed) throw std::runtime_error("not within tolerance");
    std::cout << "Comparison with nrgchain passed.\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Comparison with nrgchain failed: " << error.what() << '\n';
    return 1;
  }
}
