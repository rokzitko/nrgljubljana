// discretize_hybridization through the installed headers only.
//
// Usage: workflow_chain <case directory>...
// The first case has one channel. Each case directory holds param-mixchain, the Gamma files, and 1/chain.dat and 2/chain.dat from mixchain --Nz 2.
// The chains of z = 1/2 and 1 must equal those files exactly (they are written with 18 digits). Then come the
// errors, the state restoration and the log routing.

#include <nrgljubljana/solver_chain.hpp>
#include <nrgljubljana/solver_preparation.hpp>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using NRG::Workflow::ChainMatrix;
using NRG::Workflow::HybridizationFunction;
using NRG::Workflow::WilsonChain;
unsigned assertions = 0;

void check(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &file) {
  std::ifstream in(file);
  check(bool(in), "Cannot read " + file.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

double number(const std::string &text) {
  char *end         = nullptr;
  const auto result = std::strtod(text.c_str(), &end);
  check(end != text.c_str() && *end == '\0', "Not a number: " + text);
  return result;
}

std::vector<std::pair<double, double>> pairs(const fs::path &file) {
  std::vector<std::pair<double, double>> result;
  std::istringstream lines(read(file));
  for (std::string line; std::getline(lines, line);) {
    if (line.empty() || line[0] == '#') continue;
    std::istringstream fields(line);
    std::string x, y;
    fields >> x >> y;
    result.emplace_back(number(x), number(y));
  }
  return result;
}

// Gamma from the files Gamma_ij-re.dat and, if present, Gamma_ij-im.dat (absent ones are zero)
HybridizationFunction read_gamma(const fs::path &directory, const unsigned n) {
  HybridizationFunction gamma;
  gamma.channels = n;
  const auto name = [&](unsigned i, unsigned j, const char *part) {
    return directory / ("Gamma_" + std::to_string(i + 1) + std::to_string(j + 1) + "-" + part + ".dat");
  };
  const bool complex = n > 1 && fs::exists(name(0, 1, "im"));
  const auto grid    = pairs(name(0, 0, "re"));
  for (const auto &p : grid) gamma.omega.push_back(p.first);
  gamma.real.assign(grid.size() * n * n, 0.0);
  if (complex) gamma.imaginary.assign(grid.size() * n * n, 0.0);
  for (unsigned i = 0; i < n; i++)
    for (unsigned j = 0; j < n; j++) {
      const auto re = pairs(name(i, j, "re"));
      for (std::size_t k = 0; k < grid.size(); k++) gamma.real[(k * n + i) * n + j] = re.at(k).second;
      if (complex && fs::exists(name(i, j, "im"))) {
        const auto im = pairs(name(i, j, "im"));
        for (std::size_t k = 0; k < grid.size(); k++) gamma.imaginary[(k * n + i) * n + j] = im.at(k).second;
      }
    }
  return gamma;
}

// chain.dat: the key=value pairs of the comment lines, the blocks line, and the elements
struct ChainFile {
  std::map<std::string, std::string> keys;
  std::string blocks;
  std::map<std::string, std::vector<double>> elements; // "E 3 1 2" -> value or Re, Im
};

ChainFile read_chain(const fs::path &file) {
  ChainFile result;
  std::istringstream lines(read(file));
  for (std::string line; std::getline(lines, line);) {
    std::istringstream fields(line);
    std::string token;
    if (line.rfind("# blocks= ", 0) == 0) {
      result.blocks = line.substr(10);
    } else if (line[0] == '#') {
      while (fields >> token)
        if (const auto eq = token.find('='); eq != std::string::npos) result.keys[token.substr(0, eq)] = token.substr(eq + 1);
    } else {
      std::string name, n, i, j;
      fields >> name >> n >> i >> j;
      auto &values = result.elements[name + " " + n + " " + i + " " + j];
      while (fields >> token) values.push_back(number(token));
    }
  }
  return result;
}

std::string site(const std::optional<unsigned> &s) { return s ? std::to_string(*s) : "none"; }

void compare_matrix(const ChainFile &file, const std::string &name, const unsigned n, const ChainMatrix &m, const bool complex,
                    const std::string &where) {
  for (unsigned i = 0; i < m.rows; i++)
    for (unsigned j = 0; j < m.cols; j++) {
      const auto key = name + " " + std::to_string(n) + " " + std::to_string(i + 1) + " " + std::to_string(j + 1);
      const auto it  = file.elements.find(key);
      check(it != file.elements.end(), where + ": " + key + " missing in chain.dat");
      const auto index = i * m.cols + j;
      std::ostringstream message;
      message << std::setprecision(18) << where << ": " << key << " differs: " << m.real[index];
      check(it->second.size() == (complex ? 2u : 1u) && it->second[0] == m.real[index], message.str());
      if (complex) check(it->second[1] == m.imaginary[index], where + ": Im " + key + " differs");
    }
}

void compare(const WilsonChain &chain, const fs::path &file) {
  const auto where   = file.string();
  const auto f       = read_chain(file);
  const bool complex = f.keys.at("complex") == "1";
  check(!chain.V.real.empty() && chain.V.imaginary.empty() != complex, where + ": arithmetic");
  check(number(f.keys.at("z")) == chain.z && number(f.keys.at("Lambda")) == chain.Lambda
          && number(f.keys.at("bandrescale")) == chain.bandrescale,
        where + ": z, Lambda or bandrescale");
  check(f.keys.at("channels") == std::to_string(chain.channels) && f.keys.at("Nmax") == std::to_string(chain.nmax)
          && f.keys.at("gauge") == chain.gauge,
        where + ": channels, Nmax or gauge");
  check(f.keys.at("levels") == std::to_string(chain.levels) && f.keys.at("coupled_levels") == std::to_string(chain.coupled_levels)
          && f.keys.at("theta_rank") == std::to_string(chain.theta_rank) && f.keys.at("min_rank") == std::to_string(chain.min_rank),
        where + ": levels or ranks");
  check(f.keys.at("rank_drop_site") == site(chain.rank_drop_site) && f.keys.at("continued_from_site") == site(chain.continued_from_site)
          && f.keys.at("sensitive_from_site") == site(chain.sensitive_from_site),
        where + ": diagnostic sites");
  check(number(f.keys.at("theta_condition")) == chain.theta_condition
          && number(f.keys.at("min_residual_condition")) == chain.min_residual_condition
          && number(f.keys.at("max_star_sensitivity")) == chain.max_star_sensitivity,
        where + ": diagnostic values");
  check(std::isfinite(chain.max_interval_deviation) && chain.max_interval_deviation >= 0, where + ": max_interval_deviation");
  std::string blocks;
  for (const auto &block : chain.blocks) {
    blocks += blocks.empty() ? "{" : " {";
    for (std::size_t k = 0; k < block.size(); k++) blocks += (k ? "," : "") + std::to_string(block[k] + 1);
    blocks += "}";
  }
  check(f.blocks == (chain.blocks.size() > 1 ? blocks : ""), where + ": blocks " + blocks);

  check(chain.E.size() == chain.nmax + 1 && chain.T.size() == chain.nmax + 1, where + ": sites");
  std::size_t elements = chain.V.rows * chain.V.cols;
  compare_matrix(f, "V", 0, chain.V, complex, where);
  for (unsigned n = 0; n <= chain.nmax; n++) {
    compare_matrix(f, "E", n, chain.E[n], complex, where);
    compare_matrix(f, "T", n, chain.T[n], complex, where);
    elements += chain.E[n].rows * chain.E[n].cols + chain.T[n].rows * chain.T[n].cols;
  }
  check(f.elements.size() == elements, where + ": number of elements");
}

std::map<std::string, std::string> param_section(const fs::path &directory) {
  return NRG::Workflow::read_parameter_sections(read(directory / "param-mixchain")).at("param");
}

void equals_mixchain(const fs::path &directory) {
  const auto param = param_section(directory);
  const auto gamma = read_gamma(directory, param.contains("channels") ? static_cast<unsigned>(std::stoul(param.at("channels"))) : 1u);
  const auto nmax  = static_cast<unsigned>(std::stoul(param.at("Nmax")));
  const auto log   = directory / "discretize_hybridization.log";
  fs::remove(log);
  const auto chains = NRG::Workflow::discretize_hybridization(gamma, param, {0.5, 1.0}, nmax, log.string());
  check(chains.size() == 2, "Two twists");
  compare(chains[0], directory / "1" / "chain.dat");
  compare(chains[1], directory / "2" / "chain.dat");
  check(read(log).find("# Gamma: channels=" + std::to_string(gamma.channels)) != std::string::npos, "Log written");
  std::cout << "  " << directory.filename().string() << ": equal to mixchain\n";
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

void errors_and_state(const fs::path &directory) {
  const auto param = param_section(directory);
  const auto gamma = read_gamma(directory, 1);
  const auto call  = [&](const HybridizationFunction &g, std::map<std::string, std::string> p, std::vector<double> twists,
                        unsigned nmax) { return NRG::Workflow::discretize_hybridization(g, p, twists, nmax, "chain_errors.log"); };

  // Keys that Gamma or the call give are ignored; Nmax is replaced by nmax
  auto p          = param;
  p["channels"]   = "3";
  p["dos_prefix"] = "nothing";
  p["z"]          = "0.25";
  p["Nmax"]       = "4";
  p["discretization"] = "Z";
  fs::remove("chain_errors.log");
  const auto chains = call(gamma, p, {0.5}, 6);
  check(chains.size() == 1 && chains[0].z == 0.5 && chains[0].nmax == 6 && chains[0].channels == 1, "Ignored keys");

  p["discretization"] = "Y";
  rejects([&] { fs::remove("chain_errors.log"); call(gamma, p, {0.5}, 6); }, "Z scheme");
  rejects([&] { fs::remove("chain_errors.log"); call(gamma, param, {}, 6); }, "twist");
  rejects([&] { fs::remove("chain_errors.log"); call(gamma, param, {1.5}, 6); }, "z must be in (0,1]");
  rejects([&] { fs::remove("chain_errors.log"); call(gamma, param, {0.5}, 0); }, "nmax");
  auto bad = gamma;
  bad.real.pop_back();
  rejects([&] { fs::remove("chain_errors.log"); call(bad, param, {0.5}, 6); }, "channels^2");
  bad = gamma;
  bad.imaginary.assign(3, 0.0);
  rejects([&] { fs::remove("chain_errors.log"); call(bad, param, {0.5}, 6); }, "imaginary");
  bad = gamma;
  for (std::size_t k = 0; k < bad.real.size(); k += 2) bad.real[k] = -bad.real[k];
  rejects([&] { fs::remove("chain_errors.log"); call(bad, param, {0.5}, 6); }, "positive semidefinite");

  // An existing log is not overwritten
  rejects([&] { call(gamma, param, {0.5}, 6); }, "chain_errors.log");
  fs::remove("chain_errors.log");

  // The caller's state is restored, and nothing goes to the caller's streams with a log
  const auto cwd = fs::current_path();
  std::cout << std::setprecision(3);
  call(gamma, param, {0.5}, 6);
  check(fs::current_path() == cwd && std::cout.precision() == 3, "State restored");
  fs::remove("chain_errors.log");
  std::cout << std::setprecision(6);
  std::cout << "  errors and state: ok\n";
}

// Gamma on 201 nodes in [-1, 1], with the diagonal elements given and the off-diagonal ones all equal to offdiagonal
HybridizationFunction in_memory(const std::vector<double (*)(double)> &diagonal, const double offdiagonal = 0) {
  HybridizationFunction gamma;
  const auto n   = static_cast<unsigned>(diagonal.size());
  gamma.channels = n;
  for (int k = -100; k <= 100; k++) {
    const double w = k / 100.0;
    gamma.omega.push_back(w);
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < n; j++) gamma.real.push_back(i == j ? diagonal[i](w) : offdiagonal);
  }
  return gamma;
}

double flat(double) { return 0.5; }
double quarter(double) { return 0.25; }
double semicircle(double w) { return 1 - w * w > 0 ? 0.5 * std::sqrt(1 - w * w) : 0.0; }

void adapter() {
  const std::map<std::string, std::string> param{{"Lambda", "2"}};
  const auto chain_of = [&](const HybridizationFunction &gamma) {
    return NRG::Workflow::discretize_hybridization(gamma, param, {0.5}, 10).at(0);
  };

  // One channel: theta = V^2, xi = T, zeta = E
  const auto one = chain_of(in_memory({flat}));
  const auto c   = NRG::Workflow::diagonal_coefficients(one);
  check(c.channels.size() == 1 && c.z == one.z && c.Lambda == one.Lambda && c.bandrescale == one.bandrescale, "Coefficients header");
  check(c.channels[0].theta == one.V.real[0] * one.V.real[0], "theta = V^2");
  check(c.channels[0].xi.size() == 11 && c.channels[0].zeta.size() == 11, "Sites 0..nmax");
  for (unsigned n = 0; n <= 10; n++)
    check(c.channels[0].xi[n] == one.T[n].real[0] && c.channels[0].zeta[n] == one.E[n].real[0], "xi = T, zeta = E");

  // A diagonal 2x2 Gamma: two channels, each equal to a run of its own
  const auto two  = NRG::Workflow::diagonal_coefficients(chain_of(in_memory({flat, semicircle})));
  const auto a    = NRG::Workflow::diagonal_coefficients(chain_of(in_memory({flat})));
  const auto b    = NRG::Workflow::diagonal_coefficients(chain_of(in_memory({semicircle})));
  const auto same = [](const NRG::Workflow::CoefficientChannel &x, const NRG::Workflow::CoefficientChannel &y) {
    return x.theta == y.theta && x.xi == y.xi && x.zeta == y.zeta;
  };
  check(two.channels.size() == 2 && same(two.channels[0], a.channels[0]) && same(two.channels[1], b.channels[0]),
        "Diagonal Gamma: the channels of separate runs");

  // Refused, naming the element
  rejects([&] { NRG::Workflow::diagonal_coefficients(chain_of(in_memory({flat, quarter}, 0.1))); }, "not diagonal: V(1,2)");
  auto diagonal = chain_of(in_memory({flat, semicircle}));
  auto modified = diagonal;
  modified.E[3].real[1] = 1e-300;
  rejects([&] { NRG::Workflow::diagonal_coefficients(modified); }, "not diagonal: E_3(1,2) = ");
  modified = diagonal;
  modified.T[0].real[2] = -1e-20;
  rejects([&] { NRG::Workflow::diagonal_coefficients(modified); }, "T_0(2,1)");
  modified = one;
  modified.V.imaginary = {0.1};
  for (auto &m : modified.E) m.imaginary = {0.0};
  rejects([&] { NRG::Workflow::diagonal_coefficients(modified); }, "Im V(1,1)");
  modified = one;
  modified.gauge = "nambu";
  rejects([&] { NRG::Workflow::diagonal_coefficients(modified); }, "polar gauge");
  std::cout << "  diagonal_coefficients: ok\n";
}

} // namespace

int main(int argc, char **argv) {
  try {
    check(argc >= 2, "Usage: workflow_chain <case directory>...");
    for (int i = 1; i < argc; i++) equals_mixchain(fs::absolute(argv[i]));
    errors_and_state(fs::absolute(argv[1]));
    adapter();
    std::cout << "discretize_hybridization tests passed (" << assertions << " checks).\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "discretize_hybridization tests failed: " << error.what() << '\n';
    return 1;
  }
}
