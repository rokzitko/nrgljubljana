#include <nrgljubljana/solver_chain.hpp>

#include <algorithm>
#include <complex>
#include <cstddef>
#include <functional>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

#include "../mixchain/configuration.hpp"
#include "../mixchain/discretize.hpp"
#include "../mixchain/load.hpp"

namespace NRG::Workflow {
namespace detail {
void guarded(const std::string &diagnostic_file, const std::function<void()> &function);
}
namespace {

namespace MC = NRG::MixChain;

// The element as chain.dat writes it: factor times the element, in the arithmetic of the chain.
template <typename S> ChainMatrix chain_matrix(const MC::Matrix<S> &m, const double factor) {
  ChainMatrix result;
  result.rows = static_cast<unsigned>(m.rows());
  result.cols = static_cast<unsigned>(m.cols());
  for (Eigen::Index i = 0; i < m.rows(); i++)
    for (Eigen::Index j = 0; j < m.cols(); j++) {
      const S value = MC::make_scalar<S>(factor, 0) * m(i, j);
      result.real.push_back(std::real(value));
      if constexpr (MC::is_complex_v<S>) result.imaginary.push_back(std::imag(value));
    }
  return result;
}

template <typename S> WilsonChain wilson_chain(const MC::Discretization<S> &d) {
  const auto &chain  = d.chain;
  const auto header  = MC::chain_header(d.star);
  const auto &diag   = chain.diagnostics;
  WilsonChain result;
  result.z           = header.z;
  result.Lambda      = header.Lambda;
  result.bandrescale = header.bandrescale;
  result.channels    = static_cast<unsigned>(chain.channels);
  result.nmax        = chain.Nmax;
  result.gauge       = MC::chain_gauge_name(chain.gauge);
  result.V           = chain_matrix(chain.V, 1.0);
  for (const auto &m : chain.E) result.E.push_back(chain_matrix(m, header.bandrescale));
  for (const auto &m : chain.T) result.T.push_back(chain_matrix(m, header.bandrescale));
  for (const auto &block : chain.blocks) result.blocks.emplace_back(block.begin(), block.end());
  result.levels                 = diag.levels;
  result.coupled_levels         = diag.coupled_levels;
  result.theta_rank             = diag.theta_rank;
  result.min_rank               = diag.min_rank;
  result.rank_drop_site         = diag.rank_drop_site;
  result.continued_from_site    = MC::first_continued_site(chain, header.untabulated_to - header.untabulated_from);
  result.theta_condition        = diag.theta_condition;
  result.min_residual_condition = diag.min_residual_condition;
  result.max_star_sensitivity   = diag.max_star_sensitivity;
  result.sensitive_from_site    = diag.sensitive_from_site;
  for (const auto &star_diagnostics : d.star.diagnostics)
    result.max_interval_deviation = std::max(result.max_interval_deviation, star_diagnostics.max_interval_deviation);
  return result;
}

template <typename S>
std::vector<WilsonChain> run(const MC::GammaTable &table, const MC::Configuration &configuration, const std::vector<double> &twists) {
  const auto input = MC::gamma_input<S>(table, configuration.gamma, std::cout);
  std::vector<WilsonChain> result;
  for (const auto &d : MC::discretize(input, configuration, twists, std::cout)) result.push_back(wilson_chain(d));
  return result;
}

std::vector<WilsonChain> compute(const HybridizationFunction &gamma, const std::map<std::string, std::string> &param_section,
                                 const std::vector<double> &twists, const unsigned nmax) {
  const auto n = static_cast<std::size_t>(gamma.channels);
  if (n == 0) throw std::invalid_argument("Gamma must have at least one channel.");
  if (gamma.real.size() != gamma.omega.size() * n * n)
    throw std::invalid_argument("Gamma: real must hold channels^2 values per frequency.");
  if (!gamma.imaginary.empty() && gamma.imaginary.size() != gamma.real.size())
    throw std::invalid_argument("Gamma: imaginary must be empty or of the size of real.");
  if (twists.empty()) throw std::invalid_argument("At least one twist is required.");
  if (nmax == 0) throw std::invalid_argument("nmax must be greater than 0.");

  auto entries = param_section;
  for (const auto *key : {"channels", "dos_prefix", "z"}) entries.erase(key);
  entries["Nmax"] = std::to_string(nmax);
  if (const auto it = entries.find("discretization"); it != entries.end() && (it->second.empty() || it->second.front() != 'Z'))
    throw std::invalid_argument("discretization=" + it->second + ": the chain is built in the Z scheme only.");
  const auto configuration = MC::read_configuration(entries, true, true);

  MC::GammaTable table;
  table.channels = static_cast<int>(n);
  table.omega    = gamma.omega;
  table.complex  = !gamma.imaginary.empty();
  table.gamma.resize(gamma.real.size());
  for (std::size_t k = 0; k < gamma.real.size(); k++)
    table.gamma[k] = {gamma.real[k], table.complex ? gamma.imaginary[k] : 0.0};

  return table.complex ? run<std::complex<double>>(table, configuration, twists) : run<double>(table, configuration, twists);
}

} // namespace

std::vector<WilsonChain> discretize_hybridization(const HybridizationFunction &gamma,
                                                  const std::map<std::string, std::string> &param_section,
                                                  const std::vector<double> &twists, const unsigned nmax,
                                                  const std::string &diagnostic_file) {
  std::vector<WilsonChain> result;
  detail::guarded(diagnostic_file, [&] { result = compute(gamma, param_section, twists, nmax); });
  return result;
}

} // namespace NRG::Workflow
