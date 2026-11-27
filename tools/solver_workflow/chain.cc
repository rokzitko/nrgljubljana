#include <nrgljubljana/solver_chain.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <iomanip>
#include <optional>
#include <sstream>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <numeric>

#include "../mixchain/configuration.hpp"
#include "../mixchain/discretize.hpp"
#include "../mixchain/load.hpp"
#include "../instantiate/workflow.hpp"

namespace NRG::Workflow {
namespace detail {
void guarded(const std::string &diagnostic_file, const std::function<void()> &function);
void prepared(const std::string &directory, const std::string &diagnostic_file, const std::function<void()> &function);
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

// The element (i,j) of m, which must vanish exactly unless keep; name and site describe it in the error.
double element(const ChainMatrix &m, const unsigned i, const unsigned j, const bool keep, const std::string &name,
               const std::optional<unsigned> site) {
  const auto where = [&] {
    return name + (site ? "_" + std::to_string(*site) : std::string()) + "(" + std::to_string(i + 1) + "," + std::to_string(j + 1) + ")";
  };
  const auto index = static_cast<std::size_t>(i) * m.cols + j;
  if (!m.imaginary.empty() && m.imaginary.at(index) != 0.0) {
    std::ostringstream message;
    message << std::setprecision(17) << "The chain is not diagonal: Im " << where() << " = " << m.imaginary[index] << " is not zero.";
    throw std::invalid_argument(message.str());
  }
  const auto value = m.real.at(index);
  if (!keep && value != 0.0) {
    std::ostringstream message;
    message << std::setprecision(17) << "The chain is not diagonal: " << where() << " = " << value
            << " is not zero. Channels that Gamma couples, or a diagonal Gamma discretized with split_blocks=false, "
               "need a chain of their own.";
    throw std::invalid_argument(message.str());
  }
  return value;
}

std::array<double, 4> qs_block(const ChainMatrix &matrix, const std::string &name) {
  if (matrix.rows != 2 || matrix.cols != 2 || matrix.real.size() != 4 ||
      (!matrix.imaginary.empty() && matrix.imaginary.size() != 4))
    throw std::invalid_argument(name + " must be a real 2 by 2 Wilson block.");
  if (std::any_of(matrix.imaginary.begin(), matrix.imaginary.end(), [](double x) { return x != 0.0; }) ||
      std::any_of(matrix.real.begin(), matrix.real.end(), [](double x) { return !std::isfinite(x); }))
    throw std::invalid_argument(name + " must contain finite real Wilson coefficients.");
  const auto &a = matrix.real;
  const auto scale = std::max({std::abs(a[0]), std::abs(a[1]), std::abs(a[2]), std::abs(a[3])});
  if (std::abs(a[1] - a[2]) > 32 * std::numeric_limits<double>::epsilon() * scale)
    throw std::invalid_argument(name + " is not symmetric in the polar gauge.");
  const auto cross = std::midpoint(a[1], a[2]);
  return {a[0], cross, cross, a[3]};
}

} // namespace

ChainCoefficients diagonal_coefficients(const WilsonChain &chain) {
  if (chain.gauge != "polar") throw std::invalid_argument("diagonal_coefficients requires the polar gauge, not " + chain.gauge + ".");
  const auto n = chain.channels;
  const auto square = [n](const ChainMatrix &m, const std::string &name) {
    if (m.rows != n || m.cols != n || m.real.size() != static_cast<std::size_t>(n) * n
        || !(m.imaginary.empty() || m.imaginary.size() == m.real.size()))
      throw std::invalid_argument("diagonal_coefficients: " + name + " is not " + std::to_string(n) + "x" + std::to_string(n) + ".");
  };
  if (chain.E.size() != chain.nmax + 1 || chain.T.size() != chain.nmax + 1)
    throw std::invalid_argument("diagonal_coefficients: expected " + std::to_string(chain.nmax + 1) + " sites.");
  square(chain.V, "V");
  for (unsigned s = 0; s <= chain.nmax; s++) {
    square(chain.E[s], "E_" + std::to_string(s));
    square(chain.T[s], "T_" + std::to_string(s));
  }

  ChainCoefficients result;
  result.z           = chain.z;
  result.Lambda      = chain.Lambda;
  result.bandrescale = chain.bandrescale;
  result.channels.resize(n);
  for (unsigned i = 0; i < n; i++)
    for (unsigned j = 0; j < n; j++) {
      const auto v = element(chain.V, i, j, i == j, "V", std::nullopt);
      if (i == j) result.channels[i].theta = v * v;
    }
  for (unsigned s = 0; s <= chain.nmax; s++)
    for (unsigned i = 0; i < n; i++)
      for (unsigned j = 0; j < n; j++) {
        const auto e = element(chain.E[s], i, j, i == j, "E", s);
        const auto t = element(chain.T[s], i, j, i == j, "T", s);
        if (i == j) {
          result.channels[i].zeta.push_back(e);
          result.channels[i].xi.push_back(t);
        }
      }
  return result;
}

std::vector<WilsonChain> discretize_hybridization(const HybridizationFunction &gamma,
                                                  const std::map<std::string, std::string> &param_section,
                                                  const std::vector<double> &twists, const unsigned nmax,
                                                  const std::string &diagnostic_file) {
  std::vector<WilsonChain> result;
  detail::guarded(diagnostic_file, [&] { result = compute(gamma, param_section, twists, nmax); });
  return result;
}

ChainCoefficients nrgchain_coefficients(const InitializerOptions &options, const std::string &diagnostic_file) {
  ChainCoefficients result;
  detail::prepared(options.calculation_directory, diagnostic_file, [&] {
    const auto computed = detail::nrgchain_coefficients(options.parameter_file, options.nmax);
    result.z            = computed.z;
    result.Lambda       = computed.Lambda;
    result.bandrescale  = computed.bandrescale;
    for (const auto &channel : computed.data.channels) result.channels.push_back({channel.theta, channel.xi, channel.zeta});
  });
  return result;
}

InitializerResult instantiate(const InitializerOptions &options, const ChainCoefficients &coefficients,
                              const std::string &diagnostic_file) {
  NRG::Tools::NrgChain::WilsonData wilson;
  for (const auto &channel : coefficients.channels) wilson.channels.push_back({channel.theta, channel.xi, channel.zeta});
  InitializerResult result{};
  detail::prepared(options.calculation_directory, diagnostic_file, [&] {
    const auto energy = detail::instantiate(options.parameter_file, options.template_directory, options.nmax, wilson,
                                            coefficients.z, coefficients.Lambda, coefficients.bandrescale);
    result = InitializerResult{std::filesystem::absolute("data").string(), energy};
  });
  return result;
}

InitializerResult instantiate(const InitializerOptions &options, const WilsonChain &chain,
                              const std::string &diagnostic_file) {
  if (chain.channels != 2 || chain.gauge != "polar")
    throw std::invalid_argument("Matrix instantiate requires a two-channel Wilson chain in the polar gauge.");
  if (options.nmax == 0 || options.nmax > 998 || chain.nmax < options.nmax ||
      chain.E.size() != static_cast<std::size_t>(chain.nmax) + 1 || chain.T.size() != chain.E.size())
    throw std::invalid_argument("Matrix Wilson coefficients do not cover the requested sites 0..nmax.");

  detail::QSMixing mixing{qs_block(chain.V, "V"), {}, {}};
  NRG::Tools::NrgChain::WilsonData wilson;
  wilson.channels.resize(2);
  for (const unsigned channel : {0u, 1u}) {
    const auto v = mixing.coupling[3 * channel];
    if (v < 0.0) throw std::invalid_argument("The polar Wilson coupling must have nonnegative diagonal elements.");
    wilson.channels[channel].theta = v * v;
  }
  for (unsigned site = 0; site <= options.nmax; ++site) {
    const auto e = qs_block(chain.E[site], "E_" + std::to_string(site));
    const auto t = qs_block(chain.T[site], "T_" + std::to_string(site));
    for (const unsigned channel : {0u, 1u}) {
      wilson.channels[channel].xi.push_back(t[3 * channel]);
      wilson.channels[channel].zeta.push_back(e[3 * channel]);
    }
    mixing.hopping.push_back(t[1]);
    mixing.onsite.push_back(e[1]);
  }

  InitializerResult result{};
  detail::prepared(options.calculation_directory, diagnostic_file, [&] {
    const auto energy = detail::instantiate(options.parameter_file, options.template_directory, options.nmax, wilson,
                                            chain.z, chain.Lambda, chain.bandrescale, &mixing);
    result = InitializerResult{std::filesystem::absolute("data").string(), energy};
  });
  return result;
}

} // namespace NRG::Workflow
