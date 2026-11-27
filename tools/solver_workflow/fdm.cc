#include <nrgljubljana/solver_fdm.hpp>
#include <nrgljubljana/solver_output.hpp>

#include "../common/parallel.hpp"

#include <Eigen/LU>
#include <gsl/gsl_integration.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <memory>
#include <numbers>
#include <stdexcept>

namespace NRG::Workflow {
namespace {

using Complex = std::complex<double>;
using Matrix = Eigen::Matrix<Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
namespace fs = std::filesystem;

std::size_t square_size(unsigned n) {
  if (n == 0 || n > std::numeric_limits<std::size_t>::max() / n)
    throw std::invalid_argument("Invalid FDM orbital count.");
  return static_cast<std::size_t>(n) * n;
}

std::map<std::string, double> read_expectations(const std::string &directory) {
  const auto path = fs::path(directory) / "customfdm";
  std::ifstream file(path);
  if (!file) throw std::runtime_error("Could not open " + path.string());
  const std::string contents{std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
  if (file.bad()) throw std::runtime_error("Could not read " + path.string());
  return read_final_values(contents);
}

// Integrate the piecewise linear spectrum against the thermal kernel. On a
// logarithmic mesh, mass-lumped cells can be much wider than the kernel.
void thermal_convolution(NRG::Numerics::Table &density, double width, std::size_t jobs) {
  if (!(std::isfinite(width) && width > 0.0)) throw std::invalid_argument("Invalid FDM thermal width.");
  const auto &mesh = density.x;
  const auto source = density.y;
  const auto thermal_tail = [](double u) {
    const double exp_negative = std::exp(-std::abs(u));
    return (u < 0.0 ? 1.0 : exp_negative) / (1.0 + exp_negative);
  };
  constexpr double cutoff = 36.0;
  constexpr std::size_t order = 8;
  std::array<double, order> nodes{}, weights{};
  const std::unique_ptr<gsl_integration_glfixed_table, decltype(&gsl_integration_glfixed_table_free)>
      quadrature(gsl_integration_glfixed_table_alloc(order), gsl_integration_glfixed_table_free);
  if (!quadrature) throw std::bad_alloc{};
  for (std::size_t i = 0; i < order; ++i)
    gsl_integration_glfixed_point(0., 1., i, &nodes[i], &weights[i], quadrature.get());
  NRG::Tools::for_each_index(mesh.size(), jobs, [&](const std::size_t i) {
    const double x = mesh[i];
    const auto scaled = [&](double y) { return (y - x) / width; };
    // Continue the endpoint values across the two ends of the sampled window.
    double value = source.front() * thermal_tail(-scaled(mesh.front())) +
                   source.back() * thermal_tail(scaled(mesh.back()));
    const auto first = std::lower_bound(mesh.begin(), mesh.end(), x - cutoff * width) - mesh.begin();
    for (auto j = std::max<decltype(first)>(1, first); j < std::ssize(mesh) && mesh[j - 1] < x + cutoff * width; ++j) {
      const auto a = std::max(-cutoff, scaled(mesh[j - 1])), b = std::min(cutoff, scaled(mesh[j]));
      const auto pieces = std::max(1, static_cast<int>(std::ceil((b - a) / 2.0)));
      const auto step = (b - a) / pieces;
      for (int p = 0; p < pieces; ++p)
        for (std::size_t k = 0; k < order; ++k) {
          const auto u = a + (p + nodes[k]) * step;
          const auto t = std::clamp((x + width * u - mesh[j - 1]) / (mesh[j] - mesh[j - 1]), 0., 1.);
          const auto q = thermal_tail(std::abs(u));
          value += step * weights[k] * q * (1.0 - q) * std::lerp(source[j - 1], source[j], t);
        }
    }
    density.y[i] = value;
  });
}

void broaden_pair(FdmMatrix &dst, std::size_t stride, std::size_t offset, const FdmResult &result,
                  const std::vector<std::string> &directories, const NRG::Numerics::BroadeningOptions &options,
                  const std::string &kind, const std::string &first, const std::string &second) {
  const auto name = spectrum_file_name(kind, first + "-" + second, "FDM") + ".bin";
  std::vector<std::string> paths;
  for (const auto &directory : directories) paths.push_back((fs::path(directory) / name).string());
  auto primary = options;
  primary.final_derfd_gamma.reset();
  auto density = NRG::Numerics::broaden_files(paths, result.omega, primary);
  if (options.final_derfd_gamma)
    thermal_convolution(density, *options.final_derfd_gamma * options.temperature, options.jobs);
  const auto kk = NRG::Numerics::kramers_kronig(density, options.jobs);
  for (std::size_t k = 0; k < result.omega.size(); ++k)
    dst[k * stride + offset] = -std::numbers::pi * Complex{kk.y[k], density.y[k]};
}

} // namespace

FdmMatrix symmetric_self_energy(unsigned orbitals, std::span<const Complex> hartree,
                               std::span<const Complex> green, std::span<const Complex> left,
                               std::span<const Complex> right, std::span<const Complex> composite) {
  const auto size = square_size(orbitals);
  const auto map = [size, orbitals](auto values) {
    if (values.size() != size || std::any_of(values.begin(), values.end(), [](auto v) {
          return !std::isfinite(v.real()) || !std::isfinite(v.imag());
        }))
      throw std::invalid_argument("Wrong dimensions or nonfinite entries in the FDM estimator.");
    return Eigen::Map<const Matrix>(values.data(), orbitals, orbitals);
  };

  const Eigen::FullPivLU<Matrix> factor(map(green));
  if (!factor.isInvertible()) throw std::runtime_error("Singular impurity Green function in the FDM estimator.");
  const Matrix sigma = map(hartree) + map(composite) - map(left) * factor.solve(map(right));
  if (!sigma.allFinite()) throw std::runtime_error("Nonfinite self-energy in the FDM estimator.");
  return {sigma.data(), sigma.data() + size};
}

FdmResult collect_fdm(const std::vector<std::string> &directories, const FdmOperators &operators,
                      const std::vector<double> &omega, const NRG::Numerics::BroadeningOptions &broadening) {
  if (directories.empty() || omega.empty() || operators.electron.size() > std::numeric_limits<unsigned>::max())
    throw std::invalid_argument("FDM collection needs a mesh, orbital operators and twist directories.");
  FdmResult result;
  result.omega = omega;
  result.orbitals = static_cast<unsigned>(operators.electron.size());
  const auto size = square_size(result.orbitals);
  if (operators.interaction.size() != result.orbitals || operators.hartree.size() != size ||
      (!operators.spin.empty() && operators.spin.size() != result.orbitals) ||
      omega.size() > std::numeric_limits<std::size_t>::max() / size)
    throw std::invalid_argument("Incompatible FDM operator or mesh dimensions.");

  result.expectations = read_expectations(directories.front());
  for (std::size_t twist = 1; twist < directories.size(); ++twist) {
    const auto values = read_expectations(directories[twist]);
    if (values.size() != result.expectations.size()) throw std::runtime_error("Inconsistent FDM expectation columns.");
    for (const auto &[key, value] : values) {
      const auto it = result.expectations.find(key);
      if (it == result.expectations.end()) throw std::runtime_error("Unexpected FDM expectation: " + key);
      it->second += value;
    }
  }
  for (auto &[name, value] : result.expectations) value /= static_cast<double>(directories.size());
  for (const auto &name : operators.hartree) {
    const auto it = result.expectations.find(name);
    if (it == result.expectations.end()) throw std::runtime_error("Missing Hartree operator: " + name);
    result.hartree.emplace_back(it->second);
  }

  for (auto *matrix : {&result.green, &result.left, &result.right, &result.composite, &result.self_energy})
    matrix->resize(omega.size() * size);
  if (!operators.spin.empty()) result.spin.resize(omega.size() * size);

  for (unsigned i = 0; i < result.orbitals; ++i)
    for (unsigned j = 0; j < result.orbitals; ++j) {
      const auto index = static_cast<std::size_t>(i) * result.orbitals + j;
      const auto store = [&](FdmMatrix &matrix, const std::string &kind, const auto &a, const auto &b) {
        broaden_pair(matrix, size, index, result, directories, broadening, kind, a[i], b[j]);
      };
      store(result.green, "specd", operators.electron, operators.electron);
      store(result.left, "specd", operators.interaction, operators.electron);
      store(result.right, "specd", operators.electron, operators.interaction);
      store(result.composite, "specd", operators.interaction, operators.interaction);
      if (!operators.spin.empty()) store(result.spin, "spect", operators.spin, operators.spin);
    }

  const auto matrix = [size](const FdmMatrix &values, std::size_t k) {
    return std::span(values).subspan(k * size, size);
  };
  for (std::size_t k = 0; k < omega.size(); ++k) {
    try {
      const auto sigma = symmetric_self_energy(result.orbitals, result.hartree, matrix(result.green, k),
                                               matrix(result.left, k), matrix(result.right, k), matrix(result.composite, k));
      std::ranges::copy(sigma, result.self_energy.begin() + k * size);
    } catch (const std::exception &e) {
      throw std::runtime_error("FDM self-energy at mesh index " + std::to_string(k) + ": " + e.what());
    }
  }
  return result;
}

} // namespace NRG::Workflow
