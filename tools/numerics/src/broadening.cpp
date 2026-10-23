#include <nrgljubljana/numerics.hpp>

#include <cmath>
#include <cstdint>
#include <fstream>
#include <map>
#include <stdexcept>

#include "../../broaden/broadening_core.hpp"

namespace NRG::Numerics {
namespace {

void positive_width(const double value) {
  if (!std::isfinite(value) || value <= 0.0)
    throw std::invalid_argument("Broadening widths and temperature must be finite and positive.");
}

void finite_weights(const std::vector<double> &weights) {
  for (const auto weight : weights)
    if (!std::isfinite(weight)) throw std::overflow_error("Nonfinite broadened spectrum or averaged weight.");
}

} // namespace

SpectrumBins load_spectrum_bins(const std::string &filename, const SpectrumLayout layout) {
  std::size_t columns, component;
  switch (layout) {
    case SpectrumLayout::real_pairs: columns = 2; component = 1; break;
    case SpectrumLayout::complex_real: columns = 3; component = 1; break;
    case SpectrumLayout::complex_imaginary: columns = 3; component = 2; break;
    default: throw std::invalid_argument("Unknown spectrum layout.");
  }
  std::ifstream input(filename, std::ios::binary | std::ios::ate);
  if (!input) throw std::runtime_error("Cannot open binary spectrum: " + filename);
  const auto end = input.tellg();
  if (end < std::streampos(0)) throw std::runtime_error("Cannot determine binary spectrum size: " + filename);
  const auto bytes = static_cast<std::uintmax_t>(static_cast<std::streamoff>(end));
  if (bytes % (columns * sizeof(double)) != 0)
    throw std::invalid_argument("Binary spectrum contains an incomplete record: " + filename);
  std::vector<double> buffer;
  if (bytes > static_cast<std::uintmax_t>(std::numeric_limits<std::streamsize>::max())
      || bytes / sizeof(double) > buffer.max_size())
    throw std::length_error("Binary spectrum is too large: " + filename);
  buffer.resize(static_cast<std::size_t>(bytes / sizeof(double)));
  input.seekg(0, std::ios::beg);
  if (!input) throw std::runtime_error("Cannot seek binary spectrum: " + filename);
  if (bytes != 0) input.read(reinterpret_cast<char *>(buffer.data()), static_cast<std::streamsize>(bytes));
  if (!input) throw std::runtime_error("Cannot read complete binary spectrum: " + filename);
  input.close();
  if (!input) throw std::runtime_error("Cannot close binary spectrum: " + filename);
  for (const auto value : buffer)
    if (!std::isfinite(value)) throw std::invalid_argument("Binary spectrum components must be finite: " + filename);
  SpectrumBins result;
  result.reserve(buffer.size() / columns);
  for (std::size_t i = 0; i < buffer.size(); i += columns) result.push_back({buffer[i], buffer[i + component]});
  return result;
}

Table broaden(const std::vector<SpectrumBins> &spectra, const std::vector<double> &mesh, const BroadeningOptions &options) {
  if (spectra.empty()) throw std::invalid_argument("At least one spectrum contribution is required.");
  positive_width(options.alpha);
  positive_width(options.temperature);
  positive_width(options.omega0_ratio);
  const auto omega0 = options.omega0_ratio * options.temperature;
  positive_width(omega0);
  for (const auto gamma : {options.final_gaussian_gamma, options.final_derfd_gamma}) {
    if (gamma) {
      positive_width(*gamma);
      positive_width(*gamma * options.temperature);
    }
  }
  if (!std::isfinite(options.accumulation)) throw std::invalid_argument("Accumulation must be finite.");
  if (!std::isfinite(options.filter_low) || !std::isfinite(options.filter_high)
      || options.filter_low < 0.0 || options.filter_high < options.filter_low)
    throw std::invalid_argument("Filters must be finite with 0 <= low <= high.");
  if (options.jobs == 0) throw std::invalid_argument("Worker count must be positive.");
  auto mode = Broadening::CrossoverMode::output_frequency;
  bool gaussian = false;
  switch (options.kernel) {
    case BroadeningKernel::hybrid_output_frequency: break;
    case BroadeningKernel::hybrid_peak_frequency: mode = Broadening::CrossoverMode::peak_frequency; break;
    case BroadeningKernel::gaussian: gaussian = true; break;
    default: throw std::invalid_argument("Unknown broadening kernel.");
  }
  if (mesh.empty()) throw std::invalid_argument("Broadening mesh must not be empty.");
  for (std::size_t i = 0; i < mesh.size(); ++i) {
    if (!std::isfinite(mesh[i]) || (i != 0 && !(mesh[i - 1] < mesh[i])))
      throw std::invalid_argument("Broadening mesh must be finite and strictly increasing.");
    if (!gaussian && mesh[i] == 0.0)
      throw std::invalid_argument("Hybrid broadening requires nonzero destination frequencies.");
  }

  std::map<double, double> spectrum;
  for (const auto &bins : spectra) {
    for (const auto &bin : bins)
      if (!std::isfinite(bin.frequency) || !std::isfinite(bin.weight))
        throw std::invalid_argument("Spectrum frequencies and weights must be finite.");
    Broaden::Core::accumulate(spectrum, bins.size(),
                             [&](const std::size_t i) { return bins[i].frequency; },
                             [&](const std::size_t i) { return bins[i].weight; });
  }
  std::vector<double> frequencies, weights;
  Broaden::Core::average(spectrum, spectra.size(), frequencies, weights);
  finite_weights(weights);
  Broaden::Core::filter(frequencies, weights, options.filter_low, options.filter_high,
                       options.keep_nonnegative, options.keep_negative);
  Table result{mesh, Broaden::Core::broaden(frequencies, weights, mesh, options.alpha, omega0, gaussian, mode,
                                          options.accumulation, options.jobs)};
  finite_weights(result.y);
  if (options.final_gaussian_gamma) {
    Broaden::Core::convolve(mesh, result.y, *options.final_gaussian_gamma * options.temperature,
                           options.jobs, Broaden::Core::gaussian_kernel);
    finite_weights(result.y);
  }
  if (options.final_derfd_gamma) {
    Broaden::Core::convolve(mesh, result.y, *options.final_derfd_gamma * options.temperature,
                           options.jobs, Broaden::Core::derfd_kernel);
    finite_weights(result.y);
  }
  return result;
}

Table broaden_files(const std::vector<std::string> &filenames, const std::vector<double> &mesh,
                    const BroadeningOptions &options, const SpectrumLayout layout) {
  if (filenames.empty()) throw std::invalid_argument("At least one spectrum file is required.");
  std::vector<SpectrumBins> spectra;
  spectra.reserve(filenames.size());
  for (const auto &filename : filenames) spectra.push_back(load_spectrum_bins(filename, layout));
  return broaden(spectra, mesh, options);
}

} // namespace NRG::Numerics
