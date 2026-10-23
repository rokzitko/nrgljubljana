#include <nrgljubljana/numerics.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <future>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <random>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <gsl/gsl_errno.h>

namespace {
using namespace NRG::Numerics;
namespace fs = std::filesystem;
constexpr auto pi = std::numbers::pi;
constexpr auto infinity = std::numeric_limits<double>::infinity();
constexpr auto nan = std::numeric_limits<double>::quiet_NaN();
const std::vector<double> mesh{-2.0, -0.7, -0.2, -0.04, 0.03, 0.1, 0.4, 1.1, 2.0};
const std::vector<SpectrumBins> spectra{
  {{0.4, -0.7}, {-0.3, 0.2}, {0.02, -0.15}, {-0.3, 0.1}, {0.0, 0.05}},
  {{-0.8, -0.4}, {0.4, 0.2}, {0.1, 0.3}}, {}};

void check(const bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

void near(const double actual, const double expected, const double atol = 2e-14, const double rtol = 2e-13) {
  if (!std::isfinite(actual) || !std::isfinite(expected) || std::abs(actual - expected) > atol + rtol * std::abs(expected)) {
    std::cerr << std::setprecision(17) << "actual=" << actual << " expected=" << expected << '\n';
    throw std::runtime_error("Broadening oracle mismatch.");
  }
}

void near(const Table &actual, const Table &expected) {
  check(actual.x == expected.x && actual.y.size() == expected.y.size(), "Output mesh or length changed.");
  for (std::size_t i = 0; i < actual.y.size(); ++i) near(actual.y[i], expected.y[i]);
}

template<typename Exception = std::exception, typename Function>
void rejects(Function function) {
  try { function(); } catch (const Exception &) { return; }
  throw std::runtime_error("Expected an exception.");
}

struct TemporaryDirectory {
  fs::path path;
  TemporaryDirectory() {
    std::random_device random;
    do { path = fs::temp_directory_path() / ("nrg-broadening-" + std::to_string(random())); }
    while (!fs::create_directory(path));
  }
  ~TemporaryDirectory() { std::error_code error; fs::remove_all(path, error); }
};

void binary(const fs::path &path, const std::vector<double> &values, const bool trailing_byte = false) {
  std::ofstream out(path, std::ios::binary);
  out.exceptions(std::ios::failbit | std::ios::badbit);
  if (!values.empty()) out.write(reinterpret_cast<const char *>(values.data()), values.size() * sizeof(double));
  if (trailing_byte) out.put('x');
  out.close();
}

// Independent single-peak oracle; no production kernel headers in this test.
double kernel(const double output, const double peak, const BroadeningOptions &o) {
  const auto gauss = [](double x, double y, double width) {
    const auto d = (x - y) / width;
    return std::exp(-d * d) / (width * std::sqrt(pi));
  };
  if (o.kernel == BroadeningKernel::gaussian) return gauss(output, peak, o.alpha);
  double x = output, y = peak;
  if (x > o.accumulation && y > o.accumulation) { x -= o.accumulation; y -= o.accumulation; }
  if (x < -o.accumulation && y < -o.accumulation) { x += o.accumulation; y += o.accumulation; }
  double logarithmic = 0.0;
  if (y != 0.0 && x * y > 0.0) {
    const auto d = std::log(x / y) / o.alpha - o.alpha / 4.0;
    logarithmic = std::exp(-d * d) / (o.alpha * std::abs(x) * std::sqrt(pi));
  }
  const auto omega0 = o.omega0_ratio * o.temperature;
  const auto frequency = std::abs(o.kernel == BroadeningKernel::hybrid_peak_frequency ? peak : output);
  const auto h = frequency > omega0 ? 1.0 : std::exp(-std::pow(std::log(frequency / omega0) / o.alpha, 2));
  return logarithmic * h + gauss(output, peak, omega0) * (1.0 - h);
}

void kernels() {
  for (const auto kind : {BroadeningKernel::hybrid_output_frequency, BroadeningKernel::hybrid_peak_frequency,
                          BroadeningKernel::gaussian}) {
    for (const auto accumulation : {0.0, 0.07}) {
      BroadeningOptions o{.alpha = 0.4, .temperature = 0.2, .omega0_ratio = 1.0, .kernel = kind, .accumulation = accumulation};
      for (const auto peak : {-0.5, -0.01, 0.0, 0.08, 0.4}) {
        const auto result = broaden({{{peak, -0.7}}}, mesh, o);
        check(result.x == mesh, "Broadening changed mesh labels.");
        for (std::size_t i = 0; i < mesh.size(); ++i) near(result.y[i], -0.7 * kernel(mesh[i], peak, o));
      }
      const auto odd = broaden({{{-0.4, -1.0}, {0.4, 1.0}}}, {-1.0, -0.1, 0.1, 1.0}, o);
      check(odd.y[0] < 0.0 && odd.y[3] > 0.0, "Signed weights were clipped.");
      near(odd.y[0], -odd.y[3]);
      near(odd.y[1], -odd.y[2]);
    }
  }
  const BroadeningOptions gaussian{.alpha = 0.3, .temperature = 0.2, .kernel = BroadeningKernel::gaussian};
  near(broaden({{{0.0, -2.0}}}, {0.0}, gaussian).y[0], -2.0 / (0.3 * std::sqrt(pi)));
  BroadeningOptions output{.alpha = 0.4, .temperature = 0.2, .omega0_ratio = 1};
  auto peak = output;
  peak.kernel = BroadeningKernel::hybrid_peak_frequency;
  check(std::abs(broaden({{{0.5, 1}}}, {0.1}, output).y[0] - broaden({{{0.5, 1}}}, {0.1}, peak).y[0]) > 0.01,
        "Hybrid crossover modes are indistinguishable.");
}

void averaging() {
  BroadeningOptions o{.alpha = 0.4, .temperature = 0.2, .kernel = BroadeningKernel::gaussian};
  const std::vector<SpectrumBins> merged{{{-0.8, -0.4 / 3}, {-0.3, (0.2 + 0.1) / 3}, {0, 0.05 / 3},
                                         {0.02, -0.15 / 3}, {0.1, 0.3 / 3}, {0.4, (-0.7 + 0.2) / 3}}};
  near(broaden(spectra, mesh, o), broaden(merged, mesh, o));
  const auto result = broaden({{{0.4, 1e16}}, {{0.4, -1e16}}, {{0.4, 1}}}, mesh, o);
  near(result, broaden({{{0.4, 1.0 / 3}}}, mesh, o));
  const auto reordered = broaden({{{0.4, 1e16}}, {{0.4, 1}}, {{0.4, -1e16}}}, mesh, o);
  check(reordered.y == std::vector<double>(mesh.size(), 0.0), "Averaging did not retain input summation order.");
  const auto adjacent = std::nextafter(1.0, 2.0);
  auto narrow = o;
  narrow.alpha = 1e-16;
  const auto exact_keys = broaden({{{1.0, 1.0}, {adjacent, -1.0}}}, {1.0, adjacent}, narrow);
  near(exact_keys.y[0], kernel(1.0, 1.0, narrow) - kernel(1.0, adjacent, narrow));
  check(exact_keys.y[0] > 1e15, "Distinct floating-point keys were coalesced.");
  for (const auto &empty : {std::vector<SpectrumBins>{{}}, std::vector<SpectrumBins>{{}, {}}})
    check(broaden(empty, mesh, o).y == std::vector<double>(mesh.size(), 0.0), "Empty contribution is not zero.");

  const std::vector<SpectrumBins> boundaries{{{-0.5, -2}, {-0.2, -1}, {0.0, -3}, {0.2, 4}, {0.5, 5}}};
  o.filter_low = 0.2;
  o.filter_high = 0.5;
  near(broaden(boundaries, mesh, o), broaden({{{-0.5, -2}, {-0.2, -1}, {0.2, 4}, {0.5, 5}}}, mesh, o));
  o.keep_negative = false;
  near(broaden(boundaries, mesh, o), broaden({{{0.2, 4}, {0.5, 5}}}, mesh, o));
  o.keep_nonnegative = false;
  check(broaden(boundaries, mesh, o).y == std::vector<double>(mesh.size(), 0.0), "Sign filters did not remove all weights.");
  o.keep_negative = true;
  o.filter_low = 0;
  near(broaden(boundaries, mesh, o), broaden({{{-0.5, -2}, {-0.2, -1}}}, mesh, o));
  o.keep_nonnegative = true;
  o.keep_negative = false;
  near(broaden(boundaries, mesh, o), broaden({{{0.0, -3}, {0.2, 4}, {0.5, 5}}}, mesh, o));
}

Table convolution_oracle(const Table &source, const double sigma, const bool derfd) {
  auto result = source;
  for (std::size_t i = 1; i + 1 < source.x.size(); ++i) {
    if (std::abs(source.x[i]) >= 100 * sigma) continue;
    result.y[i] = 0;
    for (std::size_t j = 1; j + 1 < source.x.size(); ++j) {
      const auto d = (source.x[i] - source.x[j]) / sigma;
      const auto k = derfd ? 1.0 / (2 * sigma * (1 + std::cosh(d))) : std::exp(-d * d / 2) / (sigma * std::sqrt(2 * pi));
      const auto width = (source.x[j + 1] + source.x[j]) / 2 - (source.x[j - 1] + source.x[j]) / 2;
      result.y[i] += source.y[j] * k * width;
    }
  }
  return result;
}

void convolution() {
  const std::vector<double> grid{-20, -10, -9.999, -0.8, -0.2, 0.05, 0.2, 0.9, 9.999, 10, 20};
  const std::vector<SpectrumBins> bins{{{-10, -0.5}, {-0.3, 0.2}, {0.1, -0.7}, {10, 1}, {20, 5}}};
  BroadeningOptions o{.alpha = 0.4, .temperature = 0.2, .kernel = BroadeningKernel::gaussian};
  const auto primary = broaden(bins, grid, o);
  for (const bool gaussian : {false, true}) {
    for (const bool derfd : {false, true}) {
      o.final_gaussian_gamma = gaussian ? std::optional<double>(0.5) : std::nullopt;
      o.final_derfd_gamma = derfd ? std::optional<double>(0.5) : std::nullopt;
      auto expected = primary;
      if (gaussian) expected = convolution_oracle(expected, 0.1, false);
      if (derfd) expected = convolution_oracle(expected, 0.1, true);
      const auto actual = broaden(bins, grid, o);
      near(actual, expected);
      for (const auto i : {0, 1, 9, 10}) check(actual.y[i] == primary.y[i], "Endpoint or cutoff boundary changed.");
      o.jobs = 4;
      check(broaden(bins, grid, o).y == actual.y, "Parallel convolution changed arithmetic.");
      o.jobs = 1;
    }
  }
  const auto reversed = convolution_oracle(convolution_oracle(primary, 0.1, true), 0.1, false);
  check(std::abs(broaden(bins, grid, o).y[4] - reversed.y[4]) > 1e-6, "Convolution fixture does not distinguish pass order.");
  for (const auto &small : {std::vector<double>{0.1}, std::vector<double>{-0.1, 0.1}}) {
    const auto final = broaden(bins, small, o);
    auto no_final = o;
    no_final.final_gaussian_gamma.reset();
    no_final.final_derfd_gamma.reset();
    check(final.y == broaden(bins, small, no_final).y, "One/two-point endpoints changed.");
  }
}

void parser() {
  TemporaryDirectory temp;
  const auto path = temp.path / "spectrum.bin";
  binary(path, {0.4, -2, -0.3, 1, 0.4, 0.25});
  auto bins = load_spectrum_bins(path.string());
  check(bins.size() == 3 && bins[0].frequency == 0.4 && bins[0].weight == -2 && bins[2].weight == 0.25,
        "Pair parser reordered, merged or clipped bins.");
  const BroadeningOptions o{.alpha = 0.4, .temperature = 0.2};
  const auto empty = temp.path / "empty.bin";
  binary(empty, {});
  check(load_spectrum_bins(empty.string()).empty(), "Empty file was not accepted.");
  near(broaden_files({path.string(), empty.string()}, mesh, o), broaden({bins, {}}, mesh, o));
  binary(path, {0.4, -2, 3, -0.3, 1, -4, 0.4, 0.25, -0.5});
  const auto real = load_spectrum_bins(path.string(), SpectrumLayout::complex_real);
  const auto imag = load_spectrum_bins(path.string(), SpectrumLayout::complex_imaginary);
  check(real.size() == 3 && real[0].weight == -2 && real[2].weight == 0.25, "Complex real component mismatch.");
  check(imag.size() == 3 && imag[0].weight == 3 && imag[1].weight == -4, "Complex imaginary component mismatch.");
  for (const auto layout : {SpectrumLayout::complex_real, SpectrumLayout::complex_imaginary})
    near(broaden_files({path.string(), empty.string()}, mesh, o, layout),
         broaden({layout == SpectrumLayout::complex_real ? real : imag, {}}, mesh, o));
  binary(path, {0.1, 1}, true);
  rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string()); });
  binary(path, {0.1});
  rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string()); });
  binary(path, {0.1, 1});
  rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string(), SpectrumLayout::complex_real); });
  for (const auto value : {nan, infinity, -infinity}) {
    for (std::size_t component = 0; component < 3; ++component) {
      std::vector<double> record{0.1, 1, -1};
      record[component] = value;
      binary(path, record);
      for (const auto layout : {SpectrumLayout::complex_real, SpectrumLayout::complex_imaginary})
        rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string(), layout); });
      if (component < 2) {
        record.resize(2);
        binary(path, record);
        rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string()); });
      }
    }
  }
  rejects<std::runtime_error>([&] { load_spectrum_bins((temp.path / "missing").string()); });
  rejects([&] { load_spectrum_bins(temp.path.string()); });
  rejects<std::invalid_argument>([&] { load_spectrum_bins(path.string(), static_cast<SpectrumLayout>(99)); });
}

void validation() {
  const BroadeningOptions valid{.alpha = 0.4, .temperature = 0.2};
  rejects<std::invalid_argument>([&] { broaden({}, mesh, valid); });
  rejects<std::invalid_argument>([&] { broaden_files({}, mesh, valid); });
  for (const auto &grid : {std::vector<double>{}, {0.0}, {-0.1, 0.0, 0.1}, {0.1, 0.1}, {0.2, 0.1}, {nan}, {infinity}})
    rejects<std::invalid_argument>([&] { broaden(spectra, grid, valid); });
  for (const auto value : {nan, infinity, -infinity}) {
    rejects<std::invalid_argument>([&] { broaden({{{value, 1}}}, mesh, valid); });
    rejects<std::invalid_argument>([&] { broaden({{{0.1, value}}}, mesh, valid); });
  }
  for (const auto member : {&BroadeningOptions::alpha, &BroadeningOptions::temperature, &BroadeningOptions::omega0_ratio}) {
    for (const auto value : {0.0, -1.0, nan, infinity}) {
      auto o = valid;
      o.*member = value;
      rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
    }
  }
  for (const auto member : {&BroadeningOptions::final_gaussian_gamma, &BroadeningOptions::final_derfd_gamma}) {
    for (const auto value : {0.0, -1.0, nan, infinity}) {
      auto o = valid;
      o.*member = value;
      rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
    }
    for (const auto tiny : {false, true}) {
      auto o = valid;
      o.temperature = tiny ? 1e-300 : 1e300;
      o.*member = o.temperature;
      rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
    }
  }
  for (const auto value : {1e-300, 1e300}) {
    auto o = valid;
    o.temperature = o.omega0_ratio = value;
    rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  }
  for (const auto member : {&BroadeningOptions::accumulation, &BroadeningOptions::filter_low, &BroadeningOptions::filter_high}) {
    for (const auto value : {nan, infinity}) {
      auto o = valid;
      o.*member = value;
      rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
    }
  }
  auto o = valid;
  o.filter_low = -1;
  rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  o = valid;
  o.filter_high = -1;
  rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  o.filter_low = 1;
  o.filter_high = 0.5;
  rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  o = valid;
  o.jobs = 0;
  rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  o = valid;
  o.kernel = static_cast<BroadeningKernel>(99);
  rejects<std::invalid_argument>([&] { broaden(spectra, mesh, o); });
  const auto maximum = std::numeric_limits<double>::max();
  rejects<std::overflow_error>([&] { broaden({{{0.1, maximum}, {0.1, maximum}}}, mesh, valid); });
  o = valid;
  o.kernel = BroadeningKernel::gaussian;
  o.alpha = std::numeric_limits<double>::denorm_min();
  rejects<std::overflow_error>([&] { broaden({{{0.1, 1}}}, {0.1}, o); });
  o = valid;
  o.final_derfd_gamma = 1;
  o.jobs = 3;
  // Source-cell arithmetic would overflow; workers must propagate exceptions.
  rejects<std::overflow_error>([&] { broaden({{}}, {-1, 0.1, 1, maximum / 2, maximum}, o); });
}

void gsl_handler(const char *, const char *, int, int) {}

void concurrency() {
  BroadeningOptions o{.alpha = 0.4, .temperature = 0.2, .omega0_ratio = 1,
                      .final_gaussian_gamma = 0.2, .final_derfd_gamma = 0.3};
  const auto original_cwd = fs::current_path();
  const auto expected = broaden(spectra, mesh, o);
  TemporaryDirectory temp;
  std::vector<std::string> files;
  for (std::size_t z = 0; z < spectra.size(); ++z) {
    std::vector<double> values;
    for (const auto &bin : spectra[z]) { values.push_back(bin.frequency); values.push_back(bin.weight); }
    const auto path = temp.path / (std::to_string(z) + ".bin");
    binary(path, values);
    files.push_back(path.string());
  }
  const auto previous_handler = gsl_set_error_handler(gsl_handler);
  std::ostringstream output, errors;
  auto *cout_buffer = std::cout.rdbuf(output.rdbuf());
  auto *cerr_buffer = std::cerr.rdbuf(errors.rdbuf());
  try {
    std::vector<std::future<void>> calls;
    for (std::size_t job = 1; job <= 8; ++job) {
      calls.push_back(std::async(std::launch::async, [&, job] {
        auto options = o;
        options.jobs = job;
        for (int repeat = 0; repeat < 12; ++repeat) {
          check(broaden(spectra, mesh, options).y == expected.y, "Repeated/threaded broadening changed arithmetic.");
          check(broaden_files(files, mesh, options).y == expected.y, "Concurrent file loading/broadening changed arithmetic.");
          rejects<std::invalid_argument>([&] { broaden(spectra, {0.0}, options); });
        }
      }));
    }
    for (auto &call : calls) call.get();
  } catch (...) {
    std::cout.rdbuf(cout_buffer);
    std::cerr.rdbuf(cerr_buffer);
    gsl_set_error_handler(previous_handler);
    throw;
  }
  std::cout.rdbuf(cout_buffer);
  std::cerr.rdbuf(cerr_buffer);
  const auto after_handler = gsl_set_error_handler(previous_handler);
  check(after_handler == gsl_handler, "Broadening changed the GSL error handler.");
  check(output.str().empty() && errors.str().empty(), "Native broadening logged output.");
  check(fs::current_path() == original_cwd, "Native broadening changed cwd.");
}

void cli_fixtures(const fs::path &root) {
  for (int scenario = 0; scenario < 10; ++scenario) {
    const auto directory = root / std::to_string(scenario);
    fs::create_directory(directory);
    BroadeningOptions o{.alpha = 0.4, .temperature = 0.2};
    auto layout = SpectrumLayout::real_pairs;
    std::ofstream arguments(directory / "arguments.txt");
    arguments.exceptions(std::ios::failbit | std::ios::badbit);
    if (scenario == 0) { o.final_derfd_gamma = 0.3; arguments << "-x\n0.3\n"; }
    if (scenario >= 1) o.omega0_ratio = 1.0;
    if (scenario == 1 || scenario == 7) { o.kernel = BroadeningKernel::hybrid_peak_frequency; arguments << "-n\n"; }
    if (scenario == 2 || scenario == 8) { o.kernel = BroadeningKernel::gaussian; arguments << "-g\n"; }
    if (scenario == 3) { o.accumulation = 0.07; arguments << "-a\n0.07\n"; }
    if (scenario == 4) { o.keep_negative = false; o.filter_low = 0.02; o.filter_high = 0.4; arguments << "-P\n-l\n0.02\n-h\n0.4\n"; }
    if (scenario == 5) { o.keep_nonnegative = false; arguments << "-N\n"; }
    if (scenario == 6) { layout = SpectrumLayout::complex_real; arguments << "-2\n"; }
    if (scenario == 7) { layout = SpectrumLayout::complex_imaginary; arguments << "-3\n"; }
    if (scenario >= 8) {
      o.final_gaussian_gamma = 0.2;
      o.final_derfd_gamma = 0.3;
      arguments << "-x\n0.3\n-f\n0.2\n"; // CLI option order must not change pass order.
    }
    arguments << "-L\nmesh.dat\nspec.bin\n3\n0.4\n0.2\n" << std::setprecision(17) << o.omega0_ratio << '\n';
    arguments.close();
    std::vector<std::string> files;
    for (std::size_t z = 0; z < spectra.size(); ++z) {
      const auto subdir = directory / std::to_string(z + 1);
      fs::create_directory(subdir);
      std::vector<double> values;
      for (const auto &bin : spectra[z]) {
        values.push_back(bin.frequency);
        values.push_back(bin.weight);
        if (layout != SpectrumLayout::real_pairs) values.push_back(-2 * bin.weight);
      }
      const auto file = subdir / "spec.bin";
      binary(file, values);
      files.push_back(file.string());
    }
    const auto expected = broaden_files(files, mesh, o, layout);
    for (const auto jobs : {1, 3, 20}) {
      o.jobs = jobs;
      check(broaden_files(files, mesh, o, layout).y == expected.y, "File API is not deterministic across worker counts.");
    }
    std::ofstream grid(directory / "mesh.dat"), table(directory / "native.dat");
    grid.exceptions(std::ios::failbit | std::ios::badbit);
    table.exceptions(std::ios::failbit | std::ios::badbit);
    grid << std::setprecision(18);
    table << std::setprecision(18);
    for (std::size_t i = 0; i < mesh.size(); ++i) {
      grid << mesh[i] << " 0\n";
      table << mesh[i] << ' ' << expected.y[i] << '\n';
    }
    grid.close();
    table.close();
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    if (argc < 2) throw std::invalid_argument("Expected a test group.");
    const std::string group = argv[1];
    if (group == "kernels") kernels();
    else if (group == "averaging") averaging();
    else if (group == "convolution") convolution();
    else if (group == "parser") parser();
    else if (group == "validation") validation();
    else if (group == "concurrency") concurrency();
    else if (group == "cli_fixtures" && argc == 3) cli_fixtures(argv[2]);
    else throw std::invalid_argument("Unknown broadening test group.");
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
