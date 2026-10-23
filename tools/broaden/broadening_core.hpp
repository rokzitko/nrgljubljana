// Numerical broadening path shared by the CLI and the native numerics API.
#ifndef NRG_TOOLS_BROADENING_CORE_HPP
#define NRG_TOOLS_BROADENING_CORE_HPP

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <limits>
#include <map>
#include <mutex>
#include <numeric>
#include <stdexcept>
#include <thread>
#include <vector>

#include "../../c++/broadening.hpp"

namespace NRG::Broaden::Core {

using std::size_t;
inline const double R_1_SQRT2PI = 1.0 / std::sqrt(2.0 * M_PI);

inline double gaussian_kernel(const double x, const double y, const double sigma) {
  const auto d = (x - y) / sigma;
  return R_1_SQRT2PI * std::exp(-d * d / 2.0) / sigma;
}

// -d/dw f_FD(x-y), with effective temperature sigma.
inline double derfd_kernel(const double x, const double y, const double sigma) {
  const auto d = (x - y) / sigma;
  return 1.0 / ((1.0 + std::cosh(d)) * 2.0 * sigma);
}

template<typename FNC>
void for_each_index(const size_t count, const size_t requested_workers, FNC &&calculate) {
  if (count == 0) return;
  if (requested_workers == 0) throw std::invalid_argument("Worker count must be positive.");
  const auto worker_count = std::min(count, requested_workers);
  if (worker_count == 1) {
    for (size_t index = 0; index < count; ++index) calculate(index);
    return;
  }

  std::atomic<size_t> next_index{0};
  std::atomic<bool> stopped{false};
  std::exception_ptr failure;
  auto failure_index = std::numeric_limits<size_t>::max();
  std::mutex failure_mutex;
  auto run_worker = [&] {
    while (!stopped.load(std::memory_order_relaxed)) {
      const auto index = next_index.fetch_add(1, std::memory_order_relaxed);
      if (index >= count) break;
      try {
        calculate(index);
      } catch (...) {
        {
          const std::lock_guard lock(failure_mutex);
          if (!failure || index < failure_index) {
            failure = std::current_exception();
            failure_index = index;
          }
        }
        stopped.store(true, std::memory_order_relaxed);
      }
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  try {
    for (size_t worker = 0; worker < worker_count; ++worker) workers.emplace_back(run_worker);
  } catch (...) {
    stopped.store(true, std::memory_order_relaxed);
    for (auto &worker : workers) worker.join();
    throw;
  }
  for (auto &worker : workers) worker.join();
  if (failure) std::rethrow_exception(failure);
}

template<typename Frequency, typename Weight>
void accumulate(std::map<double, double> &spectrum, const size_t count, Frequency frequency, Weight weight) {
  for (size_t i = 0; i < count; ++i) spectrum[frequency(i)] += weight(i);
}

inline double average(std::map<double, double> &spectrum, const size_t count,
                      std::vector<double> &frequencies, std::vector<double> &weights) {
  auto sum = 0.0;
  for (auto &entry : spectrum) {
    const auto weight = (entry.second /= static_cast<double>(count));
    frequencies.push_back(entry.first);
    weights.push_back(weight);
    sum += weight;
  }
  return sum;
}

inline void filter(const std::vector<double> &frequencies, std::vector<double> &weights,
                   const double low, const double high, const bool nonnegative, const bool negative) {
  for (size_t j = 0; j < frequencies.size(); ++j) {
    if (std::abs(frequencies[j]) < low) weights[j] = 0.0;
    if (std::abs(frequencies[j]) > high) weights[j] = 0.0;
  }
  for (size_t j = 0; j < frequencies.size(); ++j) {
    if (!nonnegative && frequencies[j] >= 0) weights[j] = 0.0;
    if (!negative && frequencies[j] < 0) weights[j] = 0.0;
  }
}

inline double kernel(const double output, const double peak, const double alpha, const double omega0,
                     const bool gaussian, const Broadening::CrossoverMode mode, const double accumulation) {
  return gaussian ? Broadening::gaussian(output, peak, alpha)
                  : Broadening::hybrid_kernel(output, peak, alpha, omega0, mode, accumulation);
}

inline std::vector<double> broaden(const std::vector<double> &frequencies, const std::vector<double> &weights,
                                   const std::vector<double> &mesh, const double alpha, const double omega0,
                                   const bool gaussian, const Broadening::CrossoverMode mode,
                                   const double accumulation, const size_t jobs) {
  std::vector<double> result(mesh.size());
  for_each_index(mesh.size(), jobs, [&](const size_t index) {
    const auto m = mesh[index];
    // Keep the CLI's reduction grouping, including zero and signed weights.
    result[index] = std::transform_reduce(weights.begin(), weights.end(), frequencies.begin(), 0.0, std::plus<>(),
                                         [=](const auto weight, const auto freq) {
                                           return weight * kernel(m, freq, alpha, omega0, gaussian, mode, accumulation);
                                         });
  });
  return result;
}

template<typename S, typename T, typename FNC>
void convolve(const std::vector<S> &mesh, std::vector<T> &a, const double sigma,
              const size_t jobs, FNC kernel, const double cutoff_ratio = 100) {
  const auto nr_mesh = mesh.size();
  if (nr_mesh == 0) throw std::runtime_error("Output mesh is empty.");
  if (sigma <= 0.0) throw std::invalid_argument("Convolution width must be greater than 0.");
  auto b(a); // source
  for_each_index(nr_mesh, jobs, [&](const size_t i) {
    if (i == 0 || i + 1 == nr_mesh) {
      a[i] = b[i];
      return;
    }
    const auto x = mesh[i];
    if (std::abs(x) < cutoff_ratio * sigma) {
      auto sum = 0.0;
      // Interior source-cell quadrature, deliberately without normalization.
      for (size_t j = 1; j < nr_mesh - 1; j++) {
        const auto y = mesh[j];
        const auto width = (mesh[j + 1] + mesh[j]) / 2.0 - (mesh[j - 1] + mesh[j]) / 2.0;
        const auto kw = kernel(x, y, sigma) * width;
        if (!std::isfinite(kw)) throw std::overflow_error("Nonfinite convolution kernel weight.");
        sum += b[j] * kw;
      }
      a[i] = sum;
    } else {
      a[i] = b[i];
    }
  });
}

} // namespace NRG::Broaden::Core

#endif
