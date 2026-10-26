#ifndef _tools_common_parallel_hpp_
#define _tools_common_parallel_hpp_

#include <algorithm>
#include <atomic>
#include <cctype>
#include <charconv>
#include <cstddef>
#include <exception>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <thread>
#include <vector>

namespace NRG::Tools {

inline auto parse_worker_count(const std::string_view value, const std::string_view source) {
  auto trim = [](std::string_view token) {
    while (!token.empty() && std::isspace(static_cast<unsigned char>(token.front()))) token.remove_prefix(1);
    while (!token.empty() && std::isspace(static_cast<unsigned char>(token.back()))) token.remove_suffix(1);
    return token;
  };
  auto first_value = std::size_t{0};
  auto first_token = true;
  auto offset = std::size_t{0};
  while (offset <= value.size()) {
    const auto comma = value.find(',', offset);
    const auto token = trim(value.substr(offset, comma == std::string_view::npos ? comma : comma - offset));
    std::size_t parsed = 0;
    const auto [end, error] = std::from_chars(token.data(), token.data() + token.size(), parsed);
    if (token.empty() || error != std::errc{} || end != token.data() + token.size() || parsed == 0)
      throw std::invalid_argument(std::string(source) + " must be a comma-separated list of positive integers: "
                                  + std::string(value));
    if (first_token) {
      first_value = parsed;
      first_token = false;
    }
    if (comma == std::string_view::npos) break;
    offset = comma + 1;
  }
  return first_value;
}

template<typename FNC>
void for_each_index(const std::size_t count, const std::size_t requested_workers, FNC &&calculate) {
  if (count == 0) return;
  if (requested_workers == 0) throw std::invalid_argument("Worker count must be positive.");
  const auto worker_count = std::min(count, requested_workers);
  if (worker_count == 1) {
    for (std::size_t index = 0; index < count; ++index) calculate(index);
    return;
  }

  std::atomic<std::size_t> next_index{0};
  std::atomic<bool> stopped{false};
  std::exception_ptr failure;
  auto failure_index = std::numeric_limits<std::size_t>::max();
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
    for (std::size_t worker = 0; worker < worker_count; ++worker) workers.emplace_back(run_worker);
  } catch (...) {
    stopped.store(true, std::memory_order_relaxed);
    for (auto &worker : workers) worker.join();
    throw;
  }
  for (auto &worker : workers) worker.join();
  if (failure) std::rethrow_exception(failure);
}

} // namespace NRG::Tools

#endif
