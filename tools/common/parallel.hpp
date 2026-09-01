#ifndef _tools_common_parallel_hpp_
#define _tools_common_parallel_hpp_

#include <cctype>
#include <charconv>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

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

} // namespace NRG::Tools

#endif
