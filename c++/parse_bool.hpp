#ifndef _parse_bool_hpp_
#define _parse_bool_hpp_

#include <cctype>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>

namespace NRG {

inline bool parse_bool(std::string_view value) {
  const auto original = value;
  while (!value.empty() && std::isspace(static_cast<unsigned char>(value.front()))) value.remove_prefix(1);
  while (!value.empty() && std::isspace(static_cast<unsigned char>(value.back()))) value.remove_suffix(1);

  const auto equals_ignore_case = [value](const std::string_view expected) {
    if (value.size() != expected.size()) return false;
    for (std::size_t i = 0; i < value.size(); ++i) {
      if (std::tolower(static_cast<unsigned char>(value[i])) !=
          std::tolower(static_cast<unsigned char>(expected[i])))
        return false;
    }
    return true;
  };

  if (value == "1" || equals_ignore_case("true") || equals_ignore_case("yes")) return true;
  if (value == "0" || equals_ignore_case("false") || equals_ignore_case("no")) return false;
  throw std::runtime_error("Invalid boolean value [" + std::string(original) +
                           "]; expected yes/no, true/false, or 1/0.");
}

} // namespace NRG

#endif
