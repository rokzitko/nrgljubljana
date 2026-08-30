#ifndef _tools_common_parser_hpp_
#define _tools_common_parser_hpp_

#include <cerrno>
#include <charconv>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <ios>
#include <istream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>

namespace NRG::Tools {

inline std::string trim_whitespace(const std::string &value) {
  const auto first = value.find_first_not_of(" \t\n\r\f\v");
  if (first == std::string::npos) return {};
  const auto last = value.find_last_not_of(" \t\n\r\f\v");
  return value.substr(first, last - first + 1);
}

inline auto parse_parameter_double(const std::string_view value, const std::string_view name) {
  const std::string text(value);
  char *end = nullptr;
  errno = 0;
  const auto result = std::strtod(text.c_str(), &end);
  const bool underflowed_to_zero = errno == ERANGE && result == 0.0;
  if (underflowed_to_zero || end == text.c_str() || end != text.c_str() + text.size() || !std::isfinite(result))
    throw std::invalid_argument(std::string(name) + " must be a finite representable number: " + text);
  return result;
}

inline auto parse_parameter_int(const std::string_view value, const std::string_view name) {
  int result = 0;
  const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), result);
  if (error != std::errc{} || end != value.data() + value.size())
    throw std::invalid_argument(std::string(name) + " must be an integer: " + std::string(value));
  return result;
}

inline bool find_block(std::ifstream &stream, const std::string &name) {
  const std::string target = "[" + name + "]";
  stream.clear();
  stream.seekg(0, std::ios::beg);
  while (stream) {
    std::string line;
    std::getline(stream, line);
    if (stream && target == line) break;
  }
  return bool(stream);
}

template<typename Map>
void parse_key_value_block(std::ifstream &stream, Map &params) {
  while (stream) {
    std::string line;
    std::getline(stream, line);
    if (!stream) break;
    line = trim_whitespace(line);
    if (line.empty()) continue;
    if (line[0] == '[') break;
    if (line[0] == '#') continue;
    const auto pos_eq = line.find_first_of('=');
    if (pos_eq == std::string::npos) continue;
    const auto keyword = trim_whitespace(line.substr(0, pos_eq));
    const auto value   = trim_whitespace(line.substr(pos_eq + 1));
    params.insert_or_assign(keyword, value);
  }
}

template<typename Map, typename T, typename Converter>
T get_or_default(const Map &params, const std::string &keyword, const T &def, Converter convert) {
  if (const auto it = params.find(keyword); it != params.end())
    return convert(it->second);
  return def;
}

} // namespace NRG::Tools

#endif
