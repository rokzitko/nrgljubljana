#ifndef _tools_common_parser_hpp_
#define _tools_common_parser_hpp_

#include <fstream>
#include <string>

namespace NRG::Tools {

inline std::string trim_whitespace(const std::string &value) {
  const auto first = value.find_first_not_of(" \t\n\r\f\v");
  if (first == std::string::npos) return {};
  const auto last = value.find_last_not_of(" \t\n\r\f\v");
  return value.substr(first, last - first + 1);
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
