// parameter_sections.hpp - Line syntax of parameter files and a strict reader of all their sections

#ifndef _parameter_sections_hpp_
#define _parameter_sections_hpp_

#include <cctype>
#include <cstddef>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>

namespace NRG {

// One line of a parameter file. The line is trimmed of surrounding whitespace first.
struct ParameterLine {
  enum class Kind {
    skip,             // blank line or whole-line # comment
    header,           // [name]; name is the text between the brackets, verbatim
    malformed_header, // starts with [ but does not end with ]
    entry,            // key=value, split at the first =; key and value trimmed
    other             // anything else (no =)
  };
  Kind kind;
  std::string name; // section name (header) or key (entry)
  std::string value;
};

inline std::string_view trim_parameter_text(std::string_view s) {
  const auto space = [](const char c) { return std::isspace(static_cast<unsigned char>(c)) != 0; };
  while (!s.empty() && space(s.front())) s.remove_prefix(1);
  while (!s.empty() && space(s.back())) s.remove_suffix(1);
  return s;
}

inline ParameterLine parameter_line(std::string_view raw) {
  const auto line = trim_parameter_text(raw);
  if (line.empty() || line.front() == '#') return {ParameterLine::Kind::skip, {}, {}};
  if (line.front() == '[') {
    if (line.size() < 2 || line.back() != ']') return {ParameterLine::Kind::malformed_header, {}, {}};
    return {ParameterLine::Kind::header, std::string(line.substr(1, line.size() - 2)), {}};
  }
  const auto pos_eq = line.find('=');
  if (pos_eq == std::string_view::npos) return {ParameterLine::Kind::other, {}, {}};
  return {ParameterLine::Kind::entry, std::string(trim_parameter_text(line.substr(0, pos_eq))),
          std::string(trim_parameter_text(line.substr(pos_eq + 1)))};
}

using ParameterSection  = std::map<std::string, std::string>;
using ParameterSections = std::map<std::string, ParameterSection>;

// All sections of parameter-file text, strictly: every line is blank, a whole-line # comment, a header line that is
// exactly [name] (non-empty name, no whitespace around or inside the brackets), or key=value inside a section. A
// section given twice, an empty key and a key given twice in a section are errors too. Text accepted here is read
// the same way by the solver (one section at a time) and by the tools. Throws std::runtime_error naming the line.
inline ParameterSections read_parameter_sections(std::string_view text) {
  ParameterSections sections;
  ParameterSection *current = nullptr;
  std::size_t number = 0;
  while (!text.empty()) {
    const auto eol = text.find('\n');
    const auto raw = text.substr(0, eol);
    text = eol == std::string_view::npos ? std::string_view{} : text.substr(eol + 1);
    ++number;
    const auto fail = [&](const std::string &what) {
      throw std::runtime_error("line " + std::to_string(number) + ": " + what + ": \"" + std::string(trim_parameter_text(raw)) + "\"");
    };
    const auto line = parameter_line(raw);
    switch (line.kind) {
      case ParameterLine::Kind::skip: break;
      case ParameterLine::Kind::malformed_header: fail("malformed section header"); break;
      case ParameterLine::Kind::header: {
        if (line.name.empty()) fail("empty section name");
        if (raw != "[" + line.name + "]" || trim_parameter_text(line.name) != line.name)
          fail("whitespace in a section header, expected [name]");
        const auto [it, inserted] = sections.try_emplace(line.name);
        if (!inserted) fail("section [" + line.name + "] given twice");
        current = &it->second;
        break;
      }
      case ParameterLine::Kind::entry:
        if (current == nullptr) fail("key before the first section");
        if (line.name.empty()) fail("empty key");
        if (!current->try_emplace(line.name, line.value).second) fail("key \"" + line.name + "\" given twice");
        break;
      case ParameterLine::Kind::other: fail("expected key=value, [section] or a # comment"); break;
    }
  }
  return sections;
}

} // namespace NRG

#endif
