#include <nrgljubljana/solver_preparation.hpp>
#include "../common/adapt_parameters.hpp"

#include <charconv>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <locale>
#include <optional>
#include <regex>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

namespace NRG::Workflow {
namespace detail {
void validate_initializer_parameters(const std::map<std::string, std::string> &parameters);
}
namespace {

std::string trim(const std::string &text) {
  const auto first = text.find_first_not_of(" \t\r\n\f\v");
  if (first == std::string::npos) return {};
  return text.substr(first, text.find_last_not_of(" \t\r\n\f\v") - first + 1);
}

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error("NRG template: " + message);
}

std::vector<std::string> fields(const std::string &text) {
  std::istringstream input(text);
  input.imbue(std::locale::classic());
  std::vector<std::string> result;
  for (std::string field; input >> field;) result.push_back(field);
  return result;
}

std::size_t count(const std::string &text, const std::string &context) {
  std::size_t value = 0;
  const auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
  require(error == std::errc{} && end == text.data() + text.size(), "invalid unsigned integer in " + context);
  return value;
}

void real(const std::string &text, const std::string &context) {
  char *end = nullptr;
  const auto value = std::strtod(text.c_str(), &end);
  require(end != text.c_str() && end == text.c_str() + text.size() && std::isfinite(value), "invalid real number in " + context);
}

// Lexical bundle lookup only; never ask the filesystem to resolve a dependency.
std::string relative_name(const std::string &name) {
  require(!name.empty() && name.front() != '/' && name.back() != '/' && name.find('\0') == std::string::npos,
          "unsafe template path: " + name);
  std::string result;
  std::size_t start = 0;
  while (start < name.size()) {
    const auto end = name.find('/', start);
    const auto part = name.substr(start, end == std::string::npos ? end : end - start);
    require(part != "..", "unsafe template path: " + name);
    if (!part.empty() && part != ".") {
      if (!result.empty()) result += '/';
      result += part;
    } else {
      require(end != std::string::npos, "template path has no filename: " + name);
    }
    if (end == std::string::npos) break;
    start = end + 1;
  }
  require(!result.empty(), "template path has no filename: " + name);
  return result;
}

using Qn = std::pair<int, int>;

std::vector<Qn> quantum_numbers(const std::string &text, std::size_t expected) {
  std::istringstream input(text);
  input.imbue(std::locale::classic());
  std::vector<Qn> result;
  for (std::size_t i = 0; i < expected; ++i) {
    int charge = 0, spin = 0;
    require(bool(input >> charge >> spin), "invalid QS quantum numbers: " + text);
    require(spin > 0, "QS spin multiplicity must be positive: " + text);
    result.emplace_back(charge, spin);
  }
  std::string extra;
  require(!(input >> extra), "extra QS quantum numbers: " + text);
  return result;
}

struct DataLines {
  std::istringstream input;
  double scale = 0;

  explicit DataLines(const std::string &text) : input(text) {}

  std::optional<std::string> next() {
    for (std::string line; std::getline(input, line);) {
      line = trim(line);
      if (line.empty()) continue;
      if (line.front() != '#') return line;
      std::istringstream comment(line);
      comment.imbue(std::locale::classic());
      std::string hash, keyword;
      if (comment >> hash >> keyword && hash == "#" && keyword == "SCALE") {
        if (!(comment >> scale)) scale = 0;
      }
    }
    return std::nullopt;
  }

  std::string take(const std::string &context) {
    const auto line = next();
    require(line.has_value(), "truncated data.in while reading " + context);
    return *line;
  }
};

// Structural counterpart of instantiate/matrix_evaluator.hpp. That evaluator
// needs model variables and a generated chain; inspection deliberately needs neither.
class MatrixShape {
  const std::string &text;
  const std::string source;
  std::size_t pos = 0;

  void check(bool condition) const { require(condition, "malformed matrix in " + source); }

  char peek() {
    while (pos < text.size()) {
      if (std::isspace(static_cast<unsigned char>(text[pos]))) { ++pos; continue; }
      if (text[pos] != '#') return text[pos];
      while (pos < text.size() && text[pos] != '\n') ++pos;
    }
    return '\0';
  }

  bool consume(char ch) {
    if (peek() != ch) return false;
    ++pos;
    return true;
  }

  void scalar() {
    struct Group {
      char close;
      unsigned arity; // Zero denotes grouping, not a function call.
      unsigned arguments = 1;
    };
    std::vector<Group> closing; // Explicit stack avoids recursion on input depth.
    bool operand = true, callable = false;
    std::string identifier;
    while (true) {
      const char ch = peek();
      if (ch == '\0' || (closing.empty() && (ch == ',' || ch == '}'))) break;
      if (std::isdigit(static_cast<unsigned char>(ch)) || ch == '.') {
        check(operand);
        char *end = nullptr;
        const auto start = text.c_str() + pos;
        const auto value = std::strtod(start, &end);
        check(end != start && std::isfinite(value));
        pos += static_cast<std::size_t>(end - start);
        operand = callable = false;
      } else if (std::isalpha(static_cast<unsigned char>(ch)) || ch == '_') {
        check(operand);
        const auto start = pos++;
        while (pos < text.size() && (std::isalnum(static_cast<unsigned char>(text[pos])) || text[pos] == '_')) ++pos;
        identifier = text.substr(start, pos - start);
        require(identifier != "parse", "unsupported parse directive in " + source);
        operand = false;
        callable = true;
      } else if (ch == '(' || ch == '[') {
        check(callable || (operand && ch == '('));
        unsigned arity = 0;
        if (callable) {
          if (identifier == "coefxi" || identifier == "coefzeta") arity = 2;
          else if (identifier == "Sqrt" || identifier == "sqrt" || identifier == "exp" || identifier == "log" ||
                   identifier == "Conjugate" || identifier == "gammaPolCh") arity = 1;
          else require(false, "unsupported function " + identifier + " in " + source);
        }
        closing.push_back({ch == '(' ? ')' : ']', arity});
        ++pos;
        operand = true;
        callable = false;
      } else if (ch == ')' || ch == ']') {
        check(!operand && !closing.empty() && closing.back().close == ch);
        check(closing.back().arity == 0 || closing.back().arguments == closing.back().arity);
        closing.pop_back();
        ++pos;
        operand = callable = false;
      } else if (ch == ',') {
        check(!operand && !closing.empty() && closing.back().arguments < closing.back().arity);
        ++closing.back().arguments;
        ++pos;
        operand = true;
        callable = false;
      } else if (ch == '+' || ch == '-' || ch == '*' || ch == '/') {
        check(!operand || ch == '+' || ch == '-');
        ++pos;
        operand = true;
        callable = false;
      } else {
        check(false); // Includes shell escapes, quoted filenames and nested lists.
      }
    }
    check(!operand && closing.empty());
  }

  std::size_t row() {
    std::size_t columns = 0;
    do { scalar(); ++columns; } while (consume(','));
    check(consume('}'));
    return columns;
  }

public:
  MatrixShape(const std::string &text_, const std::string &source_) : text(text_), source(source_) {}

  void validate(std::size_t expected_rows, std::size_t expected_columns) {
    require(text.find('\0') == std::string::npos, "NUL byte in " + source);
    std::size_t rows = 1, columns = 1;
    if (!consume('{')) scalar();
    else if (peek() != '{') columns = row();
    else {
      rows = 0;
      do {
        check(consume('{'));
        const auto width = row();
        if (rows == 0) columns = width;
        check(columns == width);
        ++rows;
      } while (consume(','));
      check(consume('}'));
    }
    check(peek() == '\0');
    require(rows == expected_rows && columns == expected_columns, "matrix dimension mismatch in " + source);
  }
};

} // namespace

unsigned resolve_chain_length(const std::map<std::string, std::string> &param_section) {
  const bool has_nmax = param_section.contains("Nmax");
  const bool has_tmin = param_section.contains("Tmin");
  if (has_nmax && has_tmin) throw std::runtime_error("Specify either Nmax or Tmin, not both!");

  std::string nmax;
  if (has_nmax) {
    nmax = trim(param_section.at("Nmax"));
    auto digits = std::string_view(nmax);
    if (digits.starts_with('-')) digits.remove_prefix(1);
    if (digits.empty() || digits.find_first_not_of("0123456789") != std::string_view::npos)
      throw std::runtime_error("Nmax must be an integer: " + nmax);
  }

  const auto numeric = [&](const std::string &name, double fallback = 0.0) {
    const auto entry = param_section.find(name);
    if (entry == param_section.end()) return fallback;
    const auto text = trim(entry->second);
    static const std::regex number(R"([+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?)");
    std::istringstream input(text);
    input.imbue(std::locale::classic());
    double value = 0;
    if (!std::regex_match(text, number) || !(input >> value) || !std::isfinite(value))
      throw std::runtime_error(name + " must be a finite number: " + text);
    if (value == 0 && text.substr(0, text.find_first_of("eE")).find_first_of("123456789") != std::string::npos)
      throw std::runtime_error("Numeric underflow in " + name);
    return value;
  };

  double tmin = 0.0;
  if (has_tmin) tmin = numeric("Tmin");
  else if (param_section.contains("T") && param_section.contains("Tmin_ratio")) {
    const double temperature = numeric("T");
    const double ratio = numeric("Tmin_ratio");
    if (temperature > 0 && ratio > 0) {
      tmin = temperature * ratio;
      if (!std::isfinite(tmin) || tmin <= 0)
        throw std::runtime_error("T*Tmin_ratio must be a finite representable number");
    }
  }

  if (tmin > 0) {
    const double lambda = numeric("Lambda", 2.0);
    if (lambda <= 1.0) throw std::runtime_error("Lambda must be greater than 1");
    const double z = numeric("z", 1.0);
    if (!(0.0 < z && z <= 1.0)) throw std::runtime_error("z must satisfy 0 < z <= 1");
    const double bandrescale = numeric("bandrescale", 1.0);
    if (bandrescale <= 0.0) throw std::runtime_error("bandrescale must be positive");
    const auto discretization = param_section.contains("discretization") ? trim(param_section.at("discretization")) : "Z";
    const char scheme = discretization.empty() ? '\0' : discretization.front();
    if (scheme != 'C' && scheme != 'Y' && scheme != 'Z')
      throw std::runtime_error("Unknown discretization: " + discretization);

    const double prefactor = scheme == 'Y' ? (1.0 + 1.0 / lambda) / 2.0 : (1.0 - 1.0 / lambda) / std::log(lambda);
    const double sc0 = prefactor * std::pow(lambda, 1.0 - z);
    unsigned resolved = 0;
    while (true) {
      // Keep the scalar operation order and inclusive threshold. Reassociation,
      // logarithmic inversion or pre-multiplying bandrescale changes boundaries.
      const double scale = sc0 * std::pow(lambda, -static_cast<double>(resolved) / 2.0) * bandrescale;
      if (!std::isfinite(scale)) throw std::runtime_error("Wilson scale is not finite");
      if (scale < tmin) break;
      if (++resolved >= 999) throw std::runtime_error("Tmin resolves to Nmax >= 999, unsupported by nrginit");
    }
    if (resolved == 0) throw std::runtime_error("Tmin is above the first Wilson scale and resolves to Nmax=0");
    return resolved;
  }

  if (!has_nmax) throw std::runtime_error("Specify Nmax, a positive Tmin, or positive T and Tmin_ratio");
  // Range-check only the effective Nmax; a positive ratio cutoff can override
  // even an arbitrarily long integer, but cannot override malformed syntax.
  const auto first = nmax.find_first_not_of('0');
  nmax.erase(0, first == std::string::npos ? nmax.size() - 1 : first);
  if (nmax.size() > 3 || nmax.front() == '-' || nmax == "0" || nmax == "999")
    throw std::runtime_error("Nmax must be an integer in the range 1..998: " + nmax);
  return static_cast<unsigned>(std::stoul(nmax));
}

void validate_preparation_parameters(const std::map<std::string, std::string> &rendered_param_section) {
  try {
    const auto scheme = rendered_param_section.contains("discretization") ? trim(rendered_param_section.at("discretization")) : "Z";
    if (scheme.empty() || std::string("CYZ").find(scheme.front()) == std::string::npos)
      throw std::invalid_argument("Unknown discretization: " + scheme);
    const auto config = NRG::Tools::resolve_adapt_parameters(rendered_param_section);
    // The combined workflow uses discretize's default integral=true, which
    // overrides f_method and uses allowed_error as CQUAD's relative tolerance.
    NRG::Tools::validate_adapt_method(config.adapt, config.hardgap, config.boundary, true);
    NRG::Tools::validate_cquad_tolerances(0.0, config.allowed_error);
    if (config.adapt)
      (void)NRG::Tools::get_or_default(rendered_param_section, "loadg", false,
                                     [](const auto &value) { return NRG::parse_bool(value); });
    if (rendered_param_section.contains("Nmax")) {
      detail::validate_initializer_parameters(rendered_param_section);
    } else {
      auto resolved = rendered_param_section;
      resolved["Nmax"] = std::to_string(resolve_chain_length(rendered_param_section));
      detail::validate_initializer_parameters(resolved);
    }
  } catch (const std::invalid_argument &error) {
    throw std::runtime_error(error.what());
  }
}

TemplateInfo inspect_template_bundle(const std::map<std::string, std::string> &files) {
  const auto model = files.find("data.in");
  require(model != files.end(), "data.in is required (paths must be relative to the template root)");
  for (const auto &[name, contents] : files) {
    require(relative_name(name) == name, "noncanonical template path: " + name);
    for (auto slash = name.find('/'); slash != std::string::npos; slash = name.find('/', slash + 1))
      require(!files.contains(name.substr(0, slash)), "template path traverses a file: " + name);
  }
  const auto &text = model->second;
  require(text.find('\0') == std::string::npos, "NUL byte in data.in");

  // Match read_data_header, not DataTemplateReader's whitespace normalization:
  // CRLF/trailing whitespace on #! or symtype can end/change native parsing.
  // Sharing this with the native reader requires extracting its heavyweight header.
  std::istringstream header(text);
  header.imbue(std::locale::classic());
  int version = -1;
  std::string symmetry;
  while (header.peek() == '#') {
    header.get();
    if (header.peek() == '!') {
      header.get();
      require(bool(header >> version), "invalid template version");
    } else {
      std::string line;
      std::getline(header, line);
      require(line.find("COMPLEX") == std::string::npos, "complex templates are unsupported");
      if (line.find("symtype", 1) != std::string::npos) {
        const auto last_space = line.find_last_of(" \t");
        if (last_space != std::string::npos && last_space + 1 < line.size()) symmetry = line.substr(last_space + 1);
      }
    }
    if (header.peek() == '\n') header.get();
  }
  std::size_t channels = 0, sites = 0, subspaces = 0;
  require(version == 9 && symmetry == "QS" && bool(header >> channels >> sites >> subspaces) && channels == 1 && subspaces > 0,
          "expected a real version-9 one-channel QS template header");
  require(subspaces <= std::numeric_limits<unsigned>::max(), "subspace count exceeds TemplateInfo range");

  DataLines data(text);
  const auto dimensions = fields(data.take("template dimensions"));
  require(dimensions.size() == 3, "invalid template dimensions");
  require(count(dimensions[0], "channels") == channels && count(dimensions[1], "chain sites") == sites &&
          count(dimensions[2], "subspaces") == subspaces, "inconsistent template header");

  const auto matrix_file = [&](const std::string &name, std::size_t rows, std::size_t columns) {
    const auto path = relative_name(name);
    const auto file = files.find(path);
    require(file != files.end(), "matrix reference is not in the template bundle: " + name);
    MatrixShape(file->second, path).validate(rows, columns);
  };

  std::map<Qn, std::size_t> sizes;
  for (std::size_t i = 0; i < subspaces; ++i) {
    const auto qn = quantum_numbers(data.take("seed quantum numbers"), 1).front();
    const auto size = count(data.take("seed dimension"), "seed dimension");
    require(size > 0, "zero seed dimension");
    const auto diag = fields(data.take("DIAG reference"));
    require(diag.size() == 2 && diag[0] == "DIAG", "expected DIAG and exactly one matrix filename");
    if (i == 0)
      require(std::isfinite(data.scale) && data.scale > 0 && std::isfinite(1.0 / data.scale), "missing or invalid # SCALE before first DIAG");
    require(sizes.emplace(qn, size).second, "duplicate seed subspace");
    matrix_file(diag[1], size, size);
  }

  bool wilson = false, energy = false;
  while (const auto line = data.next()) {
    const auto type = line->front();
    // instantiate discards/regenerates this suffix, including any later blocks.
    if (type == 'z' || type == 'Z' || type == 'T') break;
    require(wilson || type == 'f', "missing leading Wilson-chain operator block");
    if (*line == "e") {
      data.take("ground-energy placeholder");
      energy = true;
      continue;
    }
    if (type == 'f') {
      require(!wilson, "more than one Wilson-chain operator block");
      const auto indices = fields(line->substr(1));
      require(indices.size() == 2, "invalid Wilson-chain operator header");
      // The native Opch reader consumes these unsigned labels without using them.
      count(indices[0], "Wilson channel");
      count(indices[1], "Wilson operator");
      wilson = true;
    } else {
      // instantiate recognizes 'v', but the subsequent native InputData does not.
      require(std::string_view("spdtogq").find(type) != std::string_view::npos, "unsupported data.in block: " + *line);
    }
    const auto elements = count(data.take("matrix-element count"), "matrix-element count");
    for (std::size_t i = 0; i < elements; ++i) {
      const auto qn = quantum_numbers(data.take("operator quantum numbers"), 2);
      require(sizes.contains(qn[0]) && sizes.contains(qn[1]), "operator references unknown seed subspace");
      const auto rows = sizes.at(qn[0]), columns = sizes.at(qn[1]);
      const auto first = data.take("operator matrix");
      if (first.front() == 'o') matrix_file(first, rows, columns);
      else {
        for (std::size_t row = 0; row < rows; ++row) {
          const auto values = fields(row == 0 ? first : data.take("inline operator matrix"));
          require(values.size() == columns, "inline operator matrix dimension mismatch");
          for (const auto &value : values) real(value, "inline operator matrix");
        }
      }
    }
  }
  require(wilson && energy, "missing Wilson-chain operator or e block in data.in");
  return {symmetry, static_cast<unsigned>(channels), static_cast<unsigned>(subspaces)};
}

} // namespace NRG::Workflow
