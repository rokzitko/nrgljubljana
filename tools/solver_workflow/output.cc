#include <nrgljubljana/solver_output.hpp>

#include <spectrum_names.hpp>

#include <cerrno>
#include <charconv>
#include <clocale>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <string_view>
#include <system_error>
#include <utility>
#include <vector>

#ifdef __APPLE__
#include <xlocale.h>
#endif

namespace NRG::Workflow {
namespace {

constexpr std::string_view whitespace = " \t\r\n\f\v";

std::string_view trim(std::string_view s) {
  const auto first = s.find_first_not_of(whitespace);
  if (first == std::string_view::npos) return {};
  const auto last = s.find_last_not_of(whitespace);
  return s.substr(first, last - first + 1);
}

std::vector<std::string_view> fields(std::string_view s) {
  std::vector<std::string_view> result;
  std::size_t pos = 0;
  while ((pos = s.find_first_not_of(whitespace, pos)) != std::string_view::npos) {
    const auto end = s.find_first_of(whitespace, pos);
    result.push_back(s.substr(pos, end == std::string_view::npos ? std::string_view::npos : end - pos));
    pos = end;
  }
  return result;
}

// The lines of a text with their numbers
class Lines {
 public:
  explicit Lines(std::string_view text) : text_(text) {}
  bool next(std::string_view &line) {
    if (text_.empty()) return false;
    const auto eol = text_.find('\n');
    line           = text_.substr(0, eol);
    text_          = eol == std::string_view::npos ? std::string_view{} : text_.substr(eol + 1);
    ++number_;
    return true;
  }
  [[noreturn]] void fail(const std::string &what) const {
    throw std::runtime_error((number_ > 0 ? "line " + std::to_string(number_) + ": " : std::string()) + what);
  }

 private:
  std::string_view text_;
  std::size_t number_ = 0;
};

template <typename T> bool parse_number(std::string_view s, T &value) {
  const auto [end, error] = std::from_chars(s.data(), s.data() + s.size(), value);
  return error == std::errc{} && end == s.data() + s.size();
}

// Floating-point from_chars is unavailable before macOS 26 with libc++.
bool parse_number(std::string_view s, double &value) {
  if (s.empty() || s.front() == '+' || whitespace.find(s.front()) != std::string_view::npos) return false;
  const auto magnitude = s.front() == '-' ? s.substr(1) : s;
  if (magnitude.starts_with("0x") || magnitude.starts_with("0X")) return false;

  struct NumericLocale {
    locale_t handle = newlocale(LC_ALL_MASK, "C", nullptr);
    ~NumericLocale() {
      if (handle) freelocale(handle);
    }
  };
  static const NumericLocale locale;
  if (!locale.handle) throw std::runtime_error("cannot create the C locale for parsing solver output");

  const std::string text(s);
  char *end = nullptr;
  const int saved_errno = errno;
  errno = 0;
  const double parsed = strtod_l(text.c_str(), &end, locale.handle);
  const bool out_of_range = errno == ERANGE && (parsed == 0.0 || std::isinf(parsed));
  errno = saved_errno;
  // Some implementations report ERANGE even for representable subnormals.
  if (end != text.c_str() + text.size() || out_of_range) return false;
  value = parsed;
  return true;
}

} // namespace

StepTable read_step_table(const std::string &text) {
  Lines in(text);
  StepTable table;
  std::string_view raw;
  while (in.next(raw)) {
    const auto line = trim(raw);
    if (line.empty()) continue;
    if (line.front() == '#') {
      if (!table.rows.empty()) in.fail("a header line after the data");
      table.columns.clear();
      for (const auto f : fields(line.substr(1))) table.columns.emplace_back(f);
      continue;
    }
    if (table.columns.empty()) in.fail("data before a header line with the column names");
    const auto values = fields(line);
    if (values.size() != table.columns.size())
      in.fail("expected " + std::to_string(table.columns.size()) + " values, found " + std::to_string(values.size()));
    auto &row = table.rows.emplace_back();
    for (const auto v : values) {
      double x = 0;
      if (!parse_number(v, x)) in.fail("not a number: \"" + std::string(v) + "\"");
      row.push_back(x);
    }
  }
  if (table.columns.empty()) throw std::runtime_error("no header line with the column names");
  return table;
}

std::map<std::string, double> read_final_values(const std::string &text) {
  const auto table = read_step_table(text);
  if (table.rows.size() != 1) throw std::runtime_error("expected exactly one row, found " + std::to_string(table.rows.size()));
  std::map<std::string, double> result;
  for (std::size_t c = 0; c < table.columns.size(); c++) result[table.columns[c]] = table.rows[0][c];
  return result;
}

std::vector<std::vector<EnergyLevel>> read_energy_flow(const std::string &text) {
  Lines in(text);
  std::vector<std::vector<EnergyLevel>> flow;
  std::vector<EnergyLevel> step;
  std::string_view raw;
  while (in.next(raw)) {
    auto line = trim(raw);
    if (line.empty()) { // the end of a step; a step without levels is a lone empty line
      flow.push_back(std::move(step));
      step.clear();
      continue;
    }
    EnergyLevel level;
    const auto first_end = line.find_first_of(whitespace);
    if (first_end == std::string_view::npos || !parse_number(line.substr(0, first_end), level.energy))
      in.fail("expected an energy, then (quantum numbers) and [degeneracy]");
    line = trim(line.substr(first_end));
    if (line.empty() || line.front() != '(')
      in.fail("expected (quantum numbers) after the energy; annotated.dat written with dumpgroups=false is not supported");
    while (!line.empty() && line.front() == '(') {
      const auto close = line.find(')');
      if (close == std::string_view::npos) in.fail("unclosed '('");
      auto &qn = level.subspaces.emplace_back();
      for (const auto v : fields(line.substr(1, close - 1))) {
        int x = 0;
        if (!parse_number(v, x)) in.fail("not an integer quantum number: \"" + std::string(v) + "\"");
        qn.push_back(x);
      }
      line = trim(line.substr(close + 1));
    }
    if (line.size() < 3 || line.front() != '[' || line.back() != ']' || !parse_number(line.substr(1, line.size() - 2), level.degeneracy)
        || level.degeneracy < 1)
      in.fail("expected [degeneracy] at the end");
    step.push_back(std::move(level));
  }
  if (!step.empty()) flow.push_back(std::move(step));
  return flow;
}

std::string spectrum_file_name(const std::string &list, const std::string &correlator, const std::string &algorithm, const int spin) {
  static const std::map<std::string, std::string> prefixes{{"specs", "corr"},   {"spect", "spin"},  {"specot", "orbspin"},
                                                           {"specd", "spec"},   {"specq", "specq"}, {"specgt", "gt"},
                                                           {"speci1t", "i1t"},  {"speci2t", "i2t"}, {"specchit", "chit"}};
  const auto prefix = prefixes.find(list);
  if (prefix == prefixes.end()) throw std::runtime_error("unknown spectrum list: " + list);
  const bool doublet = list == "specd" || list == "specgt" || list == "speci1t" || list == "speci2t";
  if (spin != 0 && !(doublet && (spin == 1 || spin == -1)))
    throw std::runtime_error("spin " + std::to_string(spin) + " for " + list + ": only the doublet lists take +1 or -1");
  if (correlator.empty() || algorithm.empty()) throw std::runtime_error("empty correlator or algorithm");
  return NRG::spectrum_file_stem(prefix->second, algorithm, correlator + NRG::spectrum_spin_suffix(spin));
}

} // namespace NRG::Workflow
