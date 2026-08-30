#ifndef NRG_TOOLS_COMMON_DIAGNOSTICS_HPP
#define NRG_TOOLS_COMMON_DIAGNOSTICS_HPP

#include <iomanip>
#include <limits>
#include <ostream>
#include <sstream>
#include <string_view>

namespace NRG::Tools {

class ConfigurationReport {
  std::ostringstream report;

  public:
  explicit ConfigurationReport(const std::string_view program) {
    report << std::boolalpha << std::setprecision(std::numeric_limits<double>::max_digits10);
    report << program << ": configuration\n";
  }

  template<typename T>
  void value(const std::string_view name, const T &value) {
    report << "  " << name << '=' << value << '\n';
  }

  template<typename T>
  void resolved(const std::string_view name, const T &value, const std::string_view reason) {
    report << "  " << name << "=auto -> " << value << " (" << reason << ")\n";
  }

  void write(std::ostream &out) const { out << report.str(); }
};

} // namespace NRG::Tools

#endif
