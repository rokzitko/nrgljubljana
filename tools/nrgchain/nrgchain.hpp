#ifndef _nrgchain_nrgchain_hpp_
#define _nrgchain_nrgchain_hpp_

#include <filesystem>
#include <map>
#include <string>
#include <vector>

namespace NRG::Tools::NrgChain {

struct WilsonChannel {
  double theta{};
  std::vector<double> xi;
  std::vector<double> zeta;
};

struct WilsonData {
  std::vector<WilsonChannel> channels;
};

enum class TableMode {
  Calculate,
  SaveOnly,
  LoadAndTridiagonalize,
};

WilsonData calculate_from_params(const std::map<std::string, std::string> &param_values,
                                 TableMode mode = TableMode::Calculate,
                                 const std::filesystem::path &output_dir = ".");
WilsonData calculate_from_file(const std::string &param_filename = "param",
                               TableMode mode = TableMode::Calculate,
                               const std::filesystem::path &output_dir = ".");

} // namespace NRG::Tools::NrgChain

#endif
