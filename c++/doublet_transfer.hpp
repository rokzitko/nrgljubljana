#ifndef NRG_DOUBLET_TRANSFER_HPP
#define NRG_DOUBLET_TRANSFER_HPP

#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <fmt/format.h>

#include "operators.hpp"

namespace NRG::doublet_transfer {

template<typename Mask>
Mask effective_mask(const std::map<std::string, Mask> &masks, const std::string &name) {
  const auto it = masks.find(name);
  return it == masks.end() ? Mask::all : it->second;
}

template<scalar S, typename Mask>
auto seed_masks(const CustomOp<S> &doublets, const Params &P, const char *flag, std::string_view flag_name,
                Mask (*infer)(const std::string &, const MatrixElements<S> &)) {
  std::map<std::string, Mask> masks;
  if (!flag || std::string_view(flag) == "0") return masks;
  if (std::string_view(flag) != "1") throw std::invalid_argument(fmt::format("{} must be 0 or 1", flag_name));

  // These consumers observe the block structure; even a zero block can create a GT/I1T/I2T row.
  const bool dense = P.h5ops || P.checksumrules || P.logletter('f') || P.logletter('r') || P.logletter('R') ||
                     P.specgt != "" || P.speci1t != "" || P.speci2t != "";
  for (const auto &[name, seed] : doublets) {
    const auto mask = infer(name, seed); // validate every seed, even for a dense run
    masks.emplace(name, dense ? Mask::all : mask);
  }
  return masks;
}

} // namespace NRG::doublet_transfer
#endif
