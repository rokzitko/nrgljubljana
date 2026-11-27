#ifndef NRG_QSTZ_TRANSFER_HPP
#define NRG_QSTZ_TRANSFER_HPP

#include <cstdint>
#include <map>
#include <stdexcept>
#include <string>
#include <string_view>
#include <fmt/format.h>

#include "operators.hpp"

namespace NRG::qstz_transfer {

template<scalar S>
Mask infer(const std::string &name, const MatrixElements<S> &seed) {
  unsigned bits = none;
  for (const auto &[key, mat] : seed) {
    // Compare coefficients directly: even a nonzero component's squared norm can underflow.
    if (!(mat.array() != S{}).any()) continue;
    const auto dtz = std::int64_t(key.first.get("TZ")) - key.second.get("TZ");
    switch (dtz) {
      case -1: bits |= minus; break;
      case  0: bits |= zero; break;
      case +1: bits |= plus; break;
      default: throw std::invalid_argument(fmt::format("QSTZ doublet {}: unexpected nonzero transfer {} (bra TZ - ket TZ)", name, dtz));
    }
  }
  return static_cast<Mask>(bits);
}

inline Mask effective_mask(const std::map<std::string, Mask> &masks, const std::string &name) {
  const auto it = masks.find(name);
  return it == masks.end() ? all : it->second;
}

template<scalar S>
auto seed_masks(const CustomOp<S> &doublets, const Params &P, const char *flag) {
  std::map<std::string, Mask> masks;
  if (!flag || std::string_view(flag) == "0") return masks;
  if (std::string_view(flag) != "1") throw std::invalid_argument("NRG_QSTZ_TRANSFER_MASK must be 0 or 1");

  // These consumers observe the block structure; even a zero block can create a GT/I1T/I2T row.
  const bool dense = P.h5ops || P.checksumrules || P.logletter('f') || P.logletter('r') || P.logletter('R') ||
                     P.specgt != "" || P.speci1t != "" || P.speci2t != "";
  for (const auto &[name, seed] : doublets) {
    const auto mask = infer(name, seed);
    masks.emplace(name, dense ? all : mask);
  }
  return masks;
}

} // namespace NRG::qstz_transfer
#endif
