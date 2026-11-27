#ifndef NRG_QSTZ_TRANSFER_HPP
#define NRG_QSTZ_TRANSFER_HPP

#include <cstdint>
#include <stdexcept>
#include <string>
#include <fmt/format.h>

#include "doublet_transfer.hpp"

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

using doublet_transfer::effective_mask;

template<scalar S>
auto seed_masks(const CustomOp<S> &doublets, const Params &P, const char *flag) {
  return doublet_transfer::seed_masks(doublets, P, flag, "NRG_QSTZ_TRANSFER_MASK", infer<S>);
}

} // namespace NRG::qstz_transfer
#endif
