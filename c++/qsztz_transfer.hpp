#ifndef NRG_QSZTZ_TRANSFER_HPP
#define NRG_QSZTZ_TRANSFER_HPP

#include <cstdint>
#include <stdexcept>
#include <string>
#include <fmt/format.h>

#include "doublet_transfer.hpp"

namespace NRG::qsztz_transfer {

// The low and high triples are independent dTZ supports for dSZ = -1 and +1.
constexpr Mask bit(std::int64_t dsz, std::int64_t dtz) {
  if ((dsz != -1 && dsz != +1) || dtz < -1 || dtz > +1) return none;
  return static_cast<Mask>(1u << ((dsz == -1 ? 0 : 3) + dtz + 1));
}

template<scalar S>
Mask infer(const std::string &name, const MatrixElements<S> &seed) {
  unsigned bits = none;
  for (const auto &[key, mat] : seed) {
    if (!(mat.array() != S{}).any()) continue; // exact, including subnormal real or imaginary components
    const auto &[bra, ket] = key;
    const auto delta = [&](const char *qn) { return std::int64_t(bra.get(qn)) - ket.get(qn); };
    const auto dq = delta("Q"), dsz = delta("SZ"), dtz = delta("TZ");
    const auto pair = bit(dsz, dtz);
    if (dq != 1 || pair == none)
      throw std::invalid_argument(fmt::format("QSZTZ doublet {}: unexpected nonzero transfer (dQ={}, dSZ={}, dTZ={}) (bra - ket)", name, dq, dsz, dtz));
    bits |= pair;
  }
  return static_cast<Mask>(bits);
}

using doublet_transfer::effective_mask;

template<scalar S>
auto seed_masks(const CustomOp<S> &doublets, const Params &P, const char *flag) {
  return doublet_transfer::seed_masks(doublets, P, flag, "NRG_QSZTZ_TRANSFER_MASK", infer<S>);
}

} // namespace NRG::qsztz_transfer
#endif
