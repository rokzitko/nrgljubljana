#ifndef NRG_DOUBLET_TRANSFER_HPP
#define NRG_DOUBLET_TRANSFER_HPP

// Optional restriction of the recalculation of named doublet operators to the transfers of quantum numbers that are
// present in their seeds (parameter doublet_transfer_mask). The symmetry type infers a mask from the seed
// (Symmetry::doublet_transfer_mask()) and recalculates only the transfers in the mask
// (Symmetry::recalc_doublet_masked()). The masks of QSTZ and QSZTZ are defined below.

#include <cstdint>
#include <functional>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <fmt/format.h>

#include "operators.hpp"
#include "params.hpp"

namespace NRG {

namespace qstz_transfer {
// Support of bra TZ - ket TZ for a QSTZ doublet.
enum Mask : unsigned { none = 0, minus = 1, zero = 2, plus = 4, all = minus | zero | plus };

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
} // namespace qstz_transfer

namespace qsztz_transfer {
enum Mask : unsigned { none = 0, all = 0x3f }; // six joint (dSZ, dTZ) pairs

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
    const auto delta = [&bra = key.first, &ket = key.second](const char *qn) { return std::int64_t(bra.get(qn)) - ket.get(qn); };
    const auto dq = delta("Q"), dsz = delta("SZ"), dtz = delta("TZ");
    const auto pair = bit(dsz, dtz);
    if (dq != 1 || pair == none)
      throw std::invalid_argument(fmt::format("QSZTZ doublet {}: unexpected nonzero transfer (dQ={}, dSZ={}, dTZ={}) (bra - ket)", name, dq, dsz, dtz));
    bits |= pair;
  }
  return static_cast<Mask>(bits);
}
} // namespace qsztz_transfer

namespace doublet_transfer {

// The mask of a named doublet; operators without an entry are recalculated for all transfers.
using Masks = std::map<std::string, unsigned>;

// Masks of the doublets, from their seeds. 'infer' returns the mask of a seed, or nullopt if the symmetry type has
// none. Every seed is validated, but no masks are kept when doublet_transfer_mask=false or when the run needs all
// blocks.
template<scalar S>
Masks seed_masks(const CustomOp<S> &doublets, const Params &P,
                 const std::function<std::optional<unsigned>(const std::string &, const MatrixElements<S> &)> &infer) {
  Masks masks;
  if (!P.doublet_transfer_mask) return masks;
  // These consumers observe the block structure; even a zero block can create a GT/I1T/I2T row.
  const bool dense = P.h5ops || P.checksumrules || P.logletter('f') || P.logletter('r') || P.logletter('R') ||
                     P.specgt != "" || P.speci1t != "" || P.speci2t != "";
  for (const auto &[name, seed] : doublets)
    if (const auto mask = infer(name, seed); mask && !dense) masks.emplace(name, *mask);
  return masks;
}

} // namespace doublet_transfer
} // namespace NRG
#endif
