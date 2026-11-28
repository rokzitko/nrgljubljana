#ifndef _symmetry_traits_hpp_
#define _symmetry_traits_hpp_

// Properties of the symmetry types that are known from the name alone, without constructing a Symmetry object. The
// table is generated from symmetries.def; see there for the meaning of the fields.

#include <algorithm>
#include <array>
#include <cstddef>
#include <string_view>

namespace NRG {

// Build sets of symmetry types (c++/CMakeLists.txt)
inline constexpr bool symmetry_set_BASIC = true;
#ifdef NRG_SYM_MORE
inline constexpr bool symmetry_set_MORE = true;
#else
inline constexpr bool symmetry_set_MORE = false;
#endif
#ifdef NRG_SYM_ALL
inline constexpr bool symmetry_set_ALL = true;
#else
inline constexpr bool symmetry_set_ALL = false;
#endif

struct SymmetryTraits {
  std::string_view name;
  bool compiled;   // included in this build
  bool field;      // Symmetry::isfield(): doublet spectra are computed per spin, named with -u/-d
  std::size_t perchannel; // number of f-matrices (Wilson-chain operator blocks) per channel in data
  std::size_t spins;      // number of spin states per orbital
  std::size_t chmin, chmax; // supported numbers of channels without substeps; chmax 0: any
  bool complex_only;
  bool substeps;
  bool polarized;
  bool pairing;
  [[nodiscard]] constexpr bool supports_channels(const std::size_t channels) const {
    return channels >= chmin && (chmax == 0 || channels <= chmax);
  }
};

inline constexpr std::array symmetry_traits{
#define NRG_SYMMETRY(name, impl, set, field, perchannel, spins, chmin, chmax, complex_only, substeps, polarized, pairing) \
  SymmetryTraits{#name, symmetry_set_##set, field, perchannel, spins, chmin, chmax, complex_only, substeps, polarized, pairing},
#include "symmetries.def"
#undef NRG_SYMMETRY
};

// The traits of a symmetry type, or nullptr if the name is unknown
inline const SymmetryTraits *find_symmetry_traits(const std::string_view name) {
  const auto it = std::find_if(symmetry_traits.begin(), symmetry_traits.end(), [&](const auto &t) { return t.name == name; });
  return it == symmetry_traits.end() ? nullptr : &*it;
}

// Number of f-matrices (Wilson-chain operator blocks) per channel in data (1 for an unknown name)
inline std::size_t perchannel(const std::string_view symtype) {
  const auto *t = find_symmetry_traits(symtype);
  return t ? t->perchannel : 1;
}

} // namespace NRG

#endif
