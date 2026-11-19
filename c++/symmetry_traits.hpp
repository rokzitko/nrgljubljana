#ifndef _symmetry_traits_hpp_
#define _symmetry_traits_hpp_

// Properties of the symmetry types that are known from the name alone, without constructing a Symmetry object.

#include <algorithm>
#include <array>
#include <cstddef>
#include <string_view>

namespace NRG {

struct SymmetryTraits {
  std::string_view name;
  bool field; // Symmetry::isfield(): doublet spectra are computed per spin, named with -u/-d
};

// Every symmetry type of mk_sym.hpp. test/unit/mk_sym checks this list against the Symmetry classes.
inline constexpr std::array<SymmetryTraits, 31> symmetry_traits{{
  {"QS", false},      {"QSZ", true},      {"ISO", false},     {"ISO2", false},    {"ISOSZ", true},
  {"SPSU2", false},   {"SPU1", true},     {"DBLQSZ", true},   {"DBLSU2", false},  {"DBLISOSZ", true},
  {"ISOLR", false},   {"ISO2LR", false},  {"ISOSZLR", true},  {"NONE", false},    {"P", false},
  {"PP", false},      {"SL", false},      {"SL3", false},     {"SPSU2LR", false}, {"SPSU2T", false},
  {"SPU1LR", true},   {"SU2", false},     {"QSLR", false},    {"QST", false},     {"QSTZ", false},
  {"QSZTZ", false},   {"QSZLR", true},    {"QJ", false},      {"U1", false},      {"SPSU2C3", false},
  {"QSC3", false},
}};

// The traits of a symmetry type, or nullptr if the name is unknown
inline const SymmetryTraits *find_symmetry_traits(const std::string_view name) {
  const auto it = std::find_if(symmetry_traits.begin(), symmetry_traits.end(), [&](const auto &t) { return t.name == name; });
  return it == symmetry_traits.end() ? nullptr : &*it;
}

// Number of f-matrices (Wilson-chain operator blocks) per channel in data
inline std::size_t perchannel(const std::string_view symtype) {
  if (symtype == "U1" || symtype == "SU2" || symtype == "DBLSU2") return 2; // spin-up and spin-down operators
  if (symtype == "NONE" || symtype == "P" || symtype == "PP") return 4;     // CR/AN and spin UP/DO
  return 1;
}

} // namespace NRG

#endif
