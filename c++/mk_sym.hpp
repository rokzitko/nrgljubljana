#ifndef _mk_sym_hpp_
#define _mk_sym_hpp_

#include <complex>
#include <cstddef>
#include <iostream>
#include <memory>
#include <ostream>
#include <string>
#include <stdexcept>
#include <type_traits>

#include <fmt/format.h>

#include "params.hpp"
#include "symmetry.hpp"
#include "outfield.hpp"

#include "symmetry_traits.hpp"

namespace NRG {

// Factory functions, defined in the sources generated from sym-instance.cc.in. Only the types of the build sets that
// are compiled in are defined, and complex-only types only for complex S.
#define NRG_SYMMETRY(name, ...) template <scalar S> std::unique_ptr<Symmetry<S>> mk_##name(const Params &P);
#include "symmetries.def"
#undef NRG_SYMMETRY

template<scalar S>
std::unique_ptr<Symmetry<S>> get(const std::string &sym_string, const Params &P)
{
#define NRG_SYMMETRY(name, impl, set, field, perchannel, spins, chmin, chmax, complex_only, ...)          \
  if constexpr (symmetry_set_##set && (!(complex_only) || std::is_same_v<S, std::complex<double>>)) \
    if (sym_string == #name) return mk_##name<S>(P);
#include "symmetries.def"
#undef NRG_SYMMETRY
  const auto *traits = find_symmetry_traits(sym_string);
  if (!traits) throw std::runtime_error("Unknown symmetry " + sym_string);
  if (!traits->compiled) throw std::runtime_error("Symmetry type " + sym_string + " is not included in this build (SYM_MORE, SYM_ALL)");
  throw std::runtime_error("Symmetry type " + sym_string + " requires complex data (options=COMPLEX)");
}

// Reject the combinations of a symmetry type with the number of channels and with the chain options that it does not
// support, before load() fills In and QN (which have P.combs+1 entries).
inline void check_symmetry_support(const Params &P, const std::string &sym_string, const size_t channels) {
  const auto *traits = find_symmetry_traits(sym_string);
  if (!traits) return; // get() throws
  if (P.substeps && !traits->substeps)
    throw std::invalid_argument("substeps=true is not supported for symmetry type " + sym_string);
  if (!P.substeps && !traits->supports_channels(channels))
    throw std::invalid_argument(fmt::format("Symmetry type {} supports {} channels, not {}.", sym_string,
                                            traits->chmax == 0 ? fmt::format("{} or more", traits->chmin)
                                            : traits->chmin == traits->chmax ? fmt::format("{}", traits->chmin)
                                                                             : fmt::format("{} to {}", traits->chmin, traits->chmax),
                                            channels));
  if (P.polarized && !traits->polarized)
    throw std::invalid_argument("polarized=true is not supported for symmetry type " + sym_string);
}

// Called immediately after parsing the information about the number of channels from the data file. This ensures
// that Invar can be parsed correctly.
template <scalar S>
std::shared_ptr<Symmetry<S>> set_symmetry(Params &P, const std::string &sym_string, const size_t channels) {
  Invar::invdim = 0;
  Invar::qntype.clear();
  Invar::names.clear();
  P.symtype = sym_string;
  P.set_channels_and_combs(channels);
  std::cout << "\nSYMMETRY TYPE: " << sym_string << std::endl;
  auto Sym = get<S>(sym_string, P);
  check_symmetry_support(P, sym_string, channels);
  Sym->load();
  Sym->erase_first();
  return Sym;
}

} // namespace

#endif
