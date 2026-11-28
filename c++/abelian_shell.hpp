// abelian_shell.hpp - Coefficient tables of abelian symmetry types, generated at run time for any number of channels.
//
// In an NRG step, the states of the added Wilson-chain site (one site per channel) are combined with the invariant
// subspaces of the previous step. For the abelian symmetry types NONE, P, U1, QSZ, and SPU1 every combination is a
// single Fock state of the added site, so all coefficients are matrix elements of fermion operators between Fock
// states. They used to be generated with Mathematica (c++/symmetry/coefnew, kept as test reference) separately for
// each channel number; this file computes the same tables for any number of channels.
//
// Conventions (identical to those of the Mathematica tables, see test/unit/symmetry/abelian_shell.cpp):
// - Modes: spin-1/2 fermions of channels 0..K-1; mode index m = 2*ch + spin with spin up (0) before spin down (1).
// - A Fock state |s> is the product of the creation operators of the occupied modes in increasing mode order acting
//   on the vacuum. Bit m of a State is the occupation of mode m.
// - A combination |i, r> of the added-site state |i> with a state |r> of the previous sites places the operators of the
//   added site to the left: |i, r> = A_i^dag |r>. Consequently an operator O acting on the previous sites has
//   <i, r|O|j, r'> = delta_ij (-1)^(n_i |O|) <r|O|r'>, where n_i is the number of electrons in |i> and |O| the
//   fermion parity of O, and the hopping terms have the signs given in hop_cr() and hop_an().
// - Combination order: the 4^K states are enumerated as in SNEG (binary counting, the first mode being the most
//   significant digit) and then stably sorted by the quantum numbers used to construct the basis of each symmetry type.

#ifndef _abelian_shell_hpp_
#define _abelian_shell_hpp_

#include <algorithm>
#include <array>
#include <bit>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <fmt/format.h>

#include "invar.hpp"

namespace NRG::AbelianShell {

using cmpl = std::complex<double>;
using State = std::uint32_t;

enum class Kind { NONE, P, U1, QSZ, SPU1 };

[[nodiscard]] inline Kind kind_from_name(const std::string &symtype) {
  if (symtype == "NONE") return Kind::NONE;
  if (symtype == "P") return Kind::P;
  if (symtype == "U1") return Kind::U1;
  if (symtype == "QSZ") return Kind::QSZ;
  if (symtype == "SPU1") return Kind::SPU1;
  throw std::invalid_argument(fmt::format("No abelian shell tables for symmetry type {}.", symtype));
}

enum class Spin : std::size_t { up = 0, down = 1 };
inline constexpr std::array<Spin, 2> both_spins{Spin::up, Spin::down};

[[nodiscard]] constexpr std::size_t mode(const std::size_t ch, const Spin s) { return 2 * ch + static_cast<std::size_t>(s); }
[[nodiscard]] constexpr bool occupied(const State s, const std::size_t m) { return ((s >> m) & 1U) != 0; }
[[nodiscard]] inline int nr_electrons(const State s) { return std::popcount(s); }
[[nodiscard]] inline int nr_electrons(const State s, const std::size_t ch, const Spin sp) { return occupied(s, mode(ch, sp)) ? 1 : 0; }
[[nodiscard]] inline int nr_electrons(const State s, const std::size_t ch) { return nr_electrons(s, ch, Spin::up) + nr_electrons(s, ch, Spin::down); }

// Elementary fermion operator c_m or c^dag_m.
struct Fermion {
  std::size_t mode;
  bool dagger;
};

// c^(dag)_m |s> = sign |s'>; sign is 0 if the result vanishes. The sign counts the occupied modes to the left of m.
[[nodiscard]] inline std::pair<int, State> apply(const Fermion &op, const State s) {
  if (occupied(s, op.mode) == op.dagger) return {0, 0};
  const auto sign = std::popcount(s & ((State{1} << op.mode) - 1)) % 2 == 0 ? 1 : -1;
  return {sign, s ^ (State{1} << op.mode)};
}

// Operator: linear combination of products of fermion operators; in each product the rightmost factor acts first.
struct Term {
  cmpl coef;
  std::vector<Fermion> ops;
};
using Operator = std::vector<Term>;

[[nodiscard]] inline cmpl matrix_element(const State bra, const Operator &op, const State ket) {
  cmpl result{};
  for (const auto &[coef, ops] : op) {
    State s = ket;
    int sign = 1;
    for (auto it = ops.rbegin(); it != ops.rend() && sign != 0; ++it) {
      const auto [sg, snew] = apply(*it, s);
      sign *= sg;
      s = snew;
    }
    if (sign != 0 && s == bra) result += coef * double(sign);
  }
  return result;
}

[[nodiscard]] inline Operator operator+(Operator a, const Operator &b) { a.insert(a.end(), b.begin(), b.end()); return a; }
[[nodiscard]] inline Operator operator*(const cmpl z, Operator a) { for (auto &t : a) t.coef *= z; return a; }

[[nodiscard]] inline Operator cdag(const std::size_t ch, const Spin s) { return {{1.0, {{mode(ch, s), true}}}}; }
[[nodiscard]] inline Operator c(const std::size_t ch, const Spin s) { return {{1.0, {{mode(ch, s), false}}}}; }
[[nodiscard]] inline Operator unity() { return {{1.0, {}}}; }
[[nodiscard]] inline Operator product(const Operator &a, const Operator &b) {
  Operator r;
  for (const auto &ta : a)
    for (const auto &tb : b) {
      Term t{ta.coef * tb.coef, ta.ops};
      t.ops.insert(t.ops.end(), tb.ops.begin(), tb.ops.end());
      r.push_back(std::move(t));
    }
  return r;
}
[[nodiscard]] inline Operator number(const std::size_t ch, const Spin s) { return product(cdag(ch, s), c(ch, s)); }
[[nodiscard]] inline Operator number(const std::size_t ch) { return number(ch, Spin::up) + number(ch, Spin::down); }

// Local operators of channel ch (SNEG conventions).
[[nodiscard]] inline Operator spinz(const std::size_t ch) { return 0.5 * (number(ch, Spin::up) + (-1.0) * number(ch, Spin::down)); }
[[nodiscard]] inline Operator spinplus(const std::size_t ch) { return product(cdag(ch, Spin::up), c(ch, Spin::down)); }
[[nodiscard]] inline Operator spinminus(const std::size_t ch) { return product(cdag(ch, Spin::down), c(ch, Spin::up)); }
[[nodiscard]] inline Operator spinx(const std::size_t ch) { return 0.5 * (spinplus(ch) + spinminus(ch)); }
[[nodiscard]] inline Operator spiny(const std::size_t ch) { return cmpl(0, -0.5) * (spinplus(ch) + (-1.0) * spinminus(ch)); }
[[nodiscard]] inline Operator isospinplus(const std::size_t ch) { return product(cdag(ch, Spin::up), cdag(ch, Spin::down)); }
[[nodiscard]] inline Operator isospinminus(const std::size_t ch) { return product(c(ch, Spin::down), c(ch, Spin::up)); }
[[nodiscard]] inline Operator isospinx(const std::size_t ch) { return 0.5 * (isospinplus(ch) + isospinminus(ch)); }
[[nodiscard]] inline Operator isospiny(const std::size_t ch) { return cmpl(0, -0.5) * (isospinplus(ch) + (-1.0) * isospinminus(ch)); }
[[nodiscard]] inline Operator isospinz(const std::size_t ch) { return 0.5 * (number(ch) + (-1.0) * unity()); }
[[nodiscard]] inline Operator charge(const std::size_t ch) { return number(ch) + (-1.0) * unity(); } // average occupancy subtracted
// Hopping between channels a and b on the added site, sum_sigma (a^dag_sigma b_sigma + h.c.)
[[nodiscard]] inline Operator hop(const std::size_t a, const std::size_t b) {
  Operator r;
  for (const auto s : both_spins) r = r + product(cdag(a, s), c(b, s)) + product(cdag(b, s), c(a, s));
  return r;
}

// Entry (i, j) with value of a table; combination indexes i, j are 1-based, as in Symmetry::make_matrix().
template<typename T> struct Entry {
  std::size_t i;
  std::size_t j;
  T value;
};

class Shell {
 public:
   Shell(const Kind kind, const std::size_t channels) : kind_(kind), channels_(channels) {
     if (channels_ == 0 || 2 * channels_ > 16) throw std::invalid_argument(fmt::format("Unsupported number of channels {}.", channels_));
     const auto n = 2 * channels_;
     const State nr_states = State{1} << n;
     for (State k = 0; k < nr_states; ++k) {
       State s = 0;
       for (std::size_t m = 0; m < n; ++m)
         if ((k >> (n - 1 - m)) & 1U) s |= State{1} << m; // first mode = most significant binary digit
       states_.push_back(s);
     }
     const auto by = [this](auto key) { std::stable_sort(states_.begin(), states_.end(), [&key](State a, State b) { return key(a) < key(b); }); };
     const auto Q = [this](State s) { return nr_electrons(s) - int(channels_); };
     const auto twoSz = [this](State s) { int r = 0; for (std::size_t ch = 0; ch < channels_; ++ch) r += nr_electrons(s, ch, Spin::up) - nr_electrons(s, ch, Spin::down); return r; };
     const auto QSz = [&](State s) { return std::pair{Q(s), twoSz(s)}; };
     switch (kind_) {
       case Kind::NONE: break;                                    // nonebasis: enumeration order
       case Kind::P: by(Q); by([&](State s) { return parity(s); }); break; // qbasis, then fermion parity
       case Kind::U1: by(Q); break;                               // qbasis
       case Kind::QSZ: by(QSz); break;                            // qszbasis
       case Kind::SPU1: by(QSz); by(twoSz); break;                // qszbasis, then Sz
     }
     for (const auto s : states_) {
       switch (kind_) {
         case Kind::NONE: qn_.emplace_back(0); in_.emplace_back(0); break;
         case Kind::P: qn_.emplace_back(parity(s)); in_.emplace_back(parity(s)); break;
         case Kind::U1: qn_.emplace_back(Q(s)); in_.emplace_back(-Q(s)); break;
         case Kind::QSZ: qn_.emplace_back(Q(s), twoSz(s)); in_.emplace_back(-Q(s), -twoSz(s)); break;
         case Kind::SPU1: qn_.emplace_back(twoSz(s)); in_.emplace_back(-twoSz(s)); break;
       }
     }
   }
   [[nodiscard]] Kind kind() const noexcept { return kind_; }
   [[nodiscard]] std::size_t channels() const noexcept { return channels_; }
   [[nodiscard]] std::size_t size() const noexcept { return states_.size(); } // number of combinations, 4^channels
   [[nodiscard]] State state(const std::size_t i) const { return states_.at(i - 1); }
   [[nodiscard]] const Invar &qn(const std::size_t i) const { return qn_.at(i - 1); } // quantum numbers of the added states
   [[nodiscard]] const Invar &in(const std::size_t i) const { return in_.at(i - 1); } // ancestor = I.combine(in)
   [[nodiscard]] int electrons(const std::size_t i) const { return nr_electrons(state(i)); }
   [[nodiscard]] int electrons(const std::size_t i, const std::size_t ch) const { return nr_electrons(state(i), ch); }
   [[nodiscard]] int electrons(const std::size_t i, const std::size_t ch, const Spin s) const { return nr_electrons(state(i), ch, s); }
   [[nodiscard]] cmpl me(const std::size_t i, const Operator &op, const std::size_t j) const { return matrix_element(state(i), op, state(j)); }

   // All nonzero matrix elements <i|op|j>, row-major; with upper=true only i <= j.
   [[nodiscard]] std::vector<Entry<cmpl>> entries(const Operator &op, const bool upper = false) const {
     std::vector<Entry<cmpl>> r;
     for (std::size_t i = 1; i <= size(); ++i)
       for (std::size_t j = upper ? i : 1; j <= size(); ++j)
         if (const auto v = me(i, op, j); v != cmpl{}) r.push_back({i, j, v});
     return r;
   }
   // Hopping term f^dag_N a_{N+1}, where f^dag_N acts on the previous sites and a on the added site:
   // <i,r| f^dag_N a |j,r'> = (-1)^{n_i} <i|a|j> <r|f^dag_N|r'>. Returns the coefficients of <r|f^dag_N|r'>.
   [[nodiscard]] std::vector<Entry<double>> hop_cr(const Operator &a, const bool upper = false) const {
     std::vector<Entry<double>> r;
     for (const auto &[i, j, v] : entries(a, upper)) r.push_back({i, j, sign(i) * v.real()});
     return r;
   }
   // Hopping term a^dag_{N+1} f_N: <i,r| a^dag f_N |j,r'> = (-1)^{n_j} <i|a^dag|j> <r|f_N|r'>.
   [[nodiscard]] std::vector<Entry<double>> hop_an(const Operator &adag, const bool upper = false) const {
     std::vector<Entry<double>> r;
     for (const auto &[i, j, v] : entries(adag, upper)) r.push_back({i, j, sign(j) * v.real()});
     return r;
   }
   // Recalculation of an operator O of the previous sites: rows (i, i, (-1)^{n_i |O|}).
   [[nodiscard]] std::vector<Entry<double>> recalc_rows(const bool fermionic) const {
     std::vector<Entry<double>> r;
     for (std::size_t i = 1; i <= size(); ++i) r.push_back({i, i, fermionic ? double(sign(i)) : 1.0});
     return r;
   }
   [[nodiscard]] int sign(const std::size_t i) const { return electrons(i) % 2 == 0 ? 1 : -1; }

 private:
   // Fermion parity label of the P symmetry type: +1 for an odd, -1 for an even charge Q = N - K measured from half
   // filling (nrgcoef convention 2(Mod[Q,2]-1/2)).
   [[nodiscard]] int parity(const State s) const { return (nr_electrons(s) - int(channels_)) % 2 != 0 ? 1 : -1; }
   Kind kind_;
   std::size_t channels_;
   std::vector<State> states_;
   std::vector<Invar> qn_, in_;
};

} // namespace NRG::AbelianShell

#endif
