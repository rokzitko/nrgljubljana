// Quantum numbers
// Rok Zitko, rok.zitko@ijs.si, 2006-2020

#ifndef _nrg_invar_hpp_
#define _nrg_invar_hpp_

#include <algorithm>
#include <array>
#include <cstddef>
#include <functional>
#include <initializer_list>
#include <ios>
#include <istream>
#include <iterator>
#include <ostream>
#include <sstream>
#include <utility>
#include <string>
#include <string_view>
#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/vector.hpp> // serialized form of Invar

#include <fmt/format.h>

#include "portabil.hpp"

namespace NRG {

using std::string_literals::operator""s;

// Conversion functions: multiplicity (integer) -> quantum number (floating point)
// WARNING: avoid using S as template variable in places where S() is used!!!
inline double S(const int SS)    { return (SS - 1.0) / 2.0;  }
inline double ISO(const int II)  { return (II - 1.0) / 2.0;  }
inline double SZ(const int SSZ)  { return (SSZ - 1.0) / 2.0; }

using InvType = std::vector<int>;

// Combination rules of the quantum numbers. A type n >= 2 is addition modulo n, with values 0..n-1 (Z_n symmetry).
constexpr int multiplicative = 0;
constexpr int additive       = 1;
constexpr int mod3           = 3; // Z_3, e.g. C_3 rotations of three channels
constexpr int mod4           = 4; // Z_4, e.g. 2S_z modulo 4 for rhombic spin anisotropy

[[nodiscard]] constexpr bool is_modulo(const int t) { return t >= 2; }

inline std::string typestr(const int t) {
  if (t == multiplicative) return "multiplicative";
  if (t == additive) return "additive";
  if (is_modulo(t)) return "mod" + std::to_string(t);
  my_assert_not_reached();
}

// Quantum numbers of an invariant subspace. At most MAX_INVDIM numbers are stored inline (the symmetry types use up to
// three), so that copying and combining Invar objects, which happens in the innermost loops of the recalculation, does
// not allocate memory. Ordering is lexicographic, as for std::vector<int>, and the serialized form is that of a
// std::vector<int>, so checkpoint and MPI archives do not depend on the in-memory storage.
class Invar {
 public:
   static constexpr size_t MAX_INVDIM = 4;
 private:
   std::array<int, MAX_INVDIM> data{};
   size_t len = 0;
   friend class boost::serialization::access;
   template <class Archive> void save(Archive &ar, [[maybe_unused]] const unsigned int version) const {
     const InvType v(data.begin(), data.begin() + len);
     ar &v;
   }
   template <class Archive> void load(Archive &ar, [[maybe_unused]] const unsigned int version) {
     InvType v;
     ar &v;
     assign(v);
   }
   BOOST_SERIALIZATION_SPLIT_MEMBER()
   void assign(const InvType &v) {
     if (v.size() > MAX_INVDIM) throw std::invalid_argument("Too many quantum numbers in Invar.");
     len = v.size();
     std::copy(v.begin(), v.end(), data.begin());
   }
   [[nodiscard]] auto qn_begin() const { return data.begin(); }
   [[nodiscard]] auto qn_end() const { return data.begin() + len; }
 public:
   inline static std::vector<int> qntype;         // must be defined before calls to Invar::combine() and Invar::invert()
   inline static std::map<std::string, int, std::less<>> names; // must be defined before calls to Invar::get()
   // invdim holds the number of quantum numbers required to specify the invariant subspaces (representations), i.e. (in
   // more fancy terms) the dimension of the Cartan subalgebra of the full symmetry algebra.
   inline static size_t invdim = 0; // 0 before initialization!
   Invar() : len(invdim) { my_assert(invdim <= MAX_INVDIM); }
   explicit Invar(const InvType &d) {
     assign(d);
     my_assert(len == invdim);
   }
   explicit Invar(const int i0) : data{i0}, len(1) {} // (int) constructor
   explicit Invar(const int i0, const int i1) : data{i0, i1}, len(2) {} // (int,int) constructor
   explicit Invar(const int i0, const int i1, const int i2) : data{i0, i1, i2}, len(3) {} // (int,int,int) constructor
   std::ostream &insertor(std::ostream &os, const std::string &delim = " "s) const {
     for (size_t i = 0; i < len; i++) os << data[i] << (i != len - 1 ? delim : "");
     return os;
   }
   friend std::ostream &operator<<(std::ostream &os, const Invar &invar) { return invar.insertor(os); }
   [[nodiscard]] auto str() const { std::ostringstream s; insertor(s); return s.str(); }
   [[nodiscard]] auto name() const { std::ostringstream s; insertor(s, "_"s); return s.str(); }
   auto & extractor(std::istream &is) {
     for (size_t i = 0; i < len; i++) {
       int qn{};
       if (is >> qn)
         data[i] = qn;
       else
         throw std::runtime_error("Failed reading quantum numbers.");
     }
     return is;
   }
   friend auto &operator>>(std::istream &is, Invar &invar) { return invar.extractor(is); }
   bool operator==(const Invar &invar2) const { return std::equal(qn_begin(), qn_end(), invar2.qn_begin(), invar2.qn_end()); }
   bool operator!=(const Invar &invar2) const { return !operator==(invar2); }
   bool operator<(const Invar &invar2)  const { return std::lexicographical_compare(qn_begin(), qn_end(), invar2.qn_begin(), invar2.qn_end()); }
   // Accessor needed, because data is private.
   [[nodiscard]] auto getqn(const size_t i) const {
     my_assert(i < len);
     return data[i];
   }
   // Quantum number addition laws for forming tensor products. For SU(2) and U(1) (additive quantum numbers) this is
   // simply addition.
   void combine(const Invar &invar2) {
     for (size_t i = 0; i < invdim; i++) {
       switch (qntype[i]) {
       case additive:
         data[i] += invar2.getqn(i); // Additive
         break;
       case multiplicative:
         my_assert(data[i] == 1 || data[i] == -1);
         data[i] *= invar2.getqn(i); // Multiplicative
         break;
       default: {
         const auto n = qntype[i];
         my_assert(is_modulo(n) && 0 <= data[i] && data[i] < n);
         data[i] = ((data[i] + invar2.getqn(i)) % n + n) % n; // addition modulo n
       }
       }
     }
   }
   // In DMNRG runs, we must perform the "inverse of the quantum number addition", i.e. find subspaces that an
   // invariant subspaces contributed *to*.
   void inverse() {
     for (size_t i = 0; i < invdim; i++)
       switch (qntype[i]) {
       case additive:
         // For SU(2) and U(1) related quantum numbers, we just flip the sign.
         data[i] = -data[i];
         break;
       case multiplicative:
         my_assert(data[i] == 1 || data[i] == -1);
         // do nothing since multiplication and division are equivalent! In principle, data[i] = 1/data[i].
         break;
       default: {
         const auto n = qntype[i];
         my_assert(is_modulo(n) && 0 <= data[i] && data[i] < n);
         data[i] = (n - data[i]) % n; // flip sign modulo n
       }
       }
   }
   [[nodiscard]] auto get(const std::string_view which) const {
     const auto i = names.find(which);
     if (i == end(names)) throw std::invalid_argument(fmt::format("{} is an unknown quantum number.", which));
     const auto index = static_cast<size_t>(i->second);
     my_assert(index < invdim);
     const auto type = qntype[index];
     const auto f    = getqn(index);
     // Sanity check
     if (type == multiplicative) my_assert(f == 1 || f == -1);
     if (is_modulo(type)) my_assert(0 <= f && f < type);
     return f;
   }
   void InvertMyParity() {
     // By convention (for QSLR, QSZLR, ISOLR, ISOSZLR), parity is the quantum number named "P"
     const auto i = names.find("P");
     if (i == end(names)) throw std::invalid_argument("Critical error: no P quantum number");
     const auto index = static_cast<size_t>(i->second);
     my_assert(index < len);
     data[index]      = -data[index];
   }
   [[nodiscard]] auto InvertParity() const {
     Invar I(*this);
     I.InvertMyParity();
     return I;
   }
};

struct InvarStructure {
  std::string name;
  int type;
};

inline void initInvar(std::initializer_list<InvarStructure> l) {
  my_assert(l.size() <= Invar::MAX_INVDIM);
  Invar::invdim = l.size();
  Invar::qntype.clear();
  Invar::qntype.reserve(Invar::invdim);
  Invar::names.clear();
  auto i = 0;
  for (const auto & [n, t]: l) {
    my_assert(t == additive || t == multiplicative || is_modulo(t));
    Invar::names[n] = i++;
    Invar::qntype.push_back(t);
    std::cout << fmt::format("{} {} {}\n", i, n, typestr(t));
  }
}

using InvarVec = std::vector<Invar>; // holds information about ancestor subspaces
using Twoinvar = std::pair<Invar, Invar>; // labels subspace bra and subspace ket for matrix elements

inline std::ostream &operator<<(std::ostream &os, const Twoinvar &p) { return os << "(" << p.first << ") (" << p.second << ")"; }

} // namespace

#endif
