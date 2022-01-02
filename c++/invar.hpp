// Quantum numbers
// Rok Zitko, rok.zitko@ijs.si, 2006-2020

#ifndef _nrg_invar_hpp_
#define _nrg_invar_hpp_

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <stdexcept>
#include <boost/serialization/vector.hpp> // for InvType = std::vector<int>
#define FMT_HEADER_ONLY
#include <fmt/format.h>
#include "portabil.hpp"

namespace NRG {

// Conversion functions: multiplicity (integer) -> quantum number (floating point)
// WARNING: avoid using S as template variable in places where S() is used!!!
inline double S(const int SS)    { return (SS - 1.0) / 2.0;  }
inline double ISO(const int II)  { return (II - 1.0) / 2.0;  }
inline double SZ(const int SSZ)  { return (SSZ - 1.0) / 2.0; }

using InvType = std::vector<int>;

constexpr int multiplicative = 0;
constexpr int additive       = 1;
constexpr int mod3           = 3;

inline auto typestr(const int t) {
  switch (t) {
    case multiplicative: return "multiplicative";
    case additive:       return "additive";
    case mod3:           return "mod3";
    default: my_assert_not_reached();
  }
}

class Invar {
 private:
   InvType data;
   friend class boost::serialization::access;
   template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &data; }
 public:
   inline static std::vector<int> qntype;         // must be defined before calls to Invar::combine() and Invar::invert()
   inline static std::map<std::string, int> names; // must be defined before calls to Invar::get()
   // invdim holds the number of quantum numbers required to specify the invariant subspaces (representations), i.e. (in
   // more fancy terms) the dimension of the Cartan subalgebra of the full symmetry algebra.
   inline static size_t invdim = 0; // 0 before initialization!
   Invar() : data(invdim) {}
   explicit Invar(const InvType &d) {
     my_assert(d.size() == data.size());
     data = d;
   }
   explicit Invar(const int i0) : data{i0} {} // (int) constructor
   explicit Invar(const int i0, const int i1) : data{i0, i1} {} // (int,int) constructor
   explicit Invar(const int i0, const int i1, const int i2) : data{i0, i1, i2} {} // (int,int,int) constructor
   std::ostream &insertor(std::ostream &os, const std::string &delim = " "s) const {
     for (size_t i = 0; i < data.size(); i++) os << data[i] << (i != data.size() - 1 ? delim : "");
     return os;
   }
   friend std::ostream &operator<<(std::ostream &os, const Invar &invar) { return invar.insertor(os); }
   [[nodiscard]] auto str() const { std::ostringstream s; insertor(s); return s.str(); }
   [[nodiscard]] auto name() const { std::ostringstream s; insertor(s, "_"s); return s.str(); }
   auto & extractor(std::istream &is) {
     for (auto &i : data) {
       int qn{};
       if (is >> qn)
         i = qn;
       else
         throw std::runtime_error("Failed reading quantum numbers.");
     }
     return is;
   }
   friend auto &operator>>(std::istream &is, Invar &invar) { return invar.extractor(is); }
   bool operator==(const Invar &invar2) const { return data == invar2.data; }
   bool operator!=(const Invar &invar2) const { return !operator==(invar2); }
   bool operator<(const Invar &invar2)  const { return data < invar2.data; }
   // Accessor needed, because data is private.
   [[nodiscard]] auto getqn(const size_t i) const {
     my_assert(i < data.size());
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
       case mod3:
         my_assert(0 <= data[i] && data[i] <= 2);
         data[i] += invar2.getqn(i); // Modulo 3 additive
         data[i] = data[i] % 3;
         break;
       default: my_assert_not_reached();
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
       case mod3:
         my_assert(0 <= data[i] && data[i] <= 2);
         // Flip sign and take modulo 3 value.
         data[i] = (3 - data[i]) % 3;
         my_assert(0 <= data[i] && data[i] <= 2);
         break;
       }
   }
   [[nodiscard]] auto get(const std::string &which) const {
     const auto i = names.find(which);
     if (i == end(names)) throw std::invalid_argument(fmt::format("{} is an unknown quantum number.", which));
     const auto index = i->second;
     my_assert(index < invdim);
     const auto type = qntype[index];
     const auto f    = getqn(index);
     // Sanity check
     switch (type) {
     case multiplicative: my_assert(f == 1 || f == -1); break;
     case mod3: my_assert(0 <= f && f <= 2); break;
     }
     return f;
   }
   void InvertMyParity() {
     // By convention (for QSLR, QSZLR, ISOLR, ISOSZLR), parity is the quantum number named "P"
     const auto i = names.find("P");
     if (i == end(names)) throw std::invalid_argument("Critical error: no P quantum number");
     const auto index = i->second;
     data[index]      = -data[index];
   }
   [[nodiscard]] auto InvertParity() const {
     Invar I(data);
     I.InvertMyParity();
     return I;
   }
};

struct InvarStructure {
  std::string name;
  int type;
};

inline void initInvar(std::initializer_list<InvarStructure> l) {
  Invar::invdim = l.size();
  auto i = 0;
  for (const auto & [n, t]: l) {
    my_assert(t == additive || t == multiplicative || t == mod3);
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
