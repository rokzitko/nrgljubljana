// This file is part of "NRG Ljubljana".
//   
// NRG Ljubljana is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//	
// NRG Ljubljana is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//			
// You should have received a copy of the GNU General Public License along
// with NRG Ljubljana; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

// Definition of the class for invariant subspaces.
// (NRG for pluggable symmetry types)
//
// Rok Zitko, rok.zitko@ijs.si, 2006-2019

#ifndef _nrg_invar_
#define _nrg_invar_

#include "rational.hpp"

using boost::rational;
using boost::rational_cast;

typedef rational<int> frac; // a fraction

// Useful constants
#define HALF frac(1, 2)
#define ONE frac(1)

// Conversion functions: multiplicity -> quantum number
inline double S(Sspin SS)    { return (SS-1.0)/2.0; }
inline double ISO(Ispin II)  { return (II-1.0)/2.0; }
inline double SZ(SZspin SSZ) { return (SSZ-1.0)/2.0; }

// For integers do not append "/1".
ostream& operator<<(ostream &os, const frac &f)
{
  os << f.numerator();
  int d = f.denominator();
  if (d != 1)
    os << '/' << d;
  return os;
}

// Allow reading both fractions and integers (i.e. fractions with
// implicit denominator 1).
istream& operator>>(istream &is, frac &f)
{
  int n, d;
  char c = 0;
  boost::detail::resetter sentry(is);

  is >> n;
  c = is.peek();
  if (c != '/') {
    f.assign(n, 1);
  } else {
    c = is.get();
    is >> std::noskipws;
    is >> d;
    if (is)
      f.assign(n, d);
  }
  return is;
}

int to_int(const frac f)
{
   const int d = f.denominator();
   if (d != 1) {
      cout << "to_int() f=" << f << endl;
      my_error("int QN expected, fractional QN found, denom=%i", d);
   }
   return f.numerator();
}

typedef std::vector<frac> InvType;
typedef std::vector<int> QNType;
typedef map<string, int> QNNameType;

const int multiplicative = 0;
const int additive = 1;
const int mod3 = 3;

// 0 = multiplicative QN, 1 = additive QN, 3 = modulo 3 additive
// 'qntype' must be defined before calls to Invar::combine()
// and Invar::invert().
QNType qntype;

// Names of quantum numbers.
// 'name' must be define before calls to Invar::get().
QNNameType name;
std::vector<string> namelist;

// For the pluggable-symmetry code (NRG_ANY), invdim holds the number
// of quantum numbers required to specify the invariant subspaces
// (representations), i.e. (in more fancy terms) the dimension of the
// Cartan subalgebra of the full symmetry algebra.

size_t invdim = 0; // 0 before initialization!

// Bug fix for serialization of STL containers on some platforms (Mac
// OS X, for instance): explicit len value. (Alternatively, we should
// propagate invdim to other MPI processes.)

class Invar {
protected:
   size_t len; // should be equal to invdim in MPI master process
   InvType data;

   friend class boost::serialization::access;
   template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & len;
      data.resize(len);
      for (size_t i = 0; i < len; i++)
	 ar & data[i];
   }
   
public:
   Invar() { 
      len = invdim;
      data = InvType(len);
   }
   
   Invar(const InvType &d) {
      my_assert(d.size() == data.size());
      data = d;
      my_assert(data.size() == len);
   }
   
   Invar(int i0) { // (int) constructor
      my_assert(invdim == 1);
      len = 1;
      data = InvType(len);
      data[0] = i0;
   }
   
   Invar(int i0, int i1) { // (int,int) constructor
      my_assert(invdim == 2);
      len = 2;
      data = InvType(len);
      data[0] = i0; data[1] = i1;
   }
   
  Invar(int i0, int i1, int i2) { // (int,int,int) constructor
     my_assert(invdim == 3);
     len = 3;
     data = InvType(len);
     data[0] = i0; data[1] = i1; data[2] = i2;
  }
	
   Invar(frac f0, frac f1) { // (frac, frac) constructor
      my_assert(invdim == 2);
      len = 2;
      data = InvType(len);
      data[0] = f0; data[1] = f1;
   }

   ostream& insertor(ostream &os) const {
      for (size_t i = 0; i < data.size(); i++)
	 os << data[i] << (i != data.size()-1 ? " " : "");
      return os;
   }

   friend ostream& operator<<(ostream &os, const Invar &invar) {
      return invar.insertor(os);
   }

   istream& extractor(istream &is) {
      for (auto & i : data) {
	 frac qn;
	 if (is >> qn)
	    i = qn;
	 else
	    my_error("Failed reading quantum numbers.");
      }
      return is;
   }

  friend istream& operator>>(istream &is, Invar &invar) {
    return invar.extractor(is);
  }

  bool operator==(const Invar &invar2) const {
    return data == invar2.data;
  }

  bool operator!=(const Invar &invar2) const {
    return !operator==(invar2);
  }

  bool operator<(const Invar &invar2) const {
    return data < invar2.data;
  }  

  // Accessor needed, because data is private.
  frac getqn(size_t i) const {
     my_assert(i < len);
     my_assert(data.size() == len);
     return data[i];
  }

  // Quantum number addition laws for forming tensor products. For
  // SU(2) and U(1) (additive quantum numbers) this is simply
  // addition.
  void combine(const Invar &invar2) {
    my_assert(qntype.size() == len);
    for (size_t i = 0; i < len; i++)
	switch (qntype[i]) {
	case additive:
	   data[i] += invar2.getqn(i); // Additive
	   break;
	case multiplicative:
	   my_assert(to_int(data[i]) == 1 || to_int(data[i]) == -1);
	   data[i] *= invar2.getqn(i); // Multiplicative
	   break;
	case mod3:
	   my_assert(0 <= data[i] && data[i] <= 2);
	   data[i] += invar2.getqn(i); // Modulo 3 additive
	   data[i] = to_int(data[i]) % 3;
	   break;
	default:
	   my_assert_not_reached();
	}
  }

  // In DMNRG runs, we must perform the "inverse of the quantum number
  // addition", i.e. find subspaces that an invariant subspaces
  // contributed *to*. 
  void inverse() {
    my_assert(qntype.size() == len);
    for (size_t i = 0; i < len; i++)
      switch (qntype[i]) {
      case additive:
	 // For SU(2) and U(1) related quantum numbers, we just flip
	 // the sign.
	 data[i] = -data[i];
	 break;
      case multiplicative:
	 my_assert(to_int(data[i]) == 1 || to_int(data[i]) == -1);
	 // do nothing since multiplication and division are
	 // equivalent! In principle, data[i] = 1/data[i].
	 break;
      case mod3:
	 my_assert(0 <= data[i] && data[i] <= 2);
	 // Flip sign and take modulo 3 value.
	 data[i] = to_int(3-data[i]) % 3;
	 my_assert(0 <= data[i] && data[i] <= 2);
	 break;
      }
  }
   
  frac get_frac(const string &which) const {
     const auto i = name.find(which);
     if (i == end(name))
	my_error("%s is an unknown quantum number.", which.c_str());
     const size_t index = i->second;
     my_assert(index < len);
     const int type = qntype[index];
     const frac f = getqn(index);
     // Sanity check
     switch (type) {
     case multiplicative:
	my_assert(to_int(f) == 1 || to_int(f) == -1);
	break;
     case mod3:
	my_assert(0 <= to_int(f) && to_int(f) <= 2);
	break;
     }
     return f;
  }
   
  // Returns an int quantum number!
  int get(const string &which) const {
     return to_int(get_frac(which));
  }

   void InvertParity() {
      // By convention (for QSLR, QSZLR, ISOLR, ISOSZLR), parity is the
      // quantum number named "P"
      const auto i = name.find("P");
      if (i == end(name))
	 my_error("Critical error: no P quantum number");
      const auto index = i->second;
//      my_assert(index < len);
      data[index] = -data[index];
   }
      
};

struct InvarStructure
{
  string qn_name;
  int qn_type;
};

std::string typestr(int t)
{
   switch (t) {
   case additive:
      return "additive";
   case multiplicative:
      return "multiplicative";
   case mod3:
      return "mod3";
   default:
      my_assert_not_reached();
   }
}

void initInvar(const InvarStructure *structure, size_t len)
{
  invdim = len; // global variable!
  cout << "invdim=" << invdim << endl;
  qntype.resize(invdim);
  namelist.resize(invdim);
  for (size_t i = 0; i < invdim; i++) {
    name[structure[i].qn_name] = i;
    namelist[i] = structure[i].qn_name;
    const int t = structure[i].qn_type;
    my_assert(t == additive || t == multiplicative || t == mod3);
    qntype[i] = t;
    cout << i << " " << structure[i].qn_name << " " 
	             << typestr(t) << endl;
  }
}

#endif // _nrg_invar_
