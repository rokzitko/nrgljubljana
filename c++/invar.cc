// Quantum numbers
// Rok Zitko, rok.zitko@ijs.si, 2006-2019

#ifndef _nrg_invar_
#define _nrg_invar_

// Conversion functions: multiplicity (integer) -> quantum number (floating point)
inline double S(Sspin SS)    { return (SS-1.0)/2.0; }
inline double ISO(Ispin II)  { return (II-1.0)/2.0; }
inline double SZ(SZspin SSZ) { return (SSZ-1.0)/2.0; }

typedef std::vector<int> InvType;
typedef std::vector<int> QNType;
typedef map<string, int> QNNameType;

const int multiplicative = 0;
const int additive = 1;
const int mod3 = 3;

// GLOBAL VARIABLE: 'qntype' must be defined before calls to
// Invar::combine() and Invar::invert().
QNType qntype;

// Names of quantum numbers. 'name' must be defined before calls to
// Invar::get().
QNNameType name;

// invdim holds the number of quantum numbers required to specify the
// invariant subspaces (representations), i.e. (in more fancy terms)
// the dimension of the Cartan subalgebra of the full symmetry
// algebra.
size_t invdim = 0; // 0 before initialization!

class Invar {
protected:
   InvType data;

private:
   friend class boost::serialization::access;
   template<class Archive> void serialize(Archive & ar, const unsigned int version) {
      ar & data;
   }
   
public:
   Invar() : data(invdim) {}
   Invar(const InvType &d) {
      my_assert(d.size() == data.size());
      data = d;
   }
   Invar(int i0) : data { i0 } { // (int) constructor
      my_assert(invdim == 1);
   }
   Invar(int i0, int i1) : data { i0,i1 } { // (int,int) constructor
      my_assert(invdim == 2);
   }
   Invar(int i0, int i1, int i2) : data { i0,i1,i2 } { // (int,int,int) constructor
      my_assert(invdim == 3);
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
	 int qn;
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
   int getqn(size_t i) const {
      my_assert(i < data.size());
      return data[i];
   }
   // Quantum number addition laws for forming tensor products. For
   // SU(2) and U(1) (additive quantum numbers) this is simply
   // addition.
   void combine(const Invar &invar2) {
      for (size_t i = 0; i < invdim; i++)
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
	 default:
	    my_assert_not_reached();
	 }
   }
   // In DMNRG runs, we must perform the "inverse of the quantum
   // number addition", i.e. find subspaces that an invariant
   // subspaces contributed *to*. 
   void inverse() {
      for (size_t i = 0; i < invdim; i++)
	 switch (qntype[i]) {
	 case additive:
	    // For SU(2) and U(1) related quantum numbers, we just flip
	    // the sign.
	    data[i] = -data[i];
	    break;
	 case multiplicative:
	    my_assert(data[i] == 1 || data[i] == -1);
	    // do nothing since multiplication and division are
	    // equivalent! In principle, data[i] = 1/data[i].
	    break;
	 case mod3:
	    my_assert(0 <= data[i] && data[i] <= 2);
	    // Flip sign and take modulo 3 value.
	    data[i] = (3-data[i]) % 3;
	    my_assert(0 <= data[i] && data[i] <= 2);
	    break;
	 }
   }
   int get(const string &which) const {
      const auto i = name.find(which);
      if (i == end(name))
	 my_error("%s is an unknown quantum number.", which.c_str());
      const size_t index = i->second;
      my_assert(index < invdim);
      const int type = qntype[index];
      const int f = getqn(index);
      // Sanity check
      switch (type) {
      case multiplicative:
	 my_assert(f == 1 || f == -1);
	 break;
      case mod3:
	 my_assert(0 <= f && f <= 2);
	 break;
      }
      return f;
   }
   void InvertParity() {
      // By convention (for QSLR, QSZLR, ISOLR, ISOSZLR), parity is the
      // quantum number named "P"
      const auto i = name.find("P");
      if (i == end(name))
	 my_error("Critical error: no P quantum number");
      const auto index = i->second;
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
  for (size_t i = 0; i < invdim; i++) {
    name[structure[i].qn_name] = i;
    const int t = structure[i].qn_type;
    my_assert(t == additive || t == multiplicative || t == mod3);
    qntype[i] = t;
    cout << i << " " << structure[i].qn_name << " " 
	             << typestr(t) << endl;
  }
}

#endif // _nrg_invar_
