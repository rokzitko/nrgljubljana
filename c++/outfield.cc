// outfield.cc - Code for flexible formatted output
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _outfield_cc_
#define _outfield_cc_

class outfield;
using Allfields = std::vector<outfield*>;

class outfield {
 private:
   const Params &P;
   string desc;
   string value{};
 public:
   explicit outfield(const Params &P_, Allfields &allfields, const string &desc_, int position = -1)
     : P(P_), desc(desc_) {
       if (position == -1)
         allfields.push_back(this);
       else {
         my_assert(position < allfields.size());
         allfields.insert(begin(allfields) + position, this);
       }
     }
   template<typename T> void setvalue(T x) {
     ostringstream tmp;
     tmp << setw(P.width_td) << setprecision(P.prec_td) << x; // XXX fmt
     value    = tmp.str();
   }
   template<typename T> void operator=(T x) { setvalue(x); }
   void putheader(ostream &F) const { F << setw(P.width_td) << desc << " "; }
   void putvalue(ostream &F) const  { F << setw(P.width_td) << value << " "; }
};

// Setup output fields that will appear in the file "td".
// Additional elements are defined in symmetry.cc.
struct TD {
  Allfields allfields;
  const Params &P;
  ofstream O;
  outfield T, E, E2, C, F, S;
  void save_values() {
    O << ' ';
    for (const auto &i : allfields) i->putvalue(O);
    O << std::endl;
  }
  void save_header() {
    O << '#';
    for (const auto &i : allfields) i->putheader(O);
    O << std::endl;
  }
  TD(const Params &P_, std::string filename) : P(P_), O(filename), 
    T(P, allfields, "T"), E(P, allfields, "<E>"), E2(P, allfields, "<E^2>"),
    C(P, allfields, "C"), F(P, allfields, "F"), S(P, allfields, "S") { 
      save_header(); 
    }
};

struct TD_FDM {
  Allfields allfields;
  const Params &P;
  ofstream O;
  outfield T, E, C, F, S;
  void save_values() {
    O << ' ';
    for (const auto &i : allfields) i->putvalue(O);
    O << std::endl;
  }
  void save_header() {
    O << '#';
    for (const auto &i : allfields) i->putheader(O);
    O << std::endl;
  }
  TD_FDM(const Params &P_, std::string filename) : P(P_), O(filename), 
    T(P, allfields, "T"), E(P, allfields, "E_fdm"),
    C(P, allfields, "C_fdm"), F(P, allfields, "F_fdm"), S(P, allfields, "S_fdm") { 
      save_header(); 
    }
};

#endif // _outfield_cc_
