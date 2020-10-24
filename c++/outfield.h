// outfield.cc - Code for flexible formatted output
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _outfield_h_
#define _outfield_h_

#include <vector>
#include <string>
#include <fstream>
#include "params.h"

class outfield;
using Allfields = std::vector<outfield*>;

class outfield {
 private:
   const Params &P;
   const std::string desc;
   std::string value{};
 public:
   explicit outfield(const Params &P_, Allfields &allfields, const std::string &desc_, int position = -1)
     : P(P_), desc(desc_) {
       if (position == -1)
         allfields.push_back(this);
       else {
         my_assert(position < allfields.size());
         allfields.insert(begin(allfields) + position, this);
       }
     }
   template<typename T> void setvalue(T x) {
     value = fmt::format("{x:>{width}.{prec}}", "x"_a=x, "prec"_a=P.prec_td, "width"_a=P.width_td);
   }
   template<typename T> void operator=(T x) { setvalue(x); }
   void putheader(std::ostream &F) const { F << std::setw(P.width_td) << desc << " "; }
   void putvalue(std::ostream &F) const  { F << std::setw(P.width_td) << value << " "; }
};

// Setup output fields that will appear in the file "td". Additional elements are defined in sym* files.
struct TD {
  Allfields allfields; // thermodynamic quantities
  const Params &P;
  const std::string filename;
  std::ofstream O;
  outfield T, E, E2, C, F, S;
  bool header_saved = false;
  void save_header() {
    O << '#';
    for (const auto &i : allfields) i->putheader(O);
    O << std::endl;
  }
  void save_values() {
    if (!O.is_open())
      O.open(filename);
    if (!header_saved) {
      save_header();
      header_saved = true;
    }
    O << ' ';
    for (const auto &i : allfields) i->putvalue(O);
    O << std::endl;
  }
  TD(const Params &P_, const std::string &filename) : P(P_), filename(filename),
    T(P, allfields, "T"), E(P, allfields, "<E>"), E2(P, allfields, "<E^2>"),
    C(P, allfields, "C"), F(P, allfields, "F"), S(P, allfields, "S") {}
};

struct TD_FDM {
  Allfields allfields;
  const Params &P;
  std::string filename;
  std::ofstream O;
  outfield T, E, C, F, S;
  bool header_saved = false;
  void save_header() {
    O << '#';
    for (const auto &i : allfields) i->putheader(O);
    O << std::endl;
  }
  void save_values() {
    if (!O.is_open())
      O.open(filename);
    if (!header_saved) {
      save_header();
      header_saved = true;
    }
    O << ' ';
    for (const auto &i : allfields) i->putvalue(O);
    O << std::endl;
  }
  TD_FDM(const Params &P_, const std::string &filename) : P(P_), filename(filename),
    T(P, allfields, "T"), E(P, allfields, "E_fdm"),
    C(P, allfields, "C_fdm"), F(P, allfields, "F_fdm"), S(P, allfields, "S_fdm") {}
};

#endif
