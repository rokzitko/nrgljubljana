// outfield.cc - Code for flexible formatted output
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _outfield_hpp_
#define _outfield_hpp_

#include <vector>
#include <string>
#include <fstream>
#include "params.hpp"

namespace NRG {

class Outfield {
 private:
   std::string desc;
   std::string value{};
   size_t prec, width;
 public:
   explicit Outfield(const std::string &desc, const size_t prec, const size_t width) : desc(desc), prec(prec), width(width) {}
   template<typename T> void set_value(const T x) {
     value = fmt::format("{x:>{width}.{prec}}", "x"_a = x, "prec"_a = prec, "width"_a = width);
   }
  // Outfield& operator=(const Outfield &o) = default;
  // Outfield& operator=(Outfield &&o) = default;
  // template<typename T> Outfield& operator=(const T x) { set_value(x); return *this; } // required for elements in std::vector
   void put_header(std::ostream &F) const { F << std::setw(width) << desc << " "; }
   void put_value(std::ostream &F) const  { F << std::setw(width) << value << " "; }
   [[nodiscard]] auto get_desc() const { return desc; }
};

class Allfields : std::vector<Outfield> {
  private:
    const size_t prec, width;
  public:
    void add(const std::string &desc, const int position = -1) {
      if (position == -1) {
        this->emplace_back(desc, prec, width);
      } else {
        my_assert(position < this->size());
        this->insert(this->begin() + position, Outfield(desc, prec, width));
      }
    }
    Allfields(const std::vector<std::string> &fields, const size_t prec, const size_t width) : prec(prec), width(width) {
      for(const auto &desc : fields) add(desc);
    }
    void save_header(std::ofstream &O) const {
      O << '#';
      for (const auto &f : *this) f.put_header(O);
      O << std::endl;
    }
    void save_values(std::ofstream &O) const {
      O << ' ';
      for (const auto &f : *this) f.put_value(O);
      O << std::endl;
    }
    template<typename T> void set(const std::string &desc, const T &x) {
      for(auto &f : *this) 
        if (f.get_desc() == desc) {
          f.set_value(x);
          return;
        }
      throw std::runtime_error("Invalid field " + desc);
    }
};

// Setup output fields that will appear in the file "td". Additional elements are defined in sym* files.
struct TD {
  const Params &P;
  Allfields allfields; // thermodynamic quantities
  const std::string filename;
  std::ofstream O; // XXX: std:optional?
  bool header_saved = false; // XXX
  void save_values() {
    if (!O.is_open())
      O.open(filename);
    if (!header_saved) {
      allfields.save_header(O);
      header_saved = true;
    }
  }
  TD(const Params &P_, const std::string &filename) 
    : P(P_), allfields({"T", "<E>", "<E^2>", "C", "F", "S"}, P.prec_td, P.width_td), filename(filename) {}
};

struct TD_FDM {
  const Params &P;
  Allfields allfields;
  std::string filename;
  std::ofstream O;
  bool header_saved = false;
  void save_values() {
    if (!O.is_open())
      O.open(filename);
    if (!header_saved) {
      allfields.save_header(O);
      header_saved = true;
    }
  }
  TD_FDM(const Params &P_, const std::string &filename)
    : P(P_), allfields({"T", "E_fdm", "C_fdm", "F_fdm", "S_fdm"}, P.prec_td, P.width_td), filename(filename) {}
};

} // namespace

#endif
