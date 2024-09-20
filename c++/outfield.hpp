// outfield.cc - Code for flexible formatted output
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _outfield_hpp_
#define _outfield_hpp_

#include <vector>
#include <string>
#include <fstream>
#include <optional>
#include <stdexcept>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include "params.hpp"

namespace NRG {

class Outfield {
 private:
   std::string desc;
   std::string value{};
   int prec, width;
 public:
   explicit Outfield(const std::string &desc, const unsigned int prec, const unsigned int width) : desc(desc), prec(prec), width(width) {}
   template<typename T> void set_value(const T x) {
     value = fmt::format("{:>{}.{}}", x, prec, width); // https://fmt.dev/latest/syntax.html
   }
   void put_header(std::ostream &F) const { F << std::setw(width) << desc << " "; }
   void put_value(std::ostream &F) const  { F << std::setw(width) << value << " "; }
   [[nodiscard]] auto get_desc() const { return desc; }
};

class Allfields : std::vector<Outfield> {
  private:
    unsigned int prec, width;
  public:
    void add(const std::string &desc, const int position = -1) {
      if (position == -1) {
        this->emplace_back(desc, prec, width);
      } else {
        my_assert(position < this->size());
        this->insert(this->begin() + position, Outfield(desc, prec, width));
      }
    }
    void add(const Allfields &more, const int position = -1) {
      if (position != -1) my_assert(position < this->size());
      this->insert(position == -1 ? this->end() : this->begin() + position, more.begin(), more.end());
    }
    Allfields(const std::vector<std::string> &fields, const size_t prec, const size_t width) : prec(prec), width(width) {
      for(const auto &desc : fields) add(desc);
    }
    void add(const std::vector<std::string> &more, const int position = -1) {
      Allfields more_fields(more, prec, width);
      this->add(more_fields, position);
    }
    void save_header(std::ostream &O) const {
      O << '#';
      for (const auto &f : *this) f.put_header(O);
      O << std::endl;
    }
    void save_values(std::ostream &O) const {
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

class TD_generic {
 private:
  std::string filename;
  std::optional<std::ofstream> O;
 public:
  Allfields allfields;
  TD_generic(const Params &P, const std::string &filename, const std::vector<std::string> &fields) : filename(filename), allfields(fields, P.prec_td, P.width_td) {}
  template<typename T> void set(const std::string &desc, const T &x) { allfields.set(desc,x); }
  void save_values() {
    if (!O) {
      O.emplace(filename);
      allfields.save_header(O.value());
    }
    allfields.save_values(O.value());
  }
};

// Setup output fields that will appear in the file "td". Additional elements are defined in sym* files.
class TD : public TD_generic {
 public:
  TD(const Params &P, const std::string &filename) : TD_generic(P, filename, {"T", "<E>", "<E^2>", "C", "F", "S"}) {}
};

class TD_FDM : public TD_generic {
 public:
  TD_FDM(const Params &P, const std::string &filename) : TD_generic(P, filename, {"T", "E_fdm", "C_fdm", "F_fdm", "S_fdm"}) {}
};

} // namespace

#endif
