// Discretization ODE solver for NRG
// ** Parsing of the parameter file

#ifndef _adapt_parser_hpp_
#define _adapt_parser_hpp_

#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>

#include "../common/parser.hpp"

namespace NRG::Adapt {

class Params : public std::map<std::string, std::string> {
 private:
 public:
   Params(const std::string &filename) {
      std::ifstream F;
      safe_open(F, filename);
      if (NRG::Tools::find_block(F, "param")) { NRG::Tools::parse_key_value_block(F, *this); }
    }
    // Return a parameter of type double, use default value if not found.
    auto P(const std::string &keyword, const double def) const {
      return NRG::Tools::get_or_default(*this, keyword, def, [](const auto &value) { return atof(value); });
    }
    auto Pint(const std::string &keyword, const int def) const {
      return NRG::Tools::get_or_default(*this, keyword, def, [](const auto &value) { return atoi(value); });
    }
    auto Pstr(const std::string &keyword, const std::string &def) const {
      return NRG::Tools::get_or_default(*this, keyword, def, [](const auto &value) { return value; });
    }
    auto Pbool(const std::string &keyword, const bool def) const {
      return NRG::Tools::get_or_default(*this, keyword, def, [](const auto &value) { return value == "true"; });
    }
};

} // namespace

#endif
