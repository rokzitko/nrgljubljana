// Discretization ODE solver for NRG
//
// ** Parsing of the parameter file
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de

#include <fstream>
#include <map>
#include <string>

#include <parse_bool.hpp>

#include "../common/parser.hpp"

std::map<std::string, std::string> params;

// Return a parameter of type double, use default value if not found.
double P(const std::string &keyword, double def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return atof(value); });
}

int Pint(const std::string &keyword, int def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return atoi(value); });
}

std::string Pstr(const std::string &keyword, std::string def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return value; });
}

bool Pbool(const std::string &keyword, bool def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return NRG::parse_bool(value); });
}

void parser(const std::string &filename) {
  std::ifstream F;
  safe_open(F, filename);

  if (NRG::Tools::find_block(F, "param")) { NRG::Tools::parse_key_value_block(F, params); }
}
