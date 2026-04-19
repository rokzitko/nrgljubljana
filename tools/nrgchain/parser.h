// Discretization ODE solver for NRG
//
// ** Parsing of the parameter file
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: parser.h,v 1.1 2009/03/20 09:53:41 rok Exp $

#include "../common/parser.hpp"

map<string, string> params;

// Return a parameter of type double, use default value if not found.
double P(const string &keyword, double def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return atof(value); });
}

int Pint(const string &keyword, int def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return atoi(value); });
}

string Pstr(const string &keyword, string def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return value; });
}

bool Pbool(const string &keyword, bool def) {
  return NRG::Tools::get_or_default(params, keyword, def, [](const auto &value) { return value == "true"; });
}

void parser(const string &filename) {
  ifstream F;
  safe_open(F, filename);

  if (NRG::Tools::find_block(F, "param")) { NRG::Tools::parse_key_value_block(F, params); }
}
