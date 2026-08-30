// Discretization ODE solver for NRG
//
// ** Input/output code
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de

#include <cstdlib>
#include <fstream>
#include <string>

#include "../common/io.hpp"

inline double atof(const std::string &s) { return std::atof(s.c_str()); }

inline int atoi(const std::string &s) { return std::atoi(s.c_str()); }

void safe_open(std::ifstream &F, const std::string &filename) {
  NRG::Tools::open_input(F, filename);
}

const int PREC = 16;

void safe_open(std::ofstream &F, const std::string &filename) {
  NRG::Tools::open_output(F, filename, PREC);
}

// Get next line from stream F, skipping empty lines and comments.
std::string getnextline(std::ifstream &F) {
  return NRG::Tools::next_data_line(F);
}
