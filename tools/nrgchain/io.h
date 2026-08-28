// Discretization ODE solver for NRG
//
// ** Input/output code
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: io.h,v 1.1 2009/03/20 09:53:41 rok Exp $

#include <cstdlib>
#include <exception>
#include <fstream>
#include <iostream>
#include <ostream>
#include <string>

#include "../common/io.hpp"

inline double atof(const std::string &s) { return std::atof(s.c_str()); }

inline int atoi(const std::string &s) { return std::atoi(s.c_str()); }

void safe_open(std::ifstream &F, const std::string &filename) {
  try {
    NRG::Tools::open_input(F, filename);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::exit(1);
  }
}

const int PREC = 16;

void safe_open(std::ofstream &F, const std::string &filename) {
  try {
    NRG::Tools::open_output(F, filename, PREC);
  } catch (const std::exception &e) {
    std::cerr << e.what() << std::endl;
    std::exit(1);
  }
}

// Get next line from stream F, skipping empty lines and comments.
std::string getnextline(std::ifstream &F) {
  return NRG::Tools::next_data_line(F);
}
