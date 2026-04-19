// Discretization ODE solver for NRG
// ** Input/output code

#ifndef _adapt_io_hpp_
#define _adapt_io_hpp_

#include <string>
#include <cstdlib>
#include <fstream>
#include <string>
#include <stdexcept>

#include "../common/io.hpp"

namespace NRG::Adapt {

inline double atof(const std::string &s) { return std::atof(s.c_str()); }
inline int atoi(const std::string &s) { return std::atoi(s.c_str()); }

inline void safe_open(std::ifstream &F, const std::string &filename) {
  NRG::Tools::open_input(F, filename);
}

inline void safe_open(std::ofstream &F, const std::string &filename, const int PREC = 16) {
  NRG::Tools::open_output(F, filename, PREC);
}

// Get next line from stream F, skipping empty lines and comments.
inline std::string getnextline(std::ifstream &F) {
  return NRG::Tools::next_data_line(F);
}

} // namespace

#endif
