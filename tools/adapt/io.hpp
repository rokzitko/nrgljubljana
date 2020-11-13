// Discretization ODE solver for NRG
// ** Input/output code

#ifndef _adapt_io_hpp_
#define _adapt_io_hpp_

#include <string>
#include <cstdlib>
#include <fstream>
#include <string>
#include <stdexcept>

namespace NRG::Adapt {

inline double atof(const std::string &s) { return std::atof(s.c_str()); }
inline int atoi(const std::string &s) { return std::atoi(s.c_str()); }

inline void safe_open(std::ifstream &F, const std::string &filename) {
  F.open(filename.c_str());
  if (!F)
    throw std::runtime_error("Can't open " + filename + " for reading.");
}

inline void safe_open(std::ofstream &F, const std::string &filename, const int PREC = 16) {
  F.open(filename.c_str());
  if (!F)
    throw std::runtime_error("Can't open " + filename + " for writing.");
  F << std::setprecision(PREC);
}

// Get next line from stream F, skipping empty lines and comments.
inline std::string getnextline(std::ifstream &F) {
  std::string line;
  while (F) {
    std::getline(F, line);
    if (!F) // bail out
      break;
    if (line.length() == 0) // skip empty lines
      continue;
    if (line[0] == '#') // skip comment lines
      continue;
    return line;
  }
  return ""; // error
}

} // namespace

#endif
