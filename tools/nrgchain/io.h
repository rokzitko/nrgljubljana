// Discretization ODE solver for NRG
//
// ** Input/output code
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: io.h,v 1.1 2009/03/20 09:53:41 rok Exp $

#include "../common/io.hpp"

inline double atof(const string &s) { return atof(s.c_str()); }

inline int atoi(const string &s) { return atoi(s.c_str()); }

void safe_open(ifstream &F, const string &filename) {
  try {
    NRG::Tools::open_input(F, filename);
  } catch (const std::exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
}

const int PREC = 16;

void safe_open(ofstream &F, const string &filename) {
  try {
    NRG::Tools::open_output(F, filename, PREC);
  } catch (const std::exception &e) {
    cerr << e.what() << endl;
    exit(1);
  }
}

// Get next line from stream F, skipping empty lines and comments.
string getnextline(ifstream &F) {
  return NRG::Tools::next_data_line(F);
}
