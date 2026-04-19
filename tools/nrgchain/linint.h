// Discretization ODE solver for NRG
//
// ** Interpolation code

#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

#include "../common/linint.hpp"

// Structures for storing tabulated data, such as rho(omega).
typedef pair<double, double> Pair;
using Vec = vector<Pair>;

struct ExitError {
  [[noreturn]] void operator()(const std::string &message) const {
    cerr << "ERROR: " << message << endl;
    exit(1);
  }
};

using LinInt = NRG::Tools::LinIntBase<Vec, ExitError>;
using IntLinInt = NRG::Tools::IntLinIntBase<Vec, ExitError>;
