// Discretization ODE solver for NRG
//
// ** Interpolation code

#include <cstdlib>
#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "../common/linint.hpp"

// Structures for storing tabulated data, such as rho(omega).
typedef std::pair<double, double> Pair;
using Vec = std::vector<Pair>;

struct ExitError {
  [[noreturn]] void operator()(const std::string &message) const {
    std::cerr << "ERROR: " << message << std::endl;
    std::exit(1);
  }
};

using LinInt = NRG::Tools::LinIntBase<Vec, ExitError>;
using IntLinInt = NRG::Tools::IntLinIntBase<Vec, ExitError>;
