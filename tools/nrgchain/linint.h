// Discretization ODE solver for NRG
//
// ** Interpolation code

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "../common/linint.hpp"

// Structures for storing tabulated data, such as rho(omega).
typedef std::pair<double, double> Pair;
using Vec = std::vector<Pair>;

struct ExitError {
  [[noreturn]] void operator()(const std::string &message) const {
    throw std::runtime_error(message);
  }
};

using LinInt = NRG::Tools::LinIntBase<Vec, ExitError>;
using IntLinInt = NRG::Tools::IntLinIntBase<Vec, ExitError>;
