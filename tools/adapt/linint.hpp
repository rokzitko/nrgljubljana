// Discretization ODE solver for NRG
// ** Interpolation code

#ifndef _adapt_linint_hpp_
#define _adapt_linint_hpp_

#include <stdexcept>
#include <utility>
#include <vector>

#include "../common/linint.hpp"

namespace NRG::Adapt {

// Structures for storing tabulated data, such as rho(omega).
using Pair = std::pair<double, double>;
using Vec = std::vector<Pair>;

struct ThrowError {
  [[noreturn]] void operator()(const std::string &message) const {
    throw std::runtime_error(message);
  }
};

using LinInt = NRG::Tools::LinIntBase<Vec, ThrowError>;
using IntLinInt = NRG::Tools::IntLinIntBase<Vec, ThrowError>;

} // namespace

#endif
