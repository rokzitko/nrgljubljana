#include <memory>

#include "nrg-general.hpp"
#include "sym-SPZ4-impl.hpp"
#include "sym-SPZ4.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPZ4(const Params &P)
{
  return std::make_unique<SymmetrySPZ4<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPZ4(const Params &P)
{
  return std::make_unique<SymmetrySPZ4<cmpl>>(P);
}

}
