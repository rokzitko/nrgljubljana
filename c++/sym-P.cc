#include "nrg-general.hpp"
#include "sym-P-impl.hpp"
#include "sym-P.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_P(const Params &P)
{
  return std::make_unique<SymmetryP<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_P(const Params &P)
{
  return std::make_unique<SymmetryP<cmpl>>(P);
}

}