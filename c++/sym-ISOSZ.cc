#include "nrg-general.hpp"
#include "sym-ISOSZ-impl.hpp"
#include "sym-ISOSZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_ISOSZ(const Params &P)
{
  return std::make_unique<SymmetryISOSZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISOSZ(const Params &P)
{
  return std::make_unique<SymmetryISOSZ<cmpl>>(P);
}

}