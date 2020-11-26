#include "nrg-general.hpp"
#include "sym-ISOSZLR-impl.hpp"
#include "sym-ISOSZLR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_ISOSZLR(const Params &P)
{
  return std::make_unique<SymmetryISOSZLR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISOSZLR(const Params &P)
{
  return std::make_unique<SymmetryISOSZLR<cmpl>>(P);
}

}