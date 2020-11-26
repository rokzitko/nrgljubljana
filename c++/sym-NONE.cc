#include "nrg-general.hpp"
#include "sym-NONE-impl.hpp"
#include "sym-NONE.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_NONE(const Params &P)
{
  return std::make_unique<SymmetryNONE<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_NONE(const Params &P)
{
  return std::make_unique<SymmetryNONE<cmpl>>(P);
}

}