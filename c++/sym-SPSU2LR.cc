#include "nrg-general.hpp"
#include "sym-SPSU2LR-impl.hpp"
#include "sym-SPSU2LR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2LR(const Params &P)
{
  return std::make_unique<SymmetrySPSU2LR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2LR(const Params &P)
{
  return std::make_unique<SymmetrySPSU2LR<cmpl>>(P);
}

}