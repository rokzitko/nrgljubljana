#include "nrg-general.hpp"
#include "sym-SPSU2-impl.hpp"
#include "sym-SPSU2.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2(const Params &P)
{
  return std::make_unique<SymmetrySPSU2<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2(const Params &P)
{
  return std::make_unique<SymmetrySPSU2<cmpl>>(P);
}

}