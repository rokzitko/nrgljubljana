#include "nrg-general.hpp"
#include "sym-SPSU2T-impl.hpp"
#include "sym-SPSU2T.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2T(const Params &P)
{
  return std::make_unique<SymmetrySPSU2T<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2T(const Params &P)
{
  return std::make_unique<SymmetrySPSU2T<cmpl>>(P);
}

}