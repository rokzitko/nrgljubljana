#include "nrg-general.hpp"
#include "sym-SPSU2T-impl.hpp"
#include "sym-SPSU2T.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2T(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2T<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2T(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2T<cmpl>>(P, allfields);
}

}