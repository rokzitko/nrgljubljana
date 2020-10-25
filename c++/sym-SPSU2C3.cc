#include "nrg-general.hpp"
#include "sym-SPSU2C3-impl.hpp"
#include "sym-SPSU2C3.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2C3(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2C3<cmpl>>(P, allfields);
}

}