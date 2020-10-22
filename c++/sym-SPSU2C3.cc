#include "nrg-general.h"
#include "sym-SPSU2C3-impl.h"
#include "sym-SPSU2C3.h" // include for consistency

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2C3(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2C3<cmpl>>(P, allfields);
}
