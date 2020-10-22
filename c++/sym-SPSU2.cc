#include "nrg-general.h"
#include "sym-SPSU2-impl.h"
#include "sym-SPSU2.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2<cmpl>>(P, allfields);
}
