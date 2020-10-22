#include "nrg-general.h"
#include "sym-SPSU2LR-impl.h"
#include "sym-SPSU2LR.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SPSU2LR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2LR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPSU2LR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPSU2LR<cmpl>>(P, allfields);
}
