#include "nrg-general.h"
#include "sym-SU2-impl.h"
#include "sym-SU2.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySU2<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySU2<cmpl>>(P, allfields);
}
