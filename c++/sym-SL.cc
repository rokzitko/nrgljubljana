#include "nrg-general.h"
#include "sym-SL-impl.h"
#include "sym-SL.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SL(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SL(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL<cmpl>>(P, allfields);
}
