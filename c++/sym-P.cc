#include "nrg-general.h"
#include "sym-P-impl.h"
#include "sym-P.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_P(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryP<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_P(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryP<cmpl>>(P, allfields);
}
