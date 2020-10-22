#include "nrg-general.h"
#include "sym-QS-impl.h"
#include "sym-QS.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QS(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQS<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QS(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQS<cmpl>>(P, allfields);
}
