#include "nrg-general.h"
#include "sym-QSTZ-impl.h"
#include "sym-QSTZ.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QSTZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSTZ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSTZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSTZ<cmpl>>(P, allfields);
}
