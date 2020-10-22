#include "nrg-general.h"
#include "sym-QSZTZ-impl.h"
#include "sym-QSZTZ.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QSZTZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZTZ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZTZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZTZ<cmpl>>(P, allfields);
}
