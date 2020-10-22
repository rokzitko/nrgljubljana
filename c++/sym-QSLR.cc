#include "nrg-general.h"
#include "sym-QSLR-impl.h"
#include "sym-QSLR.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QSLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSLR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSLR<cmpl>>(P, allfields);
}
