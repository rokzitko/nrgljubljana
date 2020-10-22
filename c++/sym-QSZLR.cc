#include "nrg-general.h"
#include "sym-QSZLR-impl.h"
#include "sym-QSZLR.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QSZLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZLR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZLR<cmpl>>(P, allfields);
}
