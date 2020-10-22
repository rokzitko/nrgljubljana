#include "nrg-general.h"
#include "sym-DBLISOSZ-impl.h"
#include "sym-DBLISOSZ.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_DBLISOSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryDBLISOSZ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_DBLISOSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryDBLISOSZ<cmpl>>(P, allfields);
}
