#include "nrg-general.h"
#include "sym-DBLSU2-impl.h"
#include "sym-DBLSU2.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_DBLSU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryDBLSU2<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_DBLSU2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryDBLSU2<cmpl>>(P, allfields);
}
