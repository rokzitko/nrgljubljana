#include "nrg-general.hpp"
#include "sym-ISOSZLR-impl.hpp"
#include "sym-ISOSZLR.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_ISOSZLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISOSZLR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISOSZLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISOSZLR<cmpl>>(P, allfields);
}
