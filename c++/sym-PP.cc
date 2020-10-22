#include "nrg-general.h"
#include "sym-PP-impl.h"
#include "sym-PP.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_PP(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryPP<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_PP(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryPP<cmpl>>(P, allfields);
}
