#include "nrg-general.hpp"
#include "sym-QJ-impl.hpp"
#include "sym-QJ.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QJ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQJ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QJ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQJ<cmpl>>(P, allfields);
}
