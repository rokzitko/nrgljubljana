#include "nrg-general.hpp"
#include "sym-ISOSZ-impl.hpp"
#include "sym-ISOSZ.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_ISOSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISOSZ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISOSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISOSZ<cmpl>>(P, allfields);
}
