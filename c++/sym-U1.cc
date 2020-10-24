#include "nrg-general.hpp"
#include "sym-U1-impl.hpp"
#include "sym-U1.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_U1(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryU1<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_U1(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryU1<cmpl>>(P, allfields);
}
