#include "nrg-general.hpp"
#include "sym-QS-impl.hpp"
#include "sym-QS.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QS(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQS<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QS(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQS<cmpl>>(P, allfields);
}
