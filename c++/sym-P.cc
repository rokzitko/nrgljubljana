#include "nrg-general.hpp"
#include "sym-P-impl.hpp"
#include "sym-P.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_P(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryP<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_P(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryP<cmpl>>(P, allfields);
}

}