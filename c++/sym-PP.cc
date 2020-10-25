#include "nrg-general.hpp"
#include "sym-PP-impl.hpp"
#include "sym-PP.hpp" // include for consistency

namespace NRG {

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

}