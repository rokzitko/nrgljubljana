#include "nrg-general.hpp"
#include "sym-PP-impl.hpp"
#include "sym-PP.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_PP(const Params &P)
{
  return std::make_unique<SymmetryPP<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_PP(const Params &P)
{
  return std::make_unique<SymmetryPP<cmpl>>(P);
}

}