#include "nrg-general.hpp"
#include "sym-QJ-impl.hpp"
#include "sym-QJ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QJ(const Params &P)
{
  return std::make_unique<SymmetryQJ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QJ(const Params &P)
{
  return std::make_unique<SymmetryQJ<cmpl>>(P);
}

}