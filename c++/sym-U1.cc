#include "nrg-general.hpp"
#include "sym-U1-impl.hpp"
#include "sym-U1.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_U1(const Params &P)
{
  return std::make_unique<SymmetryU1<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_U1(const Params &P)
{
  return std::make_unique<SymmetryU1<cmpl>>(P);
}

}