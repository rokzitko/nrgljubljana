#include "nrg-general.hpp"
#include "sym-QSZTZ-impl.hpp"
#include "sym-QSZTZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSZTZ(const Params &P)
{
  return std::make_unique<SymmetryQSZTZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZTZ(const Params &P)
{
  return std::make_unique<SymmetryQSZTZ<cmpl>>(P);
}

}