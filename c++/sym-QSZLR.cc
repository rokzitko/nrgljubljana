#include "nrg-general.hpp"
#include "sym-QSZLR-impl.hpp"
#include "sym-QSZLR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSZLR(const Params &P)
{
  return std::make_unique<SymmetryQSZLR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZLR(const Params &P)
{
  return std::make_unique<SymmetryQSZLR<cmpl>>(P);
}

}