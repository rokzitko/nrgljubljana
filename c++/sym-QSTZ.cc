#include "nrg-general.hpp"
#include "sym-QSTZ-impl.hpp"
#include "sym-QSTZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSTZ(const Params &P)
{
  return std::make_unique<SymmetryQSTZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSTZ(const Params &P)
{
  return std::make_unique<SymmetryQSTZ<cmpl>>(P);
}

}