#include "nrg-general.hpp"
#include "sym-QSLR-impl.hpp"
#include "sym-QSLR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSLR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSLR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSLR<cmpl>>(P, allfields);
}

}