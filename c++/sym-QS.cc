#include "nrg-general.hpp"
#include "sym-QS-impl.hpp"
#include "sym-QS.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QS(const Params &P)
{
  return std::make_unique<SymmetryQS<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QS(const Params &P)
{
  return std::make_unique<SymmetryQS<cmpl>>(P);
}

}
