#include "nrg-general.hpp"
#include "sym-QSC3-impl.hpp"
#include "sym-QSC3.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSC3(const Params &P)
{
  return std::make_unique<SymmetryQSC3<cmpl>>(P);
}

}