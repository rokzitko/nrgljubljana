#include "nrg-general.h"
#include "sym-QSC3-impl.h"
#include "sym-QSC3.h" // include for consistency

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSC3(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSC3<cmpl>>(P, allfields);
}
