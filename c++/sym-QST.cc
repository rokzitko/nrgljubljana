#include "nrg-general.h"
#include "sym-QST-impl.h"
#include "sym-QST.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_QST(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQST<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QST(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQST<cmpl>>(P, allfields);
}
