#include "nrg-general.hpp"
#include "sym-QST-impl.hpp"
#include "sym-QST.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QST(const Params &P)
{
  return std::make_unique<SymmetryQST<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QST(const Params &P)
{
  return std::make_unique<SymmetryQST<cmpl>>(P);
}

}