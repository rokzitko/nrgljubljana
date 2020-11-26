#include "nrg-general.hpp"
#include "sym-QSZ-impl.hpp"
#include "sym-QSZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSZ(const Params &P)
{
  return std::make_unique<SymmetryQSZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZ(const Params &P)
{
  return std::make_unique<SymmetryQSZ<cmpl>>(P);
}

}
