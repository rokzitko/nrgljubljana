#include "nrg-general.hpp"
#include "sym-QSZ-impl.hpp"
#include "sym-QSZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZ<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZ(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryQSZ<cmpl>>(P, allfields);
}

}