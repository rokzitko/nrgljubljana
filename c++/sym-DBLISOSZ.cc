#include "nrg-general.hpp"
#include "sym-DBLISOSZ-impl.hpp"
#include "sym-DBLISOSZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_DBLISOSZ(const Params &P)
{
  return std::make_unique<SymmetryDBLISOSZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_DBLISOSZ(const Params &P)
{
  return std::make_unique<SymmetryDBLISOSZ<cmpl>>(P);
}

}