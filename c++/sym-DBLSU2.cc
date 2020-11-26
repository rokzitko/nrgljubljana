#include "nrg-general.hpp"
#include "sym-DBLSU2-impl.hpp"
#include "sym-DBLSU2.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_DBLSU2(const Params &P)
{
  return std::make_unique<SymmetryDBLSU2<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_DBLSU2(const Params &P)
{
  return std::make_unique<SymmetryDBLSU2<cmpl>>(P);
}

}