#include "nrg-general.hpp"
#include "sym-SU2-impl.hpp"
#include "sym-SU2.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SU2(const Params &P)
{
  return std::make_unique<SymmetrySU2<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SU2(const Params &P)
{
  return std::make_unique<SymmetrySU2<cmpl>>(P);
}

}