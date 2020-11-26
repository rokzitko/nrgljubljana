#include "nrg-general.hpp"
#include "sym-SL-impl.hpp"
#include "sym-SL.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SL(const Params &P)
{
  return std::make_unique<SymmetrySL<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SL(const Params &P)
{
  return std::make_unique<SymmetrySL<cmpl>>(P);
}

}