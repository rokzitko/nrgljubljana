#include "nrg-general.hpp"
#include "sym-SL3-impl.hpp"
#include "sym-SL3.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SL3(const Params &P)
{
  return std::make_unique<SymmetrySL3<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SL3(const Params &P)
{
  return std::make_unique<SymmetrySL3<cmpl>>(P);
}

}