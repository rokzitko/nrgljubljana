#include "nrg-general.hpp"
#include "sym-SPU1-impl.hpp"
#include "sym-SPU1.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPU1(const Params &P)
{
  return std::make_unique<SymmetrySPU1<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPU1(const Params &P)
{
  return std::make_unique<SymmetrySPU1<cmpl>>(P);
}

}