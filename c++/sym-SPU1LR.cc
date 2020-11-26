#include "nrg-general.hpp"
#include "sym-SPU1LR-impl.hpp"
#include "sym-SPU1LR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_SPU1LR(const Params &P)
{
  return std::make_unique<SymmetrySPU1LR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPU1LR(const Params &P)
{
  return std::make_unique<SymmetrySPU1LR<cmpl>>(P);
}

}