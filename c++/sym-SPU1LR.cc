#include "nrg-general.h"
#include "sym-SPU1LR-impl.h"
#include "sym-SPU1LR.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SPU1LR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPU1LR<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SPU1LR(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySPU1LR<cmpl>>(P, allfields);
}
