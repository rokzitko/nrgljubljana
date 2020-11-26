#include "nrg-general.hpp"
#include "sym-ISOLR-impl.hpp"
#include "sym-ISOLR.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_ISOLR(const Params &P)
{
  return std::make_unique<SymmetryISOLR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISOLR(const Params &P)
{
  return std::make_unique<SymmetryISOLR<cmpl>>(P);
}

template <>
std::unique_ptr<Symmetry<double>> mk_ISO2LR(const Params &P)
{
  return std::make_unique<SymmetryISO2LR<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISO2LR(const Params &P)
{
  return std::make_unique<SymmetryISO2LR<cmpl>>(P);
}

}