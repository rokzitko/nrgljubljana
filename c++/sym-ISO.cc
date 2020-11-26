#include "nrg-general.hpp"
#include "sym-ISO-impl.hpp"
#include "sym-ISO.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_ISO(const Params &P)
{
  return std::make_unique<SymmetryISO<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISO(const Params &P)
{
  return std::make_unique<SymmetryISO<cmpl>>(P);
}

template <>
std::unique_ptr<Symmetry<double>> mk_ISO2(const Params &P)
{
  return std::make_unique<SymmetryISO2<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISO2(const Params &P)
{
  return std::make_unique<SymmetryISO2<cmpl>>(P);
}

}