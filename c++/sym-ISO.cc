#include "nrg-general.h"
#include "sym-ISO-impl.h"
#include "sym-ISO.h" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_ISO(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISO<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISO(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISO<cmpl>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<double>> mk_ISO2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISO2<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_ISO2(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetryISO2<cmpl>>(P, allfields);
}
