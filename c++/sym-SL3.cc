#include "nrg-general.hpp"
#include "sym-SL3-impl.hpp"
#include "sym-SL3.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SL3(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL3<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SL3(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL3<cmpl>>(P, allfields);
}
