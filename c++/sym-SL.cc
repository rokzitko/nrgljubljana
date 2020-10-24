#include "nrg-general.hpp"
#include "sym-SL-impl.hpp"
#include "sym-SL.hpp" // include for consistency

template <>
std::unique_ptr<Symmetry<double>> mk_SL(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL<double>>(P, allfields);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_SL(const Params &P, Allfields &allfields)
{
  return std::make_unique<SymmetrySL<cmpl>>(P, allfields);
}
