#include "nrg-general.hpp"
#include "sym-DBLQSZ-impl.hpp"
#include "sym-DBLQSZ.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_DBLQSZ(const Params &P)
{
  return std::make_unique<SymmetryDBLQSZ<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_DBLQSZ(const Params &P)
{
  return std::make_unique<SymmetryDBLQSZ<cmpl>>(P);
}

}
