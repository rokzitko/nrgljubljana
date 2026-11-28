#include <memory>

#include "nrg-general.hpp"
#include "sym-QSZ4-impl.hpp"
#include "sym-QSZ4.hpp" // include for consistency

namespace NRG {

template <>
std::unique_ptr<Symmetry<double>> mk_QSZ4(const Params &P)
{
  return std::make_unique<SymmetryQSZ4<double>>(P);
}

template <>
std::unique_ptr<Symmetry<cmpl>> mk_QSZ4(const Params &P)
{
  return std::make_unique<SymmetryQSZ4<cmpl>>(P);
}

}
