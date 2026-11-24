#ifndef NRG_INSTANTIATE_WORKFLOW_HPP
#define NRG_INSTANTIATE_WORKFLOW_HPP

#include <array>
#include <string>
#include <vector>

#include "../nrgchain/nrgchain.hpp"

namespace NRG::Workflow::detail {
// Called under the workflow preparation guard, with the calculation CWD active.
double instantiate(const std::string &parameter_file, const std::string &template_directory, unsigned nmax);
// The Wilson coefficients nrgchain computes for the parameter file (Nmax replaced by nmax), with the z, Lambda and
// bandrescale they were computed for. Intermediate files go to a staging directory that is removed.
struct ComputedWilson {
  NRG::Tools::NrgChain::WilsonData data;
  double z, Lambda, bandrescale;
};
ComputedWilson nrgchain_coefficients(const std::string &parameter_file, unsigned nmax);
struct QSMixing {
  std::array<double, 4> coupling; // polar V in the impurity-by-bath basis, before 1/sqrt(pi)
  std::vector<double> hopping, onsite; // the real off-diagonal elements of T and E
};

// The same with supplied Wilson coefficients instead of nrgchain, for the z, Lambda and bandrescale given.
double instantiate(const std::string &parameter_file, const std::string &template_directory, unsigned nmax,
                   const NRG::Tools::NrgChain::WilsonData &wilson, double z, double Lambda, double bandrescale,
                   const QSMixing *mixing = nullptr);
}

#endif
