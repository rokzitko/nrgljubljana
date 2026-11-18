// Channel-mixing discretization for NRG
// ** The stages in memory: from Gamma to the star and the Wilson chain of each z

#ifndef _mixchain_discretize_hpp_
#define _mixchain_discretize_hpp_

#include <ostream>
#include <vector>

#include "chain.hpp"
#include "chain_io.hpp"
#include "configuration.hpp"
#include "load.hpp"
#include "report.hpp"
#include "star.hpp"

namespace NRG::MixChain {

// What the chain file records of the star the chain was built from.
template<typename S> ChainFileHeader chain_header(const Star<S> &star) {
  return {star.z, star.Lambda, star.bandrescale, star.untabulated_known ? star.untabulated_from : 0.0,
          star.untabulated_known ? star.untabulated_to : 0.0};
}

// The chain of a star, with the sensitivity of the chain to the star in its diagnostics.
template<typename S> Chain<S> make_chain(const Star<S> &star, const ChainOptions &options) {
  auto chain             = build_chain(star, options);
  const auto sensitivity = star_sensitivity(star, options);
  chain.diagnostics.max_star_sensitivity      = sensitivity.largest;
  chain.diagnostics.max_star_sensitivity_site = sensitivity.largest_site;
  chain.diagnostics.sensitive_from_site       = sensitivity.from_site;
  return chain;
}

template<typename S> struct Discretization {
  Star<S> star;
  Chain<S> chain;
};

// The star and the chain for each z, without files: the z-independent setup once, then per z the star and its chain.
// The reports go to out, as the executable writes them.
template<typename S>
std::vector<Discretization<S>> discretize(const GammaInput<S> &input, const Configuration &configuration, const std::vector<double> &zs,
                                          std::ostream &out) {
  StarDiscretizer<S> discretizer(input, configuration.star);
  report_blocks(discretizer, configuration.star.split_blocks, out);
  std::vector<Discretization<S>> result;
  result.reserve(zs.size());
  for (const auto z : zs) {
    auto star = discretizer.star(z);
    report_star(star, out);
    auto chain = make_chain(star, configuration.chain);
    report_chain(chain, chain_header(star), configuration.chain.sensitivity_tolerance, out);
    result.push_back({std::move(star), std::move(chain)});
  }
  return result;
}

} // namespace NRG::MixChain

#endif
