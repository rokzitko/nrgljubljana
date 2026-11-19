// Channel-mixing discretization for NRG
// ** The configuration of the stages, from the [param] entries

#ifndef _mixchain_configuration_hpp_
#define _mixchain_configuration_hpp_

#include <cmath>
#include <map>
#include <stdexcept>
#include <string>

#include "../common/lambda.hpp"
#include "../common/representative_energy.hpp"
#include "../common/tabulated_density.hpp"
#include "branches.hpp"
#include "chain.hpp"
#include "load.hpp"
#include "mesh.hpp"
#include "parser.hpp"
#include "star.hpp"

namespace NRG::MixChain {

struct Configuration {
  GammaOptions gamma;
  StarOptions star;
  bool mmax_from_nmax{};
  ChainOptions chain;
  bool discretization_files{}; // write the chain also as one file per matrix element, beside chain.dat
};

// The keys of the star stage. cquad holds the integration settings, which are not parameters of the file.
inline void read_star_configuration(const Params &P, const NRG::Tools::CquadOptions &cquad, Configuration &configuration) {
  configuration.gamma.prefix                = P.Pstr("dos_prefix", "Gamma");
  configuration.gamma.channels              = P.Pint("channels", 1);
  configuration.gamma.bandrescale           = P.P("bandrescale", 1.0);
  configuration.gamma.hermiticity_tolerance = P.P("hermiticity_tolerance", 1e-8);

  auto &star       = configuration.star;
  star.Lambda      = NRG::Tools::LambdaCache(P.P("Lambda", 2.0));
  star.z           = P.P("z", 1.0);
  star.bandrescale = configuration.gamma.bandrescale;
  star.adapt       = P.Pbool("adapt", false);
  star.hardgap     = P.Pbool("hardgap", false);
  star.boundary    = P.P("boundary", 0.0);
  star.mesh_weight = mesh_weight_from_string(P.Pstr("mesh_weight", "frobenius"));
  star.interpolation =
    NRG::Tools::parse_density_interpolation_method(P.Pstr("density_interpolation", "linear"));
  star.branches.ordering = branch_ordering_from_string(P.Pstr("branch_ordering", "tracked"));
  star.split_blocks      = P.Pbool("split_blocks", true);
  star.allowed_error     = P.P("allowed_error", 1e-10);
  star.cquad             = cquad;

  // The star stage needs mMAX, the number of intervals. Nmax is the length of the Wilson chain, which only the chain
  // stage uses; it serves here as the customary default for mMAX.
  if (P.contains("mMAX")) {
    const auto mmax = P.Pint("mMAX", 0);
    if (mmax <= 0) throw std::invalid_argument("mMAX must be greater than 0.");
    star.mMAX = static_cast<unsigned int>(mmax);
  } else if (P.contains("Nmax")) {
    const auto nmax = P.Pint("Nmax", 0);
    if (nmax <= 0) throw std::invalid_argument("Nmax must be greater than 0.");
    star.mMAX                  = static_cast<unsigned int>(2 * nmax);
    configuration.mmax_from_nmax = true;
  } else {
    throw std::invalid_argument("Either mMAX or Nmax must be given.");
  }
}

// The keys of the chain stage.
inline void read_chain_configuration(const Params &P, Configuration &configuration) {
  if (!P.contains("Nmax")) throw std::invalid_argument("Nmax must be given to build the chain.");
  const auto nmax = P.Pint("Nmax", 0);
  if (nmax <= 0) throw std::invalid_argument("Nmax must be greater than 0.");
  configuration.chain.Nmax = static_cast<unsigned int>(nmax);

  configuration.chain.rank_tolerance = P.P("rank_tolerance", 1e-20);
  if (!(std::isfinite(configuration.chain.rank_tolerance) && configuration.chain.rank_tolerance > 0.0))
    throw std::invalid_argument("rank_tolerance must be a positive finite number.");

  configuration.chain.sensitivity_tolerance = P.P("sensitivity_tolerance", 1e-10);
  if (!(std::isfinite(configuration.chain.sensitivity_tolerance) && configuration.chain.sensitivity_tolerance > 0.0))
    throw std::invalid_argument("sensitivity_tolerance must be a positive finite number.");

  configuration.discretization_files = P.Pbool("discretization_files", false);

  configuration.chain.gauge = chain_gauge_from_string(P.Pstr("chain_gauge", "polar"));
  configuration.chain.nambu_tolerance = P.P("nambu_tolerance", 1e-8);
  if (!(std::isfinite(configuration.chain.nambu_tolerance) && configuration.chain.nambu_tolerance > 0.0))
    throw std::invalid_argument("nambu_tolerance must be a positive finite number.");
}

// The configuration of the stages asked for, from the [param] entries.
inline Configuration read_configuration(const std::map<std::string, std::string> &entries, const bool star, const bool chain,
                                        const NRG::Tools::CquadOptions &cquad = {}) {
  const Params P(entries);
  Configuration configuration;
  if (star) read_star_configuration(P, cquad, configuration);
  if (chain) read_chain_configuration(P, configuration);
  return configuration;
}

} // namespace NRG::MixChain

#endif
