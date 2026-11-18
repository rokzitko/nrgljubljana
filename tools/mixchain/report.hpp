// Channel-mixing discretization for NRG
// ** The reports of the star and the chain stage

#ifndef _mixchain_report_hpp_
#define _mixchain_report_hpp_

#include <complex>
#include <cstddef>
#include <iostream>
#include <limits>
#include <ostream>
#include <string>
#include <utility>

#include "blocks.hpp"
#include "chain.hpp"
#include "chain_io.hpp"
#include "star.hpp"
#include "types.hpp"

namespace NRG::MixChain {

// The blocks of Gamma, which do not depend on z, so they are reported once.
template<typename S> void report_blocks(const StarDiscretizer<S> &discretizer, const bool split_blocks, std::ostream &out) {
  if (!split_blocks) return;
  if (discretizer.blocks().size() > 1)
    out << "# blocks: " << blocks_name(discretizer.blocks()) << std::endl;
  else
    out << "# blocks: none (Gamma does not split)" << std::endl;
}

template<typename S> void report_star(const Star<S> &star, std::ostream &out) {
  out << "# levels=" << star.levels.size() << " complex=" << (is_complex_v<S> ? 1 : 0) << std::endl;
  for (std::size_t b = 0; b < star.diagnostics.size(); b++) {
    const auto &diagnostics = star.diagnostics[b];
    // With a single block the lines carry no prefix, as they did before blocks existed.
    const auto block = star.blocks.size() > 1 ? "block " + blocks_name({star.blocks[b]}) + ": " : std::string();
    out << "# " << block << "max_interval_deviation=" << diagnostics.max_interval_deviation
        << " at omega=" << diagnostics.max_interval_omega << std::endl;
    out << "# " << block << "max_cquad_error=" << diagnostics.max_cquad_error << std::endl;
    out << "# " << block << "crossings: " << diagnostics.crossings_pos.size() << " positive, "
        << diagnostics.crossings_neg.size() << " negative" << std::endl;

    // Where the mesh reaches below the innermost tabulated frequency, the density is the constant continuation of
    // the input: exact for a flat band, an approximation for anything with structure at low frequency.
    // These frequencies are compared with the input grid and can sit very close to an accumulation point, so they
    // are printed with every digit that distinguishes them.
    const auto precision = out.precision(std::numeric_limits<double>::max_digits10);
    for (const auto &[name, coverage] : {std::pair{"POS", &diagnostics.coverage_pos},
                                         std::pair{"NEG", &diagnostics.coverage_neg}}) {
      const auto prefix = "# " + block + name + ": ";
      if (coverage->fixed_mesh_fallback)
        out << prefix << "Gamma vanishes on this branch; fixed mesh used" << std::endl;
      if (coverage->collapsed_levels > 0)
        out << prefix << coverage->collapsed_levels
            << " representative energy levels are indistinguishable from the accumulation point "
            << coverage->accumulation_point << " in double precision" << std::endl;
      if (coverage->unresolved_intervals == 0) continue;
      out << prefix << coverage->unresolved_intervals << " of " << star.mMAX + 1
          << " intervals contain no tabulated point of the input, carrying " << coverage->unresolved_share()
          << " of the weight of this branch, the outermost being [" << coverage->unresolved_to << ", "
          << coverage->unresolved_from << "]; ";
      if (coverage->continued())
        out << "the input ends at omega=" << coverage->innermost_input
            << " and below that the density is its constant continuation" << std::endl;
      else
        out << "there the star follows the interpolant between neighbouring points" << std::endl;
    }
    out.precision(precision);
  }

  // The trace of Theta agrees with that of the integral of Gamma by construction; the off-diagonal elements only
  // approximately, because a single level per interval and branch cannot follow a rotating eigenvector.
  const auto trace_star  = std::real(star.theta.trace());
  const auto trace_exact = std::real(star.theta_exact.trace());
  out << "# tr(theta)=" << trace_star << " tr(int Gamma)=" << trace_exact << " difference="
      << trace_star - trace_exact << std::endl;
  const auto scale = star.theta_exact.norm();
  if (scale > 0.0)
    out << "# ||theta - int Gamma||/||int Gamma||=" << (star.theta - star.theta_exact).norm() / scale << std::endl;
}

template<typename S>
void report_chain(const Chain<S> &chain, const ChainFileHeader &header, const double sensitivity_tolerance,
                  std::ostream &out) {
  const auto &d = chain.diagnostics;
  out << "# chain: sites=" << chain.Nmax + 1 << " channels=" << chain.channels
      << " gauge=" << chain_gauge_name(chain.gauge) << std::endl;
  if (chain.gauge == ChainGauge::nambu)
    out << "# the Nambu structure of the blocks holds to " << d.max_nambu_deviation
        << " of their largest element" << std::endl;
  if (chain.blocks.size() > 1) out << "# blocks: " << blocks_name(chain.blocks) << std::endl;
  out << "# levels=" << d.levels << " coupled_levels=" << d.coupled_levels << " theta_rank=" << d.theta_rank
      << " theta_condition=" << d.theta_condition
      << " min_residual_condition=" << d.min_residual_condition << std::endl;
  // From the site where the chain reaches the untabulated region of the input, its coefficients rest on the constant
  // continuation of the input rather than on data. Printed in the units of the input.
  if (const auto continued = first_continued_site(chain, header.untabulated_to - header.untabulated_from)) {
    out << "# from site " << *continued << " on, the chain samples ";
    if (header.untabulated_from > 0.0) out << header.untabulated_from * header.bandrescale << " < ";
    out << "|omega| < " << header.untabulated_to * header.bandrescale
        << ", where Gamma is not tabulated (extend the input grid to lower |omega|, or lower Nmax)" << std::endl;
  }
  // How well a star in double precision determines the chain at all.
  out << "# max_star_sensitivity=" << d.max_star_sensitivity << " at site " << d.max_star_sensitivity_site << std::endl;
  if (d.sensitive_from_site)
    out << "# from site " << *d.sensitive_from_site << " on, the chain is determined by the star only to "
        << sensitivity_tolerance << " or worse (" << d.max_star_sensitivity << " at site "
        << d.max_star_sensitivity_site << "):" << std::endl
        << "# levels near an accumulation point of the mesh at a finite energy are not resolved in double precision."
        << std::endl
        << "# The coefficients from there on are not reproducible beyond that; a smaller Nmax avoids it." << std::endl;
  // A Theta of lower rank is a property of Gamma, and the chain is exact for it; a drop further down is not.
  if (d.theta_rank < chain.channels) {
    const auto decoupled = chain.channels - d.theta_rank;
    out << "# Theta has rank " << d.theta_rank << " of " << chain.channels << ": " << decoupled << " combination"
        << (decoupled == 1 ? "" : "s") << " of the impurity orbitals do" << (decoupled == 1 ? "es" : "")
        << " not couple to the bath, and its part of the chain is zero" << std::endl;
  }
  // A drop is reported for each block where it happens: what limits the chain is a property of that block's star.
  for (std::size_t b = 0; b < chain.blocks.size(); b++) {
    const auto &part = chain.block_diagnostics[b];
    if (!part.rank_drop_site) continue;
    std::cerr << "mixchain: warning: the rank of the hopping drops below " << part.theta_rank << " at site "
              << *part.rank_drop_site << " (smallest rank " << part.min_rank << "), in block "
              << blocks_name({chain.blocks[b]}) << ": the Krylov space of the star is exhausted in some direction, "
              << "and the chain is zero in it from there on. The block has " << part.coupled_levels << " of "
              << part.levels << " levels with nonzero coupling; the others lie where Gamma vanishes on the mesh or "
              << "collapsed onto an accumulation point. A mesh that follows this block (adapt=true, split_blocks=true) "
              << "or a smaller Nmax avoids it." << std::endl;
  }
}

} // namespace NRG::MixChain

#endif
