#ifndef NRGLJUBLJANA_SOLVER_CHAIN_HPP
#define NRGLJUBLJANA_SOLVER_CHAIN_HPP

#include <map>
#include <optional>
#include <string>
#include <vector>

#include <nrgljubljana/solver_workflow.hpp>

namespace NRG::Workflow {

// Gamma(omega), the N x N Hermitian positive semidefinite hybridization matrix, tabulated: for N=1 in the convention
// of the dos file of adapt and nrgchain, -Im Delta(omega+i0). See tools/mixchain/README.md.
struct HybridizationFunction {
  unsigned channels = 0;          // N
  std::vector<double> omega;      // strictly increasing, both signs
  std::vector<double> real;       // Re Gamma_ij(omega_k) at (k*N + i)*N + j, 0-based
  std::vector<double> imaginary;  // empty for a real Gamma, otherwise like real
};

// A matrix, row-major; imaginary is empty for a real chain.
struct ChainMatrix {
  unsigned rows = 0, cols = 0;
  std::vector<double> real, imaginary;
};

// The Wilson chain of one twist, as mixchain writes it to chain.dat: V as built, E and T in the units of the input
// (multiplied by bandrescale). The sites are 0..nmax.
struct WilsonChain {
  double z = 0, Lambda = 0, bandrescale = 1;
  unsigned channels = 0, nmax = 0;
  std::string gauge;                         // polar or nambu
  ChainMatrix V;                             // the impurity coupling
  std::vector<ChainMatrix> E;                // the on-site blocks E_0..E_nmax
  std::vector<ChainMatrix> T;                // the hoppings T_0..T_nmax, the last one out of the chain
  std::vector<std::vector<unsigned>> blocks; // channels discretized together, 0-based

  // Diagnostics, as the header of chain.dat records them (see tools/mixchain/README.md)
  int levels = 0, coupled_levels = 0;
  int theta_rank = 0, min_rank = 0;
  std::optional<unsigned> rank_drop_site;      // the first site where the rank of the hopping drops
  std::optional<unsigned> continued_from_site; // the first site that samples where Gamma is not tabulated
  double theta_condition = 0, min_residual_condition = 0;
  double max_star_sensitivity = 0;
  std::optional<unsigned> sensitive_from_site; // the first site not determined by the star to sensitivity_tolerance
  double max_interval_deviation = 0;           // the largest over the blocks of the star
};

// The Wilson chain of Gamma for each z in twists, without files: mixchain's star and chain stages, the z-independent
// part shared. param_section holds the [param] entries, read with mixchain's names and defaults, except that
// channels and dos_prefix are ignored (Gamma gives N), z is ignored (twists give it), Nmax is replaced by nmax and an
// explicit mMAX is kept (otherwise 2*nmax). discretization, if given, must be Z, the only scheme of mixchain.
// Complex arithmetic exactly when gamma.imaginary is given. Gamma must be positive semidefinite: an eigenvalue
// negative within psd_tolerance of the largest at a node is rounding and set to zero, anything more negative is an
// error. Under the workflow guard; a nonempty diagnostic_file receives mixchain's log, resolved against the
// caller's working directory. Throws on invalid input.
std::vector<WilsonChain> discretize_hybridization(const HybridizationFunction &gamma,
                                                  const std::map<std::string, std::string> &param_section,
                                                  const std::vector<double> &twists, unsigned nmax,
                                                  const std::string &diagnostic_file = {});

// The coefficients of one channel, as nrgchain computes them and the initializer reads them (gammaPolCh = sqrt(theta/pi),
// coefxi, coefzeta): theta, and xi and zeta for the sites 0..nmax, in the units of the input.
struct CoefficientChannel {
  double theta = 0;
  std::vector<double> xi, zeta;
};

struct ChainCoefficients {
  double z = 0, Lambda = 0, bandrescale = 1;
  std::vector<CoefficientChannel> channels; // coefficient channel k is channel k of Gamma
};

// The coefficients of a chain whose channels do not couple: theta_k = V_kk^2, xi_k = T_n(k,k), zeta_k = E_n(k,k).
// Every element left out must be exactly zero: the off-diagonal elements of V, E_n and T_n, and all imaginary parts.
// Otherwise, or for a gauge other than polar, it throws std::invalid_argument naming the element. A diagonal Gamma
// discretized with split_blocks=true (the default) gives such a chain.
ChainCoefficients diagonal_coefficients(const WilsonChain &chain);

// The coefficients nrgchain computes for options.parameter_file in options.calculation_directory, from the FSOL.dat and
// FSOLNEG.dat that discretize writes there (or analytically for band=flat): what instantiate(options) uses, with the same
// parameter checks. nmax replaces Nmax; template_directory is not used. Under the workflow guard; no files are left.
ChainCoefficients nrgchain_coefficients(const InitializerOptions &options, const std::string &diagnostic_file = {});

// instantiate with the given coefficients instead of computing them with nrgchain; everything else as the instantiate
// of solver_workflow.hpp, allowing real one- or two-channel unpolarized QS. Rungs are disabled. The coefficients must cover the sites
// 0..options.nmax (a longer chain is cut), match the coefficient channels of the template, and have the z, Lambda and
// bandrescale of the parameter file (to 1e-12 relative, with the solver's defaults 1, 2 and 1). The chain keys of the
// file (dos, band, adapt, tridiag_method, ...) are not read.
InitializerResult instantiate(const InitializerOptions &options, const ChainCoefficients &coefficients,
                              const std::string &diagnostic_file = {});

// Instantiate a real two-channel QS template with the full polar-gauge chain.
// The template uses bathCoupling(impurity, channel) = V(impurity, channel)/sqrt(pi)
// with 1-based indices, and coefrung(1,0) for the site-zero off-diagonal E.
// Set rungs=true if the off-diagonal elements of E or T are nonzero.
InitializerResult instantiate(const InitializerOptions &options, const WilsonChain &chain,
                              const std::string &diagnostic_file = {});

} // namespace NRG::Workflow

#endif
