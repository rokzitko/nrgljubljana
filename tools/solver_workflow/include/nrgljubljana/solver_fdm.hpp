#ifndef NRGLJUBLJANA_SOLVER_FDM_HPP
#define NRGLJUBLJANA_SOLVER_FDM_HPP

#include <complex>
#include <map>
#include <span>
#include <string>
#include <vector>

#include <nrgljubljana/numerics.hpp>

namespace NRG::Workflow {

struct FdmOperators {
  std::vector<std::string> electron, interaction; // d_i and [d_i,H_int] in orbital order
  std::vector<std::string> hartree; // row-major operator names for <{d_i,q_j^dagger}>
  std::vector<std::string> spin; // optional triplet operators, one per orbital
};

using FdmMatrix = std::vector<std::complex<double>>;

struct FdmResult {
  std::vector<double> omega;
  unsigned orbitals = 0;
  std::map<std::string, double> expectations;
  FdmMatrix hartree; // row-major
  FdmMatrix green, left, right, composite, self_energy, spin; // [frequency, row, column]
};

// Real, spin-symmetric FDM output (customfdm and specd/spect binary pairs).
// Read committed twist directories in the given order and average with equal weight.
// The fermionic Green function uses one spin. Correlators retain the normalization
// of the named template's operators. A missing file throws; a zero operator may
// have a valid empty bin file. Mesh requirements are those of kramers_kronig.
// A final_derfd_gamma pass convolves the piecewise-linear density with the
// derivative-of-Fermi kernel, using constant continuation of the sampled
// endpoints for this convolution.
FdmResult collect_fdm(const std::vector<std::string> &twist_directories, const FdmOperators &operators,
                      const std::vector<double> &omega, const NRG::Numerics::BroadeningOptions &broadening);

// One frequency, row-major. Throws on incompatible sizes, nonfinite values, or
// singular G. The estimator returns the raw Sigma, without causal projection.
FdmMatrix symmetric_self_energy(unsigned orbitals, std::span<const std::complex<double>> hartree,
                               std::span<const std::complex<double>> green, std::span<const std::complex<double>> left,
                               std::span<const std::complex<double>> right, std::span<const std::complex<double>> composite);

} // namespace NRG::Workflow

#endif
