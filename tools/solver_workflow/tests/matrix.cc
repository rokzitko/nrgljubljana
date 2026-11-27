// Compare the prepared QS chain against the direct spin-degenerate finite-chain resolvent.
#include <nrgljubljana/solver_chain.hpp>
#include <nrgljubljana/solver_fdm.hpp>
#include <nrgljubljana/solver_output.hpp>

#include <Eigen/Dense>

#include <algorithm>
#include <bit>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numbers>
#include <sstream>
#include <stdexcept>

namespace fs = std::filesystem;
using namespace NRG::Workflow;
using Complex = std::complex<double>;

void check(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

void close(Complex computed, Complex expected, double tolerance, const std::string &kind) {
  if (std::abs(computed - expected) > tolerance) {
    std::cerr << std::setprecision(17) << "got " << computed << ", expected " << expected << '\n';
    throw std::runtime_error(kind + ": disagrees with the independent reference.");
  }
}

void run(const fs::path &model, const fs::path &path, bool dense, double bandrescale) {
  fs::create_directories(path);
  fs::create_directories(path / "scratch");
  {
    std::ofstream param(path / "param");
    param << "[extra]\nU1=0\nU2=0\nU12=0\nJ12=0\neps1=-0.16\neps2=0.09\nV12=" << (dense ? 0.05 : 0.0)
          << "\n[param]\nsymtype=QS\nLambda=2\nz=1\nNmax=1\nbandrescale=" << bandrescale
          << "\nT=0.12\nkeep=10000\nkeepenergy=-1"
             "\ndata_has_rescaled_energies=false\nfdm=true\nfdmexpv=true\nsavebins=true\nbroaden=false\nbins=4000\n"
             "ops=n_d1 n_d2 SigmaHartree11 SigmaHartree12 SigmaHartree21 SigmaHartree22\n"
             "spect=sigma_d1-sigma_d1 sigma_d1-sigma_d2 sigma_d2-sigma_d1 sigma_d2-sigma_d2\nspecd=";
    for (const auto *first : {"A_d1", "A_d2", "self_d1", "self_d2"})
      for (const auto *second : {"A_d1", "A_d2", "self_d1", "self_d2"})
        param << first << '-' << second << ' ';
    param << "\nrungs=" << (dense ? "true" : "false") << '\n';
  }
  auto block = [dense](double a, double b, double off) {
    const auto x = dense ? off : 0.0;
    return ChainMatrix{2, 2, {a, x, x, b}, {}};
  };
  WilsonChain chain;
  chain.channels = 2;
  chain.nmax = 1;
  chain.z = 1;
  chain.Lambda = 2;
  chain.bandrescale = bandrescale;
  chain.gauge = "polar";
  chain.V = block(0.35, 0.21, 0.075);
  chain.E = {block(-0.055, 0.11, -0.065), block(0.13, -0.15, 0.035)};
  chain.T = {block(0.24, 0.17, 0.045), block(0.14, 0.09, 0.02)};

  const InitializerOptions options{path.string(), "param", model.string(), 1};
  InitializerResult initial;
  if (dense) {
    initial = instantiate(options, chain);
    auto bad = chain;
    bad.T.front().real[1] += 0.02;
    try {
      instantiate(options, bad);
      throw std::logic_error("Accepted an asymmetric Wilson chain.");
    } catch (const std::invalid_argument &) {}
  } else {
    initial = instantiate(options, diagonal_coefficients(chain));
  }

  Eigen::Matrix<double, 6, 6> h = Eigen::Matrix<double, 6, 6>::Zero();
  const auto b = [](const ChainMatrix &m) {
    return Eigen::Map<const Eigen::Matrix<double, 2, 2, Eigen::RowMajor>>(m.real.data());
  };
  h.topLeftCorner<2, 2>() << -0.16, dense ? 0.05 : 0., dense ? 0.05 : 0., 0.09;
  h.block<2, 2>(0, 2) = b(chain.V) / std::sqrt(std::numbers::pi);
  h.block<2, 2>(2, 0) = h.block<2, 2>(0, 2).transpose();
  h.block<2, 2>(2, 2) = b(chain.E[0]);
  h.block<2, 2>(4, 4) = b(chain.E[1]);
  h.block<2, 2>(2, 4) = b(chain.T[0]);
  h.block<2, 2>(4, 2) = h.block<2, 2>(2, 4).transpose();
  const Eigen::Matrix4d seed = h.topLeftCorner<4, 4>();
  close(initial.ground_energy, 2 * seed.selfadjointView<Eigen::Lower>().eigenvalues().cwiseMin(0.).sum() - b(chain.E[0]).trace(),
        2e-12, "seed ground energy");

  solve(path.string(), (path / "scratch").string(), "solve.log");
  check(fs::is_regular_file(path / "DONE"), "missing DONE");
  const Complex z{0.085, 0.27};
  const Eigen::Matrix<Complex, 6, 6> inverse = z * Eigen::Matrix<Complex, 6, 6>::Identity() - h;
  const Eigen::Matrix<Complex, 6, 6> resolvent = inverse.inverse();
  const Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, 6, 6>> eigen(h);
  check(eigen.info() == Eigen::Success, "Could not diagonalize the independent chain.");
  const auto &e = eigen.eigenvalues();
  const auto &v = eigen.eigenvectors();
  const auto occupation = [](double energy) { return 1.0 / (std::exp(energy / 0.12) + 1.0); };
  for (unsigned i = 0; i < 2; ++i)
    for (unsigned j = 0; j < 2; ++j) {
      const auto op = "A_d" + std::to_string(i + 1) + "-A_d" + std::to_string(j + 1);
      const auto name = spectrum_file_name("specd", op, "FDM") + ".bin";
      Complex g{};
      for (const auto &bin : NRG::Numerics::load_spectrum_bins((path / name).string())) g += bin.weight / (z - bin.frequency);
      close(g, resolvent(i, j), 2e-7, op);

      const auto spin = "sigma_d" + std::to_string(i + 1) + "-sigma_d" + std::to_string(j + 1);
      const auto spin_file = spectrum_file_name("spect", spin, "FDM") + ".bin";
      Complex actual{}, reference{};
      for (const auto &bin : NRG::Numerics::load_spectrum_bins((path / spin_file).string()))
        actual += bin.weight / (z - bin.frequency);
      for (int m = 0; m < 6; ++m)
        for (int n = 0; n < 6; ++n)
          reference += 0.5 * v(i, m) * v(i, n) * v(j, m) * v(j, n) *
                       (occupation(e[m]) - occupation(e[n])) / (z + e[m] - e[n]);
      close(actual, reference, 2e-7, spin);
    }
  const FdmOperators operators{{"A_d1", "A_d2"}, {"self_d1", "self_d2"},
      {"SigmaHartree11", "SigmaHartree12", "SigmaHartree21", "SigmaHartree22"}, {"sigma_d1", "sigma_d2"}};
  std::vector<double> mesh{-1.0, -0.6, -0.2, -0.05, 0.05, 0.2, 0.6, 1.0};
  const NRG::Numerics::BroadeningOptions broadening{.alpha = 0.2, .temperature = 0.12,
      .final_gaussian_gamma = std::nullopt, .final_derfd_gamma = std::nullopt};
  const auto fdm = collect_fdm({path.string()}, operators, mesh, broadening);
  for (const auto value : fdm.self_energy) close(value, 0., 1e-14, "U=0 matrix self-energy");

  auto thermal = broadening;
  thermal.final_derfd_gamma = 0.05;
  const auto smeared = collect_fdm({path.string()}, operators, mesh, thermal);
  const double width = *thermal.final_derfd_gamma * thermal.temperature;
  for (std::size_t i = 0; i < mesh.size(); ++i) {
    // Independent composite midpoint rule in the scaled thermal coordinate,
    // applied to the unbroadened-in-temperature, linearly interpolated data.
    constexpr unsigned subdivisions = 200000;
    constexpr double step_size = 72.0 / subdivisions;
    double reference = 0.0;
    for (unsigned k = 0; k < subdivisions; ++k) {
      const double u = -36.0 + step_size * (k + 0.5);
      const double x = std::clamp(mesh[i] + u * width, mesh.front(), mesh.back());
      const auto next = std::lower_bound(mesh.begin() + 1, mesh.end(), x);
      const auto j = static_cast<std::size_t>(next - mesh.begin());
      const double t = (x - mesh[j - 1]) / (mesh[j] - mesh[j - 1]);
      const double interpolated = std::lerp(fdm.green[4 * (j - 1)].imag(), fdm.green[4 * j].imag(), t);
      reference += step_size * interpolated / (2.0 * (1.0 + std::cosh(u)));
    }
    close(smeared.green[4 * i].imag(), reference, 1e-8, "FDM thermal convolution on a coarse mesh");
  }
  for (const auto value : smeared.self_energy) close(value, 0., 1e-14, "Thermally broadened U=0 self-energy");
}

void kanamori_atom(const fs::path &model, const fs::path &work) {
  const auto path = work / "kanamori-atom";
  fs::create_directories(path);
  {
    std::ofstream param(path / "param");
    param << "[extra]\nU1=0.6\nU2=0.6\nU12=0.4\nJ12=0.1\neps1=-0.54\neps2=-0.46\nV12=0\n"
             "[param]\nsymtype=QS\nLambda=2\nz=1\nNmax=1\ndata_has_rescaled_energies=false\n";
  }
  const ChainCoefficients zero{1, 2, 1, {{0, {0, 0}, {0, 0}}, {0, {0, 0}, {0, 0}}}};
  const auto result = instantiate({path.string(), "param", model.string(), 1}, zero);
  // Energies from independent particle-number sectors. In the two-electron
  // sector the triplet and singlet are at U' ± J. Pair transfer couples the
  // doublons in a block with diagonal levels 2*eps_i+U and off-diagonal J.
  const double doublon = std::hypot(0.08, 0.1);
  const std::vector<double> atom{0, -0.54, -0.54, -0.46, -0.46, -0.7, -0.7, -0.7,
                                -0.5, -0.4 - doublon, -0.4 + doublon, -0.24, -0.24, -0.16, -0.16, 0.6};
  const auto ground = *std::min_element(atom.begin(), atom.end());
  close(result.ground_energy, ground, 2e-12, "Kanamori ground state");

  std::ifstream data_file(result.data_file);
  std::ostringstream numeric;
  for (std::string line; std::getline(data_file, line);) {
    if (line.starts_with("f")) break;
    if (!line.starts_with('#')) numeric << line << '\n';
  }
  std::istringstream data(numeric.str());
  unsigned channels, chain_size, blocks;
  check(static_cast<bool>(data >> channels >> chain_size >> blocks) && channels == 2 && chain_size == 1,
        "Malformed data header.");
  std::vector<double> actual, expected;
  for (unsigned b = 0; b < blocks; ++b) {
    int charge;
    unsigned spin_multiplicity, levels;
    check(static_cast<bool>(data >> charge >> spin_multiplicity >> levels), "Malformed atomic sector.");
    for (unsigned i = 0; i < levels; ++i) {
      double energy;
      check(static_cast<bool>(data >> energy), "Missing atomic energy.");
      actual.insert(actual.end(), spin_multiplicity, energy);
    }
  }
  for (const auto energy : atom) expected.insert(expected.end(), 16, energy - ground);
  std::ranges::sort(actual);
  std::ranges::sort(expected);
  check(actual.size() == expected.size(), "Unexpected dimension of Kanamori atom and uncoupled bath.");
  for (std::size_t i = 0; i < expected.size(); ++i) close(actual[i], expected[i], 2e-12, "atomic spectrum");
}

void interacting_atom(const fs::path &model, const fs::path &work) {
  constexpr double temperature = 0.12;
  constexpr double a = -0.24, b = -0.32, ua = 0.11, ub = 0.63, hybridization = 0.075;
  const auto path = work / "mixed-atom";
  fs::create_directories(path / "scratch");
  {
    std::ofstream param(path / "param");
    param << "[extra]\nU1=" << ua << "\nU2=" << ub << "\nU12=0\nJ12=0\neps1=" << a << "\neps2=" << b
          << "\nV12=" << hybridization << "\n[param]\nsymtype=QS\nLambda=2\nz=1\nNmax=1\nT=" << temperature
          << "\nkeep=10000\nkeepenergy=-1\nbins=4000\nfdm=true\nfdmexpv=true\nsavebins=true\nbroaden=false"
             "\ndata_has_rescaled_energies=false\nspecd=A_d1-A_d2 self_d1-A_d2 A_d1-self_d2 self_d1-self_d2"
             "\nspect=sigma_d1-sigma_d2\n";
  }
  const ChainCoefficients bath{1, 2, 1, {{0, {0.22, 0}, {-0.08, 0.13}}, {0, {0.18, 0}, {0.11, -0.15}}}};
  instantiate({path.string(), "param", model.string(), 1}, bath);
  solve(path.string(), (path / "scratch").string(), "solve.log");

  const auto annihilator = [](unsigned mode) {
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Zero(16, 16);
    const auto bit = 1u << mode;
    for (unsigned state = 0; state < 16; ++state)
      if (state & bit) matrix(state ^ bit, state) = std::popcount(state & (bit - 1)) % 2 ? -1.0 : 1.0;
    return matrix;
  };
  const auto c = annihilator(0), c_down = annihilator(1), f = annihilator(2), f_down = annihilator(3);
  const Eigen::MatrixXd nc = c.transpose() * c, nc_down = c_down.transpose() * c_down;
  const Eigen::MatrixXd nf = f.transpose() * f, nf_down = f_down.transpose() * f_down;
  const Eigen::MatrixXd interaction = ua * nc * nc_down + ub * nf * nf_down;
  const Eigen::MatrixXd h = a * (nc + nc_down) + b * (nf + nf_down) +
    hybridization * (c.transpose() * f + f.transpose() * c + c_down.transpose() * f_down + f_down.transpose() * c_down) + interaction;
  const Eigen::MatrixXd qc = c * interaction - interaction * c, qf = f * interaction - interaction * f;
  const Eigen::MatrixXd sc = 0.5 * (nc - nc_down), sf = 0.5 * (nf - nf_down);
  const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(h);
  check(eig.info() == Eigen::Success, "Could not diagonalize the interacting atom.");
  const auto &e = eig.eigenvalues();
  const auto &v = eig.eigenvectors();
  const Eigen::VectorXd weights = (-(e.array() - e[0]) / temperature).exp();
  const auto partition = weights.sum();
  const Complex z{0.085, 0.27};
  const auto compare = [&](const std::string &kind, const std::string &name, const Eigen::MatrixXd &op1,
                           const Eigen::MatrixXd &op2, bool fermionic) {
    const Eigen::MatrixXd left = v.transpose() * op1 * v, right = v.transpose() * op2 * v;
    Complex expected{}, actual{};
    for (int m = 0; m < 16; ++m)
      for (int n = 0; n < 16; ++n) {
        const auto weight = weights[m] + (fermionic ? weights[n] : -weights[n]);
        expected += weight * left(m, n) * right(m, n) / (z + e[m] - e[n]);
      }
    const auto filename = spectrum_file_name(kind, name, "FDM") + ".bin";
    for (const auto &bin : NRG::Numerics::load_spectrum_bins((path / filename).string()))
      actual += bin.weight / (z - bin.frequency);
    close(actual, expected / partition, 2e-7, name);
  };
  compare("specd", "A_d1-A_d2", c, f, true);
  compare("specd", "self_d1-A_d2", qc, f, true);
  compare("specd", "A_d1-self_d2", c, qf, true);
  compare("specd", "self_d1-self_d2", qc, qf, true);
  compare("spect", "sigma_d1-sigma_d2", sc, sf, false);
}

int main(int argc, char **argv) {
  try {
    check(argc == 3, "expected template and output paths");
    const auto work = fs::absolute(argv[2]);
    fs::remove_all(work);
    run(fs::absolute(argv[1]), work / "diagonal", false, 1.0);
    run(fs::absolute(argv[1]), work / "dense", true, 1.0);
    run(fs::absolute(argv[1]), work / "diagonal-rescaled", false, 4.0);
    run(fs::absolute(argv[1]), work / "dense-rescaled", true, 4.0);
    kanamori_atom(fs::absolute(argv[1]), work);
    interacting_atom(fs::absolute(argv[1]), work);
    const FdmMatrix green{2., 0.3, 0.3, 1.}, left{1., 3., 2., 0.}, right{1., 2., 3., 0.}, zero(4);
    const auto sigma = symmetric_self_energy(2, zero, green, left, right, zero);
    close(sigma[1], -(2.0 - 6.0 * 0.3) / 1.91, 1e-13, "ordered matrix estimator");
    return 0;
  } catch (const std::exception &e) {
    std::cerr << e.what() << '\n';
    return 1;
  }
}
