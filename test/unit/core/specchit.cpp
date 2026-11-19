#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <core.hpp>

namespace {

template<typename S>
class SpecChit : public ::testing::Test {
  const std::filesystem::path previous = std::filesystem::current_path();
  Workdir directory{std::filesystem::temp_directory_path().string()};
 public:
  SpecChit() { std::filesystem::current_path(directory.get()); }
  ~SpecChit() override { std::filesystem::current_path(previous); }
};

using ScalarTypes = ::testing::Types<double, std::complex<double>>;
TYPED_TEST_SUITE(SpecChit, ScalarTypes);

Invar spin_sector(const int ss, const int orbital_t = 0) {
  InvType qns(Invar::invdim, 0);
  for (size_t i = 0; i < qns.size(); ++i)
    if (Invar::qntype[i] == multiplicative) qns[i] = 1;
  qns[Invar::names.at("SS")] = ss;
  if (Invar::names.contains("T")) qns[Invar::names.at("T")] = orbital_t;
  return Invar(qns);
}

void check_chit_output(const std::string &name, const std::vector<std::pair<double, double>> &expected) {
  const auto filename = "chit_CHIT_dens_" + name + ".dat";
  SCOPED_TRACE(filename);
  std::ifstream input(filename);
  for (const auto &[temperature, value] : expected) {
    double t{}, real{}, imaginary{};
    ASSERT_TRUE(input >> t >> real >> imaginary);
    EXPECT_NEAR(t, temperature, 1e-13);
    EXPECT_NEAR(real, value, 1e-13);
    EXPECT_NEAR(imaginary, 0.0, 1e-13);
  }
  std::string extra;
  EXPECT_FALSE(input >> extra);
  input.close();
  EXPECT_TRUE(std::filesystem::remove(filename)); // no stale output between symmetry cases
}

TYPED_TEST(SpecChit, DispatchAndCurieNormalization) { // NOLINT
  std::vector<std::pair<std::string, size_t>> symmetries{{"QS", 1}};
#ifdef NRG_SYM_MORE
  symmetries.emplace_back("SPSU2", 1);
#endif
#ifdef NRG_SYM_ALL
  for (const auto *name : {"QSLR", "SPSU2LR"}) symmetries.emplace_back(name, 2);
  for (const auto *name : {"QST", "QSTZ", "SPSU2T"}) symmetries.emplace_back(name, 3);
  if constexpr (is_complex<TypeParam>::value)
    for (const auto *name : {"QSC3", "SPSU2C3"}) symmetries.emplace_back(name, 3);
#endif
  for (const auto &[symmetry, channels] : symmetries) {
    SCOPED_TRACE(symmetry);
    for (const int ss : {2, 3}) for (const int orbital_t : {0, 1, 2}) {
      if (orbital_t != 0 && symmetry != "QST" && symmetry != "SPSU2T") continue;
      SCOPED_TRACE(ss);
      SCOPED_TRACE(orbital_t);
      Params P;
      P.absolute = true;
      P.reim = true;
      P.prec_xy = 17;
      P.specchit = "spin-spin local-local local-global global-local global-global";
      auto Sym = set_symmetry<TypeParam>(P, symmetry, channels);
      const auto I = spin_sector(ss, orbital_t);
      const auto degeneracy = ss * (2 * orbital_t + 1);
      EXPECT_EQ(Sym->mult(I), degeneracy);
      const double spin = (ss - 1) / 2.0;
      DiagInfo<TypeParam> diag;
      diag[I] = NRG::Eigen<TypeParam>(std::vector<double>{0.0}, 1.0);
      Operators<TypeParam> a;
      // <S,m|Sz|S,m> = m: the reduced element is sqrt(S(S+1)), not a Pauli matrix.
      a.opt["spin"][{I, I}] = std::sqrt(spin * (spin + 1)) * id_matrix<TypeParam>(1);
      a.opt["unused"] = {};
      a.ops["local"][{I, I}] = 2.0 * id_matrix<TypeParam>(1);
      a.opsg["global"][{I, I}] = 3.0 * id_matrix<TypeParam>(1);
      MemTime mt;
      Oprecalc<TypeParam> spectra(RUNTYPE::NRG, a, Sym, mt, P);
      ASSERT_EQ(spectra.sl.size(), 5);
      EXPECT_EQ(spectra.ops.size(), 3);
      EXPECT_EQ(spectra.ops.count({"t", "spin"}), 1);
      EXPECT_EQ(spectra.ops.count({"t", "unused"}), 0);
      for (const auto &spectrum : spectra.sl) {
        EXPECT_NE(dynamic_cast<Algo_CHIT<TypeParam> *>(spectrum.algo.get()), nullptr);
        EXPECT_EQ(spectrum.spin, 0);
        const bool triplet = &spectrum.op1 == &a.opt.at("spin");
        EXPECT_DOUBLE_EQ(spectrum.ff(I, I), triplet ? degeneracy / 3.0 : degeneracy);
      }
      Stats<TypeParam> stats(P, Sym->get_td_fields(), 0.0);
      DensMatElements<TypeParam> rho;
      const Step step(P);
      for (const double temperature : {0.2, 0.7}) {
        P.chitp = temperature;
        calc_Z(step, stats, diag, Sym->multfnc(), P);
        EXPECT_DOUBLE_EQ(stats.Zchit, degeneracy);
        spectra.sl.calc(step, diag, rho, rho, stats, mt, Sym.get(), P);
      }
      spectra.sl.save();
      const double curie = ss == 2 ? 0.25 : 2.0 / 3.0;
      check_chit_output("spin-spin", {{0.2, curie}, {0.7, curie}});
      check_chit_output("local-local", {{0.2, 4.0}, {0.7, 4.0}});
      check_chit_output("local-global", {{0.2, 6.0}, {0.7, 6.0}});
      check_chit_output("global-local", {{0.2, 6.0}, {0.7, 6.0}});
      check_chit_output("global-global", {{0.2, 9.0}, {0.7, 9.0}});
      Oprecalc<TypeParam> dm(RUNTYPE::DMNRG, a, Sym, mt, P);
      EXPECT_TRUE(dm.sl.empty());
      EXPECT_TRUE(dm.ops.empty());
    }
  }
}

TYPED_TEST(SpecChit, NonconservedLocalSpinLehmann) { // NOLINT
  // H=J S1.S2: a singlet at 0 and a triplet at J. Sz1 connects s to t0,
  // and has diagonal elements +/-1/2 in t+/t-. Thus T*chi is not <Sz1^2>=1/4.
  const double J = 0.43;
  std::vector<std::pair<std::string, size_t>> symmetries{{"QS", 1}};
#ifdef NRG_SYM_ALL
  symmetries.emplace_back("QST", 3);
  symmetries.emplace_back("SPSU2T", 3);
#endif
  for (const auto &[symmetry, channels] : symmetries) for (const bool absolute : {false, true}) {
    SCOPED_TRACE(symmetry);
    SCOPED_TRACE(absolute);
    Params P;
    P.absolute = absolute;
    P.Lambda = 2.0;
    P.reim = true;
    P.prec_xy = 17;
    P.specchit = "spin-spin";
    auto Sym = set_symmetry<TypeParam>(P, symmetry, channels);
    const Step step(P);
    const auto orbital_t = symmetry == "QS" ? 0 : 1;
    const auto singlet = spin_sector(1, orbital_t), triplet = spin_sector(3, orbital_t);
    DiagInfo<TypeParam> diag;
    diag[singlet] = NRG::Eigen<TypeParam>(std::vector<double>{0.0}, step.scale());
    diag[triplet] = NRG::Eigen<TypeParam>(std::vector<double>{J / step.scale()}, step.scale());
    TypeParam phase{1.0};
    if constexpr (is_complex<TypeParam>::value) phase = std::polar(1.0, 0.61);
    Operators<TypeParam> a;
    a.opt["spin"][{triplet, singlet}] = (phase / 2.0) * id_matrix<TypeParam>(1);
    a.opt["spin"][{singlet, triplet}] = (-std::sqrt(3.0) * conj_me(phase) / 2.0) * id_matrix<TypeParam>(1);
    a.opt["spin"][{triplet, triplet}] = (1.0 / std::sqrt(2.0)) * id_matrix<TypeParam>(1);
    MemTime mt;
    Oprecalc<TypeParam> spectra(RUNTYPE::NRG, a, Sym, mt, P);
    ASSERT_EQ(spectra.sl.size(), 1);
    Stats<TypeParam> stats(P, Sym->get_td_fields(), 0.0);
    DensMatElements<TypeParam> rho;
    std::vector<std::pair<double, double>> expected;
    for (const double temperature : {0.17, 0.63}) {
      P.chitp = temperature / step.scale();
      calc_Z(step, stats, diag, Sym->multfnc(), P);
      const auto q = std::exp(-J / temperature);
      EXPECT_NEAR(stats.Zchit, (2 * orbital_t + 1) * (1.0 + 3.0 * q), 1e-13);
      // Full four-state Lehmann sum, independent of reduced elements and chit_weight().
      const auto value = (q / 2.0 + temperature * (1.0 - q) / (2.0 * J)) / (1.0 + 3.0 * q);
      EXPECT_GT(std::abs(value - 0.25), 1e-3);
      expected.emplace_back(temperature, value);
      spectra.sl.calc(step, diag, rho, rho, stats, mt, Sym.get(), P);
    }
    spectra.sl.save();
    check_chit_output("spin-spin", expected);
  }
}

} // namespace
