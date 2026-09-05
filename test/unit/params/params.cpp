#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <list>
#include <limits>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <vector>

#include <gtest/gtest.h>

#include "test_common.hpp"
#include <h5.hpp>
#include <params.hpp>

using namespace std::string_literals;
using namespace NRG;

TEST(params, parser) {
  {
    std::list<parambase *> all;
    param<double> p{"p", "Testing parameter", "1.0", all};
    EXPECT_EQ(p.getkeyword(), "p"s);
    EXPECT_EQ(p.getdesc(), "Testing parameter"s);
    std::ostringstream ss;
    p.dump(ss);
    EXPECT_EQ(ss.str(), "p=1\n"s);
    EXPECT_EQ(p, 1.0);
    EXPECT_EQ(p.value(), 1.0);

    p.set_str("2.0");
    std::ostringstream ss2;
    p.dump(ss2);
    EXPECT_EQ(ss2.str(), "p=2 *\n"s);
    EXPECT_EQ(p, 2.0);
    EXPECT_EQ(p.value(), 2.0);

    p = 3.0;
    EXPECT_EQ(p, 3.0);
    EXPECT_EQ(p.value(), 3.0);

    EXPECT_EQ(p == 3.0, true);
    EXPECT_EQ(p == 4.0, false);

    EXPECT_THROW(param<double>("p", "another one", "1.0", all), std::runtime_error);
  }
}

TEST(params, DiagParams) { // NOLINT
  Params P;
  auto DP = DiagParams(P, 0.5);
}

TEST(params, Defaults) {
  // check important default values
  Params P;
  EXPECT_EQ(P.symtype, ""s);
  EXPECT_EQ(P.Lambda, 2.0);
  EXPECT_EQ(P.z, 1.0);
  EXPECT_EQ(P.bandrescale, 1.0);
  EXPECT_EQ(P.diag, "default"s);
  EXPECT_EQ(P.keepenergy, -1.0);
  EXPECT_EQ(P.keepmin, 0ul);
  EXPECT_EQ(P.fixeps, 1e-15);
  EXPECT_EQ(P.dm, false);
  EXPECT_EQ(P.fdm_cutoff, 1e-16);
  EXPECT_EQ(P.strategy, "kept");
  EXPECT_EQ(P.clip_tol_imag, 1e-10);
  EXPECT_EQ(P.bins, 1000ul);
  EXPECT_EQ(P.discard_trim, 1e-16);
  EXPECT_EQ(P.discard_immediately, 1e-16);
  EXPECT_EQ(P.prec_td, 10);
  EXPECT_EQ(P.prec_custom, 10);
  EXPECT_EQ(P.prec_xy, 10);
  EXPECT_EQ(P.done, true);
  EXPECT_EQ(P.calc0, true);
  EXPECT_EQ(P.lastall, false);
  EXPECT_EQ(P.lastalloverride, false);

  EXPECT_EQ(P.Ninit.value(), 0);
  EXPECT_EQ(P.Nmax, 0ul);
}

TEST(params, set_channels_and_combs) {
  Params P;
  EXPECT_EQ(P.symtype, ""s);
  P.set_channels_and_combs(1);
  EXPECT_EQ(P.channels, 1);
  EXPECT_EQ(P.coeffactor, 1);
  EXPECT_EQ(P.coefchannels, 1);
  EXPECT_EQ(P.perchannel, 1);
  EXPECT_EQ(P.spin, 2);
  EXPECT_EQ(P.combs, 4);
}

template< class InputIt, class UnaryPredicate>
bool one_of(InputIt first, InputIt last, UnaryPredicate p)
{
  return std::count(first, last, p) == 1;
}

template< class InputIt>
bool one_of(InputIt first, InputIt last)
{
  return std::count(first, last, true) == 1;
}

void check_recalc(const Params &P, const RUNTYPE runtype)
{
  std::vector v = { P.do_recalc_kept(runtype), P.do_recalc_all(runtype) };
  EXPECT_EQ(one_of(v.begin(), v.end()), true);
}

void check_recalc(const Params &P)
{
  check_recalc(P, RUNTYPE::NRG);
  check_recalc(P, RUNTYPE::DMNRG);
}

TEST(params, flags_none) {
  Params P;
  P.Nmax = 1;
  EXPECT_EQ(P.cfs_flags(), false);
  EXPECT_EQ(P.fdm_flags(), false);
  EXPECT_EQ(P.dmnrg_flags(), false);
  EXPECT_EQ(P.cfs_or_fdm_flags(), false);
  EXPECT_EQ(P.dm_flags(), false);
  EXPECT_EQ(P.keep_all_states_in_last_step(), false);
  EXPECT_EQ(P.need_rho(), false);
  EXPECT_EQ(P.need_rhoFDM(), false);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::NRG), true);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::DMNRG), true);
  check_recalc(P);
}

TEST(params, flags_cfs) {
  Params P;
  P.Nmax = 1;
  P.cfs = true;
  EXPECT_EQ(P.cfs_flags(), true);
  EXPECT_EQ(P.fdm_flags(), false);
  EXPECT_EQ(P.dmnrg_flags(), false);
  EXPECT_EQ(P.cfs_or_fdm_flags(), true);
  EXPECT_EQ(P.dm_flags(), true);
  EXPECT_EQ(P.keep_all_states_in_last_step(), true);
  EXPECT_EQ(P.need_rho(), true);
  EXPECT_EQ(P.need_rhoFDM(), false);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::NRG), true);
  EXPECT_EQ(P.do_recalc_all(RUNTYPE::DMNRG), true);
  check_recalc(P);
}

TEST(params, flags_fdm) {
  Params P;
  P.Nmax = 1;
  P.fdm = true;
  EXPECT_EQ(P.cfs_flags(), false);
  EXPECT_EQ(P.fdm_flags(), true);
  EXPECT_EQ(P.dmnrg_flags(), false);
  EXPECT_EQ(P.cfs_or_fdm_flags(), true);
  EXPECT_EQ(P.dm_flags(), true);
  EXPECT_EQ(P.keep_all_states_in_last_step(), true);
  EXPECT_EQ(P.need_rho(), false);
  EXPECT_EQ(P.need_rhoFDM(), true);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::NRG), true);
  EXPECT_EQ(P.do_recalc_all(RUNTYPE::DMNRG), true);
  check_recalc(P);
}

TEST(params, flags_fdmexpv) {
  Params P;
  P.Nmax = 1;
  P.fdmexpv = true;
  EXPECT_EQ(P.cfs_flags(), false);
  EXPECT_EQ(P.fdm_flags(), true);
  EXPECT_EQ(P.dmnrg_flags(), false);
  EXPECT_EQ(P.cfs_or_fdm_flags(), true);
  EXPECT_EQ(P.dm_flags(), true);
  EXPECT_EQ(P.keep_all_states_in_last_step(), true);
  EXPECT_EQ(P.need_rho(), false);
  EXPECT_EQ(P.need_rhoFDM(), true);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::NRG), true);
  EXPECT_EQ(P.do_recalc_all(RUNTYPE::DMNRG), true);
  check_recalc(P);
}

TEST(params, flags_dmnrg) {
  Params P;
  P.Nmax = 1;
  P.dmnrg = true;
  EXPECT_EQ(P.cfs_flags(), false);
  EXPECT_EQ(P.fdm_flags(), false);
  EXPECT_EQ(P.dmnrg_flags(), true);
  EXPECT_EQ(P.cfs_or_fdm_flags(), false);
  EXPECT_EQ(P.dm_flags(), true);
  EXPECT_EQ(P.keep_all_states_in_last_step(), false);
  EXPECT_EQ(P.need_rho(), true);
  EXPECT_EQ(P.need_rhoFDM(), false);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::NRG), true);
  EXPECT_EQ(P.do_recalc_kept(RUNTYPE::DMNRG), true);
  check_recalc(P);
}

TEST(params, scale) {
  Params P;
  EXPECT_EQ(P.Lambda, 2.0);
  EXPECT_EQ(P.discretization, "Z"s);
  const double s1 = (1.0-1.0/2.0)/std::log(2.0);
  EXPECT_LT(std::abs(P.SCALE(1) - s1), 1e-10);
  const double s2 = s1/std::sqrt(2.0);
  EXPECT_LT(std::abs(P.SCALE(2) - s2), 1e-10);

  EXPECT_LT(std::abs(P.nrg_step_scale_factor() - std::sqrt(2.0)), 1e-10);

  P.Nmax = 1;
  EXPECT_LT(std::abs(P.last_step_scale() - s1), 1e-10);

  P.absolute = true;
  EXPECT_EQ(P.nrg_step_scale_factor(), 1.0);
}

TEST(params, invalid_discretization) {
  Params P;
  P.discretization = "bad"s;
  EXPECT_THROW(P.validate(), std::invalid_argument);
}

TEST(params, E) {
  Params P;
  P.Nmax = 1;
  EXPECT_EQ(P.getEfactor(), std::sqrt(2.0));
  EXPECT_EQ(P.getE0(), 2.0);
  EXPECT_EQ(P.getEmin(), 2.0);
  EXPECT_EQ(P.getEx(), 2.0*std::sqrt(2.0));
  EXPECT_LT(std::abs(P.getEmax()-4.0), 1e-10);
}

TEST(params, validate_rejects_invalid_nmax) {
  Params P;

  EXPECT_THROW(P.validate_after_data_file(), std::invalid_argument);

  P.Ninit = 1;
  P.Nmax = 1;
  EXPECT_THROW(P.validate_after_data_file(), std::invalid_argument);

  P.Ninit = 0;
  P.Nmax = 2;
  EXPECT_NO_THROW(P.validate_after_data_file());
}

TEST(params, validate_rejects_non_positive_temperature) {
  Params P;
  P.Nmax = 2;

  P.T = 0.0;
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.T = -1.0;
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.T = 1e-3;
  EXPECT_NO_THROW(P.validate());
}

TEST(params, validate_rejects_invalid_thresholds) {
  Params P;

  P.fdm_cutoff = -1.0;
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.fdm_cutoff = 0.0;
  P.clip_tol_imag = -1.0;
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.clip_tol_imag = 0.0;
  P.fdm_cutoff = std::numeric_limits<double>::infinity();
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.fdm_cutoff = 0.0;
  P.clip_tol_imag = std::numeric_limits<double>::infinity();
  EXPECT_THROW(P.validate(), std::invalid_argument);

  P.clip_tol_imag = 0.0;
  EXPECT_NO_THROW(P.validate());
}

TEST(params, validate_rejects_partial_diagonalisation_for_floquet) {
  for (const auto &diag : {"dsyevr"s, "zheevr"s}) {
    SCOPED_TRACE(diag);
    Params P;
    P.diag = diag;
    P.diagratio = 0.5;
    P.floquet = true;
    P.extra_params["Omega"] = "1";
    EXPECT_THROW(P.validate(), std::invalid_argument);

    P.diagratio = 1.0;
    EXPECT_NO_THROW(P.validate());

    P.diagratio = 0.5;
    P.floquet = false;
    EXPECT_NO_THROW(P.validate());
  }
}

TEST(params, validate_floquet_omega) {
  Params P;
  P.floquet = true;

  for (const auto &value : {"1", "0.25", "+1e-3", "2E+2", "1e-300"}) {
    SCOPED_TRACE(value);
    P.extra_params["Omega"] = value;
    EXPECT_NO_THROW(P.validate());
    EXPECT_GT(P.floquet_omega(), 0.0);
    EXPECT_TRUE(std::isfinite(P.floquet_omega()));
  }

  P.extra_params["Omega"] = "1e-32";
  P.validate();
  EXPECT_EQ(P.floquet_omega(), 1e-32);

  for (const auto &value : {"", "0", "-0", "-1", "1junk", "1 2", "nan", "inf", "1e9999", "1e-9999",
                            "1.7976931348623159e308", "2e-324", "2.4703282292062327e-324"}) {
    SCOPED_TRACE(value);
    P.extra_params["Omega"] = value;
    EXPECT_THROW(P.validate(), std::invalid_argument);
  }

  P.extra_params.clear();
  EXPECT_THROW(P.validate(), std::invalid_argument);
}

TEST(params, ignores_floquet_omega_when_disabled) {
  Params P;
  P.extra_params["Omega"] = "not-a-number";
  EXPECT_NO_THROW(P.validate());
  EXPECT_THROW(static_cast<void>(P.floquet_omega()), std::logic_error);
}

TEST(params, floquet_keeps_existing_algorithm_paths_enabled) {
  Params P;
  P.floquet = true;
  P.extra_params["Omega"] = "0.5";
  P.dm = true;
  P.finite = true;
  P.finitemats = true;
  P.dmnrg = true;
  P.dmnrgmats = true;
  P.cfs = true;
  P.cfsgt = true;
  P.cfsls = true;
  P.fdm = true;
  P.fdmgt = true;
  P.fdmls = true;
  P.fdmmats = true;
  P.fdmexpv = true;
  P.specs = "m-m";
  P.specd = "d-d";
  P.spect = "s-s";
  P.specq = "q-q";
  P.specot = "o-o";

  EXPECT_NO_THROW(P.validate());
}

TEST(params, floquet_runtime_metadata_is_unregistered) {
  Params P;
  P.floquet = true;
  P.extra_params["Omega"] = "0.5";
  P.validate();
  P.set_floquet_mode_bounds({-2.0, 2.0});
  EXPECT_EQ(P.floquet_mode_bounds(), std::make_pair(-2.0, 2.0));

  std::ostringstream output;
  P.dump(output);
  EXPECT_EQ(output.str().find("Omega"), std::string::npos);
  EXPECT_EQ(output.str().find("mode_bounds"), std::string::npos);
}

TEST(params, h5save_stores_nlen) {
  Params P;
  P.Nmax = 7;
  P.Nlen = 3;

  H5Easy::File file("params.h5", H5Easy::File::Overwrite);
  P.h5save(file);

  EXPECT_EQ(H5Easy::load<std::vector<size_t>>(file, "params/Nmax").front(), 7UL);
  EXPECT_EQ(H5Easy::load<std::vector<size_t>>(file, "params/Nlen").front(), 3UL);

  std::remove("params.h5");
}

TEST(params, parser_accepts_extended_bool_values) {
  const auto filename = "params_bool.param";
  {
    std::ofstream file(filename);
    file << "  [param]\n"
         << "  done = no\n"
         << "calc0 = YES\n"
         << " lastall = 1\n"
         << "fdm = 0\n";
  }

  Params P(filename, "param", std::make_unique<Workdir>(".", true), true, true);
  EXPECT_FALSE(P.done);
  EXPECT_TRUE(P.calc0);
  EXPECT_TRUE(P.lastall);
  EXPECT_FALSE(P.fdm);

  std::remove(filename);
}

TEST(params, parser_accepts_threshold_values) {
  const auto filename = "params_thresholds.param";
  {
    std::ofstream file(filename);
    file << "[param]\n"
         << "fdm_cutoff = 1e-12\n"
         << "clip_tol_imag = 1e-8\n";
  }

  Params P(filename, "param", std::make_unique<Workdir>(".", true), true, true);
  EXPECT_EQ(P.fdm_cutoff, 1e-12);
  EXPECT_EQ(P.clip_tol_imag, 1e-8);

  std::remove(filename);
}

TEST(params, parser_rejects_invalid_bool_values) {
  const auto filename = "params_invalid_bool.param";
  {
    std::ofstream file(filename);
    file << "[param]\n"
         << "done = maybe\n";
  }

  EXPECT_THROW(Params(filename, "param", std::make_unique<Workdir>(".", true), true, true), std::runtime_error);
  std::remove(filename);
}

TEST(params, resume_requires_persistent_workdir) {
  Params P;
  P.resume = true;
  P.Nmax = 2;
  EXPECT_THROW(P.init_laststored(), std::invalid_argument);
}

TEST(params, persistent_workdir_requires_resume) {
  Workdir parent(".", true);
  Params P;
  P.workdir = std::make_unique<Workdir>(parent.get() + "/checkpoint", WorkdirMode::persistent_exact, true);

  EXPECT_THROW(P.init_laststored(), std::invalid_argument);
}

TEST(params, resume_finds_contiguous_checkpoint_prefix) {
  Workdir parent(".", true);
  Params P;
  P.workdir = std::make_unique<Workdir>(parent.get() + "/checkpoint", WorkdirMode::persistent_exact, true);
  P.resume = true;
  P.Ninit = 1;
  P.Nmax = 4;
  std::ofstream(P.workdir->unitaryfn(1), std::ios::binary) << "checkpoint";
  std::ofstream(P.workdir->unitaryfn(2), std::ios::binary) << "checkpoint";

  P.init_laststored();

  ASSERT_TRUE(P.laststored.has_value());
  EXPECT_EQ(P.laststored.value(), 2U);
  EXPECT_TRUE(P.resume_iteration(1));
  EXPECT_TRUE(P.resume_iteration(2));
  EXPECT_FALSE(P.resume_iteration(3));
}

TEST(params, resume_rejects_checkpoint_gaps) {
  Workdir parent(".", true);
  Params P;
  P.workdir = std::make_unique<Workdir>(parent.get() + "/checkpoint", WorkdirMode::persistent_exact, true);
  P.resume = true;
  P.Nmax = 3;
  std::ofstream(P.workdir->unitaryfn(0), std::ios::binary) << "checkpoint";
  std::ofstream(P.workdir->unitaryfn(2), std::ios::binary) << "checkpoint";

  EXPECT_THROW(P.init_laststored(), std::runtime_error);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
