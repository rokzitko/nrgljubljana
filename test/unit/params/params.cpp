#include <string>
#include <sstream>
#include <cstdio>
#include <fstream>

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
  EXPECT_EQ(P.strategy, "kept");
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
  const double s2 = s1/sqrt(2.0);
  EXPECT_LT(std::abs(P.SCALE(2) - s2), 1e-10);

  EXPECT_LT(std::abs(P.nrg_step_scale_factor() - sqrt(2.0)), 1e-10);

  P.Nmax = 1;
  EXPECT_LT(std::abs(P.last_step_scale() - s1), 1e-10);

  P.absolute = true;
  EXPECT_EQ(P.nrg_step_scale_factor(), 1.0);
}

TEST(params, E) {
  Params P;
  P.Nmax = 1;
  EXPECT_EQ(P.getEfactor(), sqrt(2.0));
  EXPECT_EQ(P.getE0(), 2.0);
  EXPECT_EQ(P.getEmin(), 2.0);
  EXPECT_EQ(P.getEx(), 2.0*sqrt(2.0));
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

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
