#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <coef.hpp>

using namespace NRG;

TEST(Coef, parse) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag_clean<double>(P, Sym);
  Coef<double> coef(P);
  std::string str = 
    "1\n"
    "0.54528747084262258072\n"
    "0.41550946829175445321\n"
    "1\n"
    "0.e-999\n"
    "0.e-998\n";
  std::istringstream ss(str);
  EXPECT_EQ(P.coefchannels, 1);
  coef.xi.read(ss, P.coefchannels);
  coef.zeta.read(ss, P.coefchannels);

  EXPECT_EQ(coef.xi.nr_tabs(), 1);
  EXPECT_EQ(coef.zeta.nr_tabs(), 1);
  EXPECT_EQ(coef.xi.max(0), 1);
  EXPECT_EQ(coef.zeta.max(0), 1);

  EXPECT_EQ(coef.xi(0,0), 0.54528747084262258072);
  EXPECT_EQ(coef.xi(1,0), 0.41550946829175445321); // first index: N, second index: alpha (channel number)
  EXPECT_EQ(coef.zeta(0,0), 0);
  EXPECT_EQ(coef.zeta(1,0), 0);
  
  coef.xi.setvalue(0, 0, 42.);
  EXPECT_EQ(coef.xi(0,0), 42.);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
