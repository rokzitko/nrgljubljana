#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include "test_clean.hpp"
#include <operators.hpp>

using namespace NRG;

TEST(Operators, MatrixElements) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag<double>(P, Sym);
  {
    std::string str =
      "1\n"
      "0 1 0 1\n"
      "1 2\n 3 4\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "2\n"
      "0 1 0 1\n"
      "1 2\n 3 4\n"
      "1 2 1 2\n"
      "1 2 3\n 4 5 6\n 7 8 9\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "1\n"
      "0 1 1 2\n"
      "1 2 3\n 4 5 6\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "1\n"
      "0 1 0 99\n"
      "1 2\n 3 4\n";
    std::istringstream ss(str);
    EXPECT_THROW(MatrixElements(ss, diag), std::runtime_error);
  }
}

TEST(Operators, Opch_empty) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag<double>(P, Sym);
  {
    auto o = Opch<double>(P);
    o.dump();
    EXPECT_EQ(o.size(), 1); // channels
    EXPECT_EQ(o[0].size(), 1); // perchannel
    EXPECT_EQ(o[0][0].size(), 0);
  }
}

TEST(Operators, Opch_clean) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag_clean<double>(P, Sym);
  diag.states_report(Sym->multfnc());
  std::string str =
    "f 0 0\n"
    "2\n"
    "1 1 0 2\n"
    "1.4142135623730951\n"
    "0 2 -1 1\n"
    "1.\n";
  std::istringstream ss(str);
  auto o = Opch<double>(ss, diag, P);
  o.dump();
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
