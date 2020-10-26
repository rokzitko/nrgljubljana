#include <gtest/gtest.h>
#include <mk_sym.hpp>
#include <symmetry.hpp>
#include <params.hpp>
#include <outfield.hpp>

using namespace NRG;

TEST(QS, mk_sym) {
  Workdir workdir("test_workdir");
  Params P("", "param", workdir, true);
  TD td(P, "td");
  auto sym = get<double>("QS", P, td.allfields);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
