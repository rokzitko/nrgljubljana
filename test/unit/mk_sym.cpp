#include <gtest/gtest.h>
#include <my_sym.hpp>

using namespace NRG;

TEST(QS, mk_sym) {
  Params P;
  Allfields allfields;
  auto get_sym("QS", P, allfields);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
