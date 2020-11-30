#include <gtest/gtest.h>
#include <complex>
#include <mk_sym.hpp>
#include <symmetry.hpp>
#include <params.hpp>
#include <outfield.hpp>

using namespace NRG;

TEST(mk_sym, QS) {
  Params P;
  auto sym = mk_QS<double>(P);
}

TEST(mk_sym, QSZ) {
  Params P;
  auto sym = mk_QSZ<double>(P);
}

TEST(mk_sym, get) {
  Params P;
  auto sym = get<double>("QS", P);
}

TEST(mk_sym, QS_complex) {
  Params P;
  auto sym = mk_QS<std::complex<double>>(P);
}

TEST(mk_sym, QSZ_complex) {
  Params P;
  auto sym = mk_QSZ<std::complex<double>>(P);
}

TEST(mk_sym, get_complex) {
  Params P;
  auto sym = get<std::complex<double>>("QS", P);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
