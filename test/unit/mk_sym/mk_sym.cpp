#include <gtest/gtest.h>
#include <complex>
#include <mk_sym.hpp>
#include <symmetry.hpp>
#include <params.hpp>
#include <outfield.hpp>
#include <symmetry_traits.hpp>

#include <stdexcept>
#include <string>

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

// symmetry_traits lists every symmetry type of get(), with its isfield(). Types not compiled in
// (SYM_MORE/SYM_ALL off, or C3 types for real S) throw from get() and are skipped.
template <typename S> int check_symmetry_traits() {
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    Params P;
    try {
      auto sym = get<S>(std::string(t.name), P);
      EXPECT_EQ(sym->isfield(), t.field) << t.name;
      ++checked;
    } catch (const std::runtime_error &) {}
  }
  return checked;
}

TEST(mk_sym, symmetry_traits) {
  EXPECT_GE(check_symmetry_traits<double>() + check_symmetry_traits<std::complex<double>>(), 4);
  EXPECT_EQ(find_symmetry_traits("XYZ"), nullptr);
  ASSERT_NE(find_symmetry_traits("QSZ"), nullptr);
  EXPECT_TRUE(find_symmetry_traits("QSZ")->field);
}

TEST(mk_sym, perchannel) {
  EXPECT_EQ(perchannel("QS"), 1u);
  EXPECT_EQ(perchannel("QSZ"), 1u);
  EXPECT_EQ(perchannel("U1"), 2u);
  EXPECT_EQ(perchannel("NONE"), 4u);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
