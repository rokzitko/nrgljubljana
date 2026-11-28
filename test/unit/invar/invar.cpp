#include <gtest/gtest.h>

#include <initializer_list>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <invar.hpp>
#include <mk_sym.hpp>
#include "test_common.hpp"

TEST(Invar, initInvar) { // NOLINT
  Params P;
  auto Sym = setup_Sym<double>(P);
  // should have Q and SS degrees of freedom
  EXPECT_EQ(Invar::invdim, 2);
  EXPECT_EQ(Invar::qntype.size(), 2);
  EXPECT_EQ(Invar::qntype[0], additive);
  EXPECT_EQ(Invar::qntype[1], additive);
  EXPECT_EQ(Invar::names.size(), 2);
  EXPECT_EQ(Invar::names["Q"], 0);
  EXPECT_EQ(Invar::names["SS"], 1);
}

TEST(Invar, initInvar_resets_global_metadata) { // NOLINT
  Params P1;
  auto Sym1 = set_symmetry<double>(P1, "QSZ", 1);
  ASSERT_EQ(Invar::qntype.size(), 2U);
  ASSERT_EQ(Invar::names.size(), 2U);
  EXPECT_TRUE(Invar::names.count("Q"));
  EXPECT_TRUE(Invar::names.count("SSZ"));

  Params P2;
  auto Sym2 = set_symmetry<double>(P2, "QS", 1);
  ASSERT_EQ(Invar::qntype.size(), 2U);
  ASSERT_EQ(Invar::names.size(), 2U);
  EXPECT_TRUE(Invar::names.count("Q"));
  EXPECT_TRUE(Invar::names.count("SS"));
  EXPECT_FALSE(Invar::names.count("SSZ"));
}

TEST(Invar, InvarQS) { // NOLINT
  using std::literals::string_literals::operator""s;

  Params P;
  auto Sym = setup_Sym<double>(P);
  {
    Invar I(1,2);
  }
  {
    std::string str = "1 2";
    std::istringstream ss(str);
    Invar I;
    ss >> I;
    EXPECT_EQ(I, Invar(1,2));
  }
  {
    std::ostringstream ss;
    Invar I(1,2);
    ss << I;
    EXPECT_EQ(ss.str(), "1 2"s);
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.str(), "1 2"s);
    EXPECT_EQ(I.name(), "1_2"s);
  }
  {
    Invar I1(1);
    EXPECT_EQ(I1.str(), "1"s);
    Invar I2(1,2);
    EXPECT_EQ(I2.str(), "1 2"s);
    Invar I3(1,2,3);
    EXPECT_EQ(I3.str(), "1 2 3"s);
  }
  {
    Invar I1(1,1);
    Invar I2(1,2);
    EXPECT_EQ(I1 == I1, true);
    EXPECT_EQ(I1 != I2, true);
    EXPECT_EQ(I1 < I2, true);
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.getqn(0), 1);
    EXPECT_EQ(I.getqn(1), 2);
  }
  {
    Invar I1(2,2);
    Invar I2(3,2);
    I1.combine(I2);
    EXPECT_EQ(I1, Invar(5, 4));
  }
  {
    Invar I(1,2);
    I.inverse();
    EXPECT_EQ(I, Invar(-1,-2)); // both additive!
  }
  {
    Invar I(1,2);
    EXPECT_EQ(I.get("Q"), 1);
    EXPECT_EQ(I.get("SS"), 2);
  }
}

TEST(Invar, InvertParity) { // NOLINT
  initInvar({
    {"Q", additive},
    {"SS", additive},
    {"P", multiplicative}
  });
  const Invar I(1,2,1);
  const auto inverted = I.InvertParity();
  EXPECT_EQ(inverted, Invar(1,2,-1));
  EXPECT_EQ(I, Invar(1,2,1));
  EXPECT_EQ(inverted.InvertParity(), I);
}

// Checkpoints and MPI messages store Invar in the format of a class with a single std::vector<int> member, which was
// its original storage. The archive format must not depend on the in-memory storage of the quantum numbers.
namespace {
struct VectorInvar {
  std::vector<int> data;
  template <class Archive> void serialize(Archive &ar, [[maybe_unused]] const unsigned int version) { ar &data; }
};
} // namespace

TEST(Invar, serialization_format) { // NOLINT
  initInvar({{"Q", additive}, {"SS", additive}, {"P", multiplicative}});
  const Invar I(-3, 4, -1);
  const VectorInvar reference{{-3, 4, -1}};
  std::ostringstream s1, s2;
  {
    boost::archive::binary_oarchive oa(s1);
    oa << I;
    boost::archive::binary_oarchive ob(s2);
    ob << reference;
  }
  EXPECT_EQ(s1.str(), s2.str());
  std::istringstream in(s2.str());
  boost::archive::binary_iarchive ia(in);
  Invar J;
  ia >> J;
  EXPECT_EQ(J, I);
}

TEST(Invar, comparison_is_lexicographic) { // NOLINT
  initInvar({{"Q", additive}, {"SS", additive}});
  EXPECT_LT(Invar(-1, 5), Invar(0, 1));
  EXPECT_LT(Invar(0, 1), Invar(0, 2));
  EXPECT_FALSE(Invar(0, 2) < Invar(0, 2));
  EXPECT_LT(Invar(1), Invar(1, 0)); // a shorter prefix orders first, as for std::vector
  EXPECT_NE(Invar(1), Invar(1, 0));
}

TEST(Invar, get_unknown_name_throws) { // NOLINT
  initInvar({{"Q", additive}, {"SS", additive}});
  EXPECT_THROW(static_cast<void>(Invar(0, 1).get("SSZ")), std::invalid_argument);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
