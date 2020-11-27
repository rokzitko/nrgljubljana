#include <gtest/gtest.h>

#include <eigen.hpp>

using namespace NRG;

TEST(Eigen, constructor) { // NOLINT
  {
    Eigen<double> e;
    EXPECT_EQ(e.getnrcomputed(), 0);
    EXPECT_EQ(e.getdim(), 0);
  }
  {
    Eigen<double> e(2,3);
    EXPECT_EQ(e.getnrcomputed(), 2);
    EXPECT_EQ(e.getdim(), 3);
  }
}

TEST(Eigen, get) { // NOLINT
  Eigen<double> e(5,8);
  e.truncate_prepare(2);
  EXPECT_EQ(e.getnrcomputed(), 5);
  EXPECT_EQ(e.getdim(), 8);
  EXPECT_EQ(e.getnrpost(), 2);
  EXPECT_EQ(e.getnrall(), 5);
  EXPECT_EQ(e.getnrkept(), 2);
  EXPECT_EQ(e.getnrdiscarded(), 3);
  EXPECT_EQ(e.all(), range0(5));
  EXPECT_EQ(e.kept(), range0(2));
  EXPECT_EQ(e.discarded(), boost::irange(2,5) );
  EXPECT_EQ(e.getnrstored(), 5);
  EXPECT_EQ(e.stored(), range0(5));
}

using EVEC = ublas::vector<double>;

TEST(Eigen, diagonal) { // NOLINT
  Eigen<double> e(3,3);
  EVEC v(3);
  v[0] = 1.0;
  v[1] = 2.0;
  v[2] = 3.0;
  e.diagonal(v);
  EXPECT_EQ(e.value_orig[0], 1.0);
  EXPECT_EQ(e.value_zero[0], 1.0);
  EXPECT_EQ(e.matrix(0,0), 1.0); // identity matrix
  EXPECT_EQ(e.matrix(0,1), 0.0);
  EXPECT_EQ(e.matrix(1,1), 1.0);
  EXPECT_EQ(e.matrix(2,2), 1.0);
  e.subtract_Egs(1.0);
  EXPECT_EQ(e.value_zero[0], 0.0);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
