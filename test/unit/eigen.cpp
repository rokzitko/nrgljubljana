#include <string>
#include <sstream>
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

template<typename T>
void VECTOR_EQ(const ublas::vector<T> &A, const ublas::vector<T> &B)
{
  EXPECT_EQ(A.size(), B.size());
  for(int i = 0; i < A.size(); i++) EXPECT_EQ(A[i], B[i]);
}

template<typename T>
void MATRIX_EQ(const ublas::matrix<T> &A, const ublas::matrix<T> &B)
{
  EXPECT_EQ(A.size1(), B.size1());
  EXPECT_EQ(A.size2(), B.size2());
  for(int i = 0; i < A.size1(); i++) 
    for(int j = 0; j < A.size2(); j++) 
      EXPECT_EQ(A(i,j), B(i,j));
}

TEST(io, read_vector) { // NOLINT
  std::string data = "5 1 2 3 4 5";
  std::istringstream ss(data);
  auto vec = read_vector<double>(ss);
  EVEC ref(5);
  ref[0] = 1; ref[1] = 2; ref[2] = 3; ref[3] = 4; ref[4] = 5;
  VECTOR_EQ(vec, ref);
}

TEST(io, read_matrix) { // NOLINT
  std::string data = "1 2 3\n 4 5 6\n";
  std::istringstream ss(data);
  auto mat = read_matrix<double>(ss, 2, 3);
  ublas::matrix<double> ref(2,3);
  ref(0,0) = 1; ref(0,1) = 2; ref(0,2) = 3; ref(1,0) = 4; ref(1,1) = 5; ref(1,2) = 6;
  MATRIX_EQ(mat, ref);
}

TEST(Diag, constructor) { // NOLINT
  std::string data = 
    "0 1\n"
    "2\n";
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
