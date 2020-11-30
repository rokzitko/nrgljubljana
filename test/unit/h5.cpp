#include <gtest/gtest.h>

#include <type_traits>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define H5_USE_BOOST
#define H5_USE_EIGEN
#include <highfive/H5File.hpp>

#include <h5.hpp>

using namespace boost::numeric;

TEST(traits, is_same_v) { // NOLINT
  {
    auto res = std::is_same_v<int,  int>; 
    EXPECT_TRUE(res);
  }
  {
    auto res = std::is_same_v<void, void>;
    EXPECT_TRUE(res);
  }
  {
    auto res = std::is_same_v<int,  void>;
    EXPECT_FALSE(res);
  }
}

template <class T>
struct is_ublas_vector : std::false_type {};
template <class T>
struct is_ublas_vector<ublas::vector<T>> : std::true_type {};

template <class T>
struct is_ublas_matrix : std::false_type {};
template <class T>
struct is_ublas_matrix<ublas::matrix<T>> : std::true_type {};

template <class T>
struct is_ublas_complex_matrix : std::false_type {};
template <class T>
struct is_ublas_complex_matrix<ublas::matrix<std::complex<T>>> : std::true_type {};

// TO DO: is_ublas_real_matrix

TEST(traits, is_ublas_vector) { // NOLINT
  EXPECT_FALSE(is_ublas_vector<int>::value);
  EXPECT_TRUE(is_ublas_vector<ublas::vector<int>>::value);
  EXPECT_TRUE(is_ublas_vector<ublas::vector<double>>::value);
}

TEST(traits, is_ublas_matrix) { // NOLINT
  EXPECT_FALSE(is_ublas_matrix<int>::value);
  EXPECT_TRUE(is_ublas_matrix<ublas::matrix<int>>::value);
  EXPECT_TRUE(is_ublas_matrix<ublas::matrix<double>>::value);
}

TEST(traits, is_ublas_complex_matrix) { // NOLINT
  EXPECT_FALSE(is_ublas_complex_matrix<int>::value);
  EXPECT_FALSE(is_ublas_complex_matrix<ublas::matrix<int>>::value);
  EXPECT_FALSE(is_ublas_complex_matrix<ublas::matrix<double>>::value);
  EXPECT_TRUE (is_ublas_complex_matrix<ublas::matrix<std::complex<double>>>::value);
}

TEST(h5dump, long_path) { // NOLINT
  H5Easy::File file("long_path.h5", H5Easy::File::Overwrite);
  const int w = 42;
  const std::string path = "/this/is/a/long/path";
  H5Easy::dump(file, path, w);
  const auto r = H5Easy::load<int>(file, path);
  EXPECT_EQ(w, r);
}

TEST(h5dump, scalar_int) { // NOLINT
  H5Easy::File file("scalar_int.h5", H5Easy::File::Overwrite);
  int w = 42;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<int>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(h5dump, scalar_size_t) { // NOLINT
  H5Easy::File file("scalar_size_t.h5", H5Easy::File::Overwrite);
  size_t w = 42;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<size_t>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(h5dump, scalar_double) { // NOLINT
  H5Easy::File file("scalar_double.h5", H5Easy::File::Overwrite);
  double w = 42.0;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<double>(file, "/path");
  EXPECT_EQ(w, r);
}

#ifdef H5_COMPLEX_IMPLEMENTED
TEST(h5dump, scalar_complex_double) { // NOLINT
  H5Easy::File file("scalar_complex_double.h5", H5Easy::File::Overwrite);
  std::complex<double> w = 42.0;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<std::complex<double>>(file, "/path");
  EXPECT_EQ(w, r);
}
#endif

TEST(h5dump, std_vector_double) { // NOLINT
  H5Easy::File file("std_vector_double.h5", H5Easy::File::Overwrite);
  std::vector w = {1.0, 2.0, 3.0};
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<std::vector<double>>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(h5dump, ublas_vector_double) { // NOLINT
  H5Easy::File file("ublas_vector_double.h5", H5Easy::File::Overwrite);
  ublas::vector<double> w(3);
  w[0] = 1.0;
  w[1] = 2.0;
  w[2] = 3.0;
  H5Easy::dump(file, "/path", w);
  // read-back as std::vector
  auto r = H5Easy::load<std::vector<double>>(file, "/path");
  EXPECT_EQ(w[2], r[2]);
}

TEST(h5dump, ublas_matrix_double2) { // NOLINT
  HighFive::File file("ublas_matrix_double2.h5", HighFive::File::Overwrite);
  ublas::matrix<double> w(2,2);
  w(0,0) = 1.0;
  w(0,1) = 2.0;
  w(1,0) = 3.0;
  w(1,1) = 4.0;
  HighFive::DataSet dataset = file.createDataSet<double>("/path", HighFive::DataSpace::From(w));
  dataset.write(w);
  ublas::matrix<double> r;
  dataset.read(r);
  EXPECT_EQ(w(1,1), r(1,1));
}

TEST(h5dump, ublas_matrix_double2_rect) { // NOLINT
  HighFive::File file("ublas_matrix_double2_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  ublas::matrix<double> w(nx,ny);
  auto cnt = 1.0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  HighFive::DataSet dataset = file.createDataSet<double>("/path", HighFive::DataSpace::From(w));
  dataset.write(w);
  ublas::matrix<double> r;
  dataset.read(r);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      EXPECT_EQ(w(x,y), r(x,y));
    }
  }
}

#ifdef H5_COMPLEX_IMPLEMENTED
TEST(h5dump, ublas_matrix_complex2_rect) { // NOLINT
  HighFive::File file("ublas_matrix_complex2_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  ublas::matrix<std::complex<double>> w(nx,ny);
  auto cnt = std::complex<double>(1.0,1.0);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  HighFive::DataSet dataset = file.createDataSet<double>("/path", HighFive::DataSpace::From(w));
  dataset.write(w);
  ublas::matrix<std::complex<double>> r;
  dataset.read(r);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      EXPECT_EQ(w(x,y), r(x,y));
    }
  }
}
#endif

TEST(h5dump, ublas_matrix_double_rect) { // NOLINT
  H5Easy::File file("ublas_matrix_double_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  ublas::matrix<double> w(nx,ny);
  auto cnt = 1.0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  H5Easy::dump(file, "/path", w);

  // read-in as ublas::matrix but using the general HighFive interface
  auto dataset = file.getDataSet("/path");
  ublas::matrix<double> r;
  dataset.read(r);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      EXPECT_EQ(w(x,y), r(x,y));
    }
  }
}

TEST(h5dump, ublas_matrix_real_part_rect) { // NOLINT
  H5Easy::File file("ublas_matrix_real_part_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  ublas::matrix<std::complex<double>> w(nx,ny);
  auto cnt = std::complex<double>(1.0,1.0);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  NRG::h5_dump_matrix(file, "/path", w);

  // read-in as ublas::matrix but using the general HighFive interface
  auto dataset = file.getDataSet("/path");
  ublas::matrix<double> r;
  dataset.read(r);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      EXPECT_EQ(w(x,y).real(), r(x,y));
    }
  }
}

template <typename T1, typename T2, int N, int M, int K, int L>
void compare_matrices(Eigen::Matrix<T1,N,M> a, Eigen::Matrix<T2,K,L> b){
    ASSERT_EQ(a.rows(), b.rows());
    ASSERT_EQ(a.cols(), b.cols());
    for(int i = 0; i < a.rows(); i++)
        for(int j = 0; j < a.cols(); j++)
            EXPECT_EQ(a(i,j), b(i,j));
}

TEST(h5dump, _eigen_dump_matrix) { // NOLINT
  H5Easy::File file("_eigen_ublas_matrix.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  Eigen::MatrixXd w(nx,ny);
  double cnt = 1.0;
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  NRG::_eigen_h5_dump_matrix(file, "/path", w);

  auto dataset = file.getDataSet("/path");
  Eigen::MatrixXd r;
  dataset.read(r);
  compare_matrices(w,r);
}

TEST(h5dump, eigen_matrix_real_part_rect) { // NOLINT
  H5Easy::File file("ublas_matrix_real_part_rect.h5", HighFive::File::Overwrite);
  const auto nx = 2;
  const auto ny = 3;
  Eigen::MatrixXcd w(nx,ny);
  auto cnt = std::complex<double>(1.0,1.0);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      w(x,y) = cnt;
      cnt = cnt + 1.0;
    }
  }
  NRG::_eigen_h5_dump_matrix(file, "/path", w);

  auto dataset = file.getDataSet("/path");
  Eigen::MatrixXd r;
  dataset.read(r);
  for (int x = 0; x < nx; x++) {
    for (int y = 0; y < ny; y++) {
      EXPECT_EQ(w(x,y).real(), r(x,y));
    }
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
