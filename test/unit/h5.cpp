#include <gtest/gtest.h>

#include <type_traits>
#include <complex>

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#define H5_USE_BOOST
#include <highfive/H5File.hpp>

#include <h5.hpp>

using namespace boost::numeric;

TEST(is_same_v, traits) { // NOLINT
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

TEST(is_ublas_vector, traits) { // NOLINT
  EXPECT_FALSE(is_ublas_vector<int>::value);
  EXPECT_TRUE(is_ublas_vector<ublas::vector<int>>::value);
  EXPECT_TRUE(is_ublas_vector<ublas::vector<double>>::value);
}

TEST(is_ublas_matrix, traits) { // NOLINT
  EXPECT_FALSE(is_ublas_matrix<int>::value);
  EXPECT_TRUE(is_ublas_matrix<ublas::matrix<int>>::value);
  EXPECT_TRUE(is_ublas_matrix<ublas::matrix<double>>::value);
}

TEST(is_ublas_complex_matrix, traits) { // NOLINT
  EXPECT_FALSE(is_ublas_complex_matrix<int>::value);
  EXPECT_FALSE(is_ublas_complex_matrix<ublas::matrix<int>>::value);
  EXPECT_FALSE(is_ublas_complex_matrix<ublas::matrix<double>>::value);
  EXPECT_TRUE (is_ublas_complex_matrix<ublas::matrix<std::complex<double>>>::value);
}

TEST(long_path, h5dump) { // NOLINT
  H5Easy::File file("long_path.h5", H5Easy::File::Overwrite);
  const int w = 42;
  const std::string path = "/this/is/a/long/path";
  H5Easy::dump(file, path, w);
  const auto r = H5Easy::load<int>(file, path);
  EXPECT_EQ(w, r);
}

TEST(scalar_int, h5dump) { // NOLINT
  H5Easy::File file("scalar_int.h5", H5Easy::File::Overwrite);
  int w = 42;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<int>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(scalar_size_t, h5dump) { // NOLINT
  H5Easy::File file("scalar_size_t.h5", H5Easy::File::Overwrite);
  size_t w = 42;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<size_t>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(scalar_double, h5dump) { // NOLINT
  H5Easy::File file("scalar_double.h5", H5Easy::File::Overwrite);
  double w = 42.0;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<double>(file, "/path");
  EXPECT_EQ(w, r);
}

#ifdef H5_COMPLEX_IMPLEMENTED
TEST(scalar_complex_double, h5dump) { // NOLINT
  H5Easy::File file("scalar_complex_double.h5", H5Easy::File::Overwrite);
  std::complex<double> w = 42.0;
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<std::complex<double>>(file, "/path");
  EXPECT_EQ(w, r);
}
#endif

TEST(std_vector_double, h5dump) { // NOLINT
  H5Easy::File file("std_vector_double.h5", H5Easy::File::Overwrite);
  std::vector w = {1.0, 2.0, 3.0};
  H5Easy::dump(file, "/path", w);
  auto r = H5Easy::load<std::vector<double>>(file, "/path");
  EXPECT_EQ(w, r);
}

TEST(ublas_vector_double, h5dump) { // NOLINT
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

TEST(ublas_matrix_double2, h5dump) { // NOLINT
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

TEST(ublas_matrix_double2_rect, h5dump) { // NOLINT
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
TEST(ublas_matrix_complex2_rect, h5dump) { // NOLINT
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

TEST(ublas_matrix_double_rect, h5dump) { // NOLINT
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

TEST(ublas_matrix_real_part_rect, h5dump) { // NOLINT
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

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
