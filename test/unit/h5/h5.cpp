#include <gtest/gtest.h>

#include <type_traits>
#include <complex>

#include <h5.hpp>

#include "compare.hpp"

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

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
