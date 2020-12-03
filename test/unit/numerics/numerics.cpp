#include <gtest/gtest.h>
#include <numerics.hpp>

using namespace NRG;

TEST(numerics, reim) {
  const auto [r, i] = reim(std::complex(1.0,2.0));
  EXPECT_EQ(r, 1.0);
  EXPECT_EQ(i, 2.0);
}

TEST(numerics, sum2){
	std::vector v = {std::pair{4,12}, std::pair{51,123}, std::pair{412,441}};
	EXPECT_EQ(sum2(v), 576);

	std::vector<std::pair<int,int>> empty_vec = {};
	EXPECT_EQ(sum2(empty_vec), 0);
}

TEST(numerics, conj_me){
  EXPECT_EQ(conj_me(std::complex(5.0,8.0)), std::complex(5.0,-8.0));
}
TEST(numerics, Zero_matrix){
  const size_t dim1 = 4;
  const size_t dim2 = 2;
  const size_t dim3 = 3;
  auto zero_m1 = Zero_matrix<double>(dim1,dim2);
  auto zero_m2 = Zero_matrix<double>(dim3);

  ASSERT_EQ(zero_m1.size1(), dim1);
  ASSERT_EQ(zero_m1.size2(), dim2);
  ASSERT_EQ(zero_m2.size1(), dim3);
  ASSERT_EQ(zero_m2.size2(), dim3);

  for(size_t i = 0; i < dim1; i++){
    for(size_t j = 0; j < dim2; j++){
      EXPECT_EQ(zero_m1(i,j), 0);
    }
  }

  for(size_t i = 0; i < dim3; i++){
    for(size_t j = 0; j < dim3; j++){
      EXPECT_EQ(zero_m2(i,j), 0);
    }
  }

}

TEST(numerics, bucket){
  std::vector v = {std::pair{4.0,12.0}, std::pair{51.0,123.0}, std::pair{412.0,441.0}};
  bucket sum(v);
  EXPECT_EQ(sum, 576);
  sum += 10;
  EXPECT_EQ(sum, 586);
}

TEST(numerics, is_even_is_odd){
  const int odd = 3;
  const int even = 4;
  EXPECT_TRUE(is_even(even));
  EXPECT_FALSE(is_odd(even));
  EXPECT_TRUE(is_odd(odd));
  EXPECT_FALSE(is_even(odd));
}

TEST(numerics, my_fcmp){
  const double a = 1.0;
  const double b = 1.0001;
  EXPECT_EQ(my_fcmp(a, b, 0.1, 0.00015), 0);
  EXPECT_EQ(my_fcmp(a, b, 0.1, 0.00001), -1);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
