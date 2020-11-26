#include <gtest/gtest.h>
#include <resample/resample.hpp>

// template<typename T1, typename T2, typename S1, typename S2>
// void compare_vector_pairs(const std::vector<std::pair<T1,T2>> &a, const std::vector<std::pair<T1,T2>> &b) {
//   ASSERT_EQ(a.size(), b.size());
//   for(int i = 0; i < a.size(); i++){
//     EXPECT_EQ(a[i], b[i]);
//   }
// }

using namespace NRG;

TEST(resample, basic){

  NRG::Resample::Resample<double> resample("txt/resample_input.txt", "txt/resample_grid.txt");
  auto const output = resample.run();
  auto const output_compare = readtable<double,double>("txt/resample_output.txt");
  
}