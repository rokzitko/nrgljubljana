#include <gtest/gtest.h>
#include <resample/resample.hpp>

template<typename T1, typename T2, typename S1, typename S2>
void compare(const std::vector<std::pair<T1,T2>> &a, const std::vector<std::pair<S1,S2>> &b) {
  ASSERT_EQ(a.size(), b.size());
  for(int i = 0; i < a.size(); i++){
    EXPECT_DOUBLE_EQ(a[i].first, b[i].first);
    EXPECT_DOUBLE_EQ(a[i].second, b[i].second);
  }
}

using namespace NRG;

TEST(resample, basic){

  NRG::Resample::Resample<double> resample("txt/resample_input.txt", "txt/resample_grid.txt");
  auto const output = resample.run();
  auto const output_compare = readtable<double,double>("txt/resample_output.txt");
  compare(*output, output_compare);

}