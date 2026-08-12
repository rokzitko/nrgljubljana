#include <gtest/gtest.h>
#include <resample/resample.hpp>

#include <cstdio>
#include <fstream>

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

TEST(resample, empty_input_throws) {
  std::vector<std::pair<double, double>> input;
  std::vector<std::pair<double, double>> grid = {{1.0, 0.0}};
  EXPECT_THROW(NRG::Resample::Resample<double>(input, grid), std::runtime_error);
}

TEST(resample, invalid_akima_input_throws) {
  std::vector<std::pair<double, double>> input = {{0.0, 1.0}, {0.0, 2.0}, {1.0, 3.0}};
  std::vector<std::pair<double, double>> grid = {{0.5, 0.0}};
  EXPECT_THROW(NRG::Resample::Resample<double>(input, grid), std::runtime_error);
}

TEST(resample, command_line_extrapolation) {
  const auto grid_filename = "resample_extrapolation_grid.dat";
  const auto default_output_filename = "resample_default_extrapolation.dat";
  const auto custom_output_filename = "resample_custom_extrapolation.dat";
  {
    std::ofstream grid_file(grid_filename);
    grid_file << "0 0\n1 0\n1.5 0\n6 0\n7 0\n";
  }

  char default_arg0[] = "resample";
  char default_arg1[] = "-e";
  char default_arg2[] = "txt/resample_input.txt";
  char default_arg3[] = "resample_extrapolation_grid.dat";
  char default_arg4[] = "resample_default_extrapolation.dat";
  char *default_argv[] = {default_arg0, default_arg1, default_arg2, default_arg3, default_arg4};

  optind = 1;
  NRG::Resample::Resample<double> default_resample(5, default_argv);
  default_resample.run();
  const auto default_output = readtable<double, double>(default_output_filename);
  const std::vector<std::pair<double, double>> default_expected = {
    {0.0, 0.1}, {1.0, 0.1}, {1.5, 0.1395833333333333}, {6.0, 2.5}, {7.0, 2.5}
  };
  compare(default_output, default_expected);

  char custom_arg0[] = "resample";
  char custom_arg1[] = "-e";
  char custom_arg2[] = "-a";
  char custom_arg3[] = "-1.25";
  char custom_arg4[] = "-b";
  char custom_arg5[] = "9.5";
  char custom_arg6[] = "txt/resample_input.txt";
  char custom_arg7[] = "resample_extrapolation_grid.dat";
  char custom_arg8[] = "resample_custom_extrapolation.dat";
  char *custom_argv[] = {
    custom_arg0, custom_arg1, custom_arg2, custom_arg3, custom_arg4,
    custom_arg5, custom_arg6, custom_arg7, custom_arg8
  };

  optind = 1;
  NRG::Resample::Resample<double> custom_resample(9, custom_argv);
  custom_resample.run();
  const auto custom_output = readtable<double, double>(custom_output_filename);
  const std::vector<std::pair<double, double>> custom_expected = {
    {0.0, -1.25}, {1.0, 0.1}, {1.5, 0.1395833333333333}, {6.0, 2.5}, {7.0, 9.5}
  };
  compare(custom_output, custom_expected);

  std::remove(grid_filename);
  std::remove(default_output_filename);
  std::remove(custom_output_filename);
}

TEST(resample, extrapolation_values_require_enable_flag) {
  char arg0[] = "resample";
  char arg1[] = "-a";
  char arg2[] = "0";
  char arg3[] = "input.dat";
  char arg4[] = "grid.dat";
  char arg5[] = "output.dat";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5};

  optind = 1;
  EXPECT_THROW(NRG::Resample::Resample<double>(6, argv), std::invalid_argument);
}
