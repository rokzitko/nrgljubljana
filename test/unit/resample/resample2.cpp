#include <gtest/gtest.h>
#include <resample/resample.hpp>

#include <array>
#include <cstddef>
#include <cstdio>
#include <fstream>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

template<typename T1, typename T2, typename S1, typename S2>
void compare(const std::vector<std::pair<T1,T2>> &a, const std::vector<std::pair<S1,S2>> &b) {
  ASSERT_EQ(a.size(), b.size());
  for(std::size_t i = 0; i < a.size(); i++){
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

TEST(resample, interpolation_methods_produce_expected_sparse_nonlinear_values) {
  using NRG::Tools::InterpolationMethod;
  const std::vector<std::pair<double, double>> input = {
    {-3.0, 4.0}, {-1.0, -2.0}, {0.5, 1.0}, {2.0, 8.0}, {5.0, -1.0}, {9.0, 6.0}
  };
  const std::vector<std::pair<double, double>> grid = {
    {-2.0, 0.0}, {-0.25, 0.0}, {1.25, 0.0}, {3.5, 0.0}, {7.0, 0.0}
  };
  struct ExpectedInterpolation {
    InterpolationMethod method;
    const char *name;
    std::array<double, 5> values;
  };
  const std::array expectations{
    ExpectedInterpolation{InterpolationMethod::linear, "linear", {1.0, -0.5, 4.5, 3.5, 2.5}},
    ExpectedInterpolation{InterpolationMethod::cspline,
                          "cspline",
                          {0.11598527171152861, -1.4892196229050279, 4.9731422993905534, 5.3880935754189947,
                           -1.006348400203148}},
    ExpectedInterpolation{InterpolationMethod::akima,
                          "akima",
                          {-0.44021739130434767, -1.0234553775743709, 4.7142223536369015, 4.241459920066359,
                           0.40394295302013417}},
    ExpectedInterpolation{InterpolationMethod::steffen, "steffen", {0.25, -1.125, 5.125, 3.5, 1.625}}
  };

  for (const auto &expectation : expectations) {
    SCOPED_TRACE(expectation.name);
    NRG::Resample::Resample<double> resample(input, grid, std::nullopt, false, 16, expectation.method);
    const auto output = resample.run();
    ASSERT_TRUE(output.has_value());
    ASSERT_EQ(output->size(), grid.size());
    for (std::size_t index = 0; index < output->size(); ++index) {
      EXPECT_DOUBLE_EQ((*output)[index].first, grid[index].first);
      EXPECT_NEAR((*output)[index].second, expectation.values[index], 1e-12);
    }
  }
}

TEST(resample, interpolation_method_minimum_sizes) {
  using NRG::Tools::InterpolationMethod;
  struct MethodMinimum {
    InterpolationMethod method;
    const char *name;
    std::size_t minimum;
  };
  const std::array cases{
    MethodMinimum{InterpolationMethod::linear, "linear", 2},
    MethodMinimum{InterpolationMethod::cspline, "cspline", 3},
    MethodMinimum{InterpolationMethod::akima, "akima", 5},
    MethodMinimum{InterpolationMethod::steffen, "steffen", 3}
  };
  const std::vector<std::pair<double, double>> grid = {{0.5, 0.0}};

  for (const auto &test_case : cases) {
    SCOPED_TRACE(test_case.name);
    EXPECT_EQ(NRG::Tools::interpolation_minimum_size(test_case.method), test_case.minimum);
    std::vector<std::pair<double, double>> input;
    for (std::size_t index = 0; index < test_case.minimum; ++index) {
      const auto x = static_cast<double>(index);
      input.emplace_back(x, x * x + static_cast<double>(index % 2));
    }
    EXPECT_NO_THROW(NRG::Resample::Resample<double>(input, grid, std::nullopt, false, 16, test_case.method));

    input.pop_back();
    try {
      NRG::Resample::Resample<double> too_small(input, grid, std::nullopt, false, 16, test_case.method);
      FAIL() << test_case.name << " accepted too few input points";
    } catch (const std::runtime_error &error) {
      const auto expected = std::string(test_case.name) + " requires at least " + std::to_string(test_case.minimum);
      EXPECT_NE(std::string(error.what()).find(expected), std::string::npos);
    }
  }
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

TEST(resample, invalid_interpolation_method_cli_throws) {
  char arg0[] = "resample";
  char arg1[] = "--interpolation";
  char arg2[] = "polynomial";
  char arg3[] = "input.dat";
  char arg4[] = "grid.dat";
  char arg5[] = "output.dat";
  char *argv[] = {arg0, arg1, arg2, arg3, arg4, arg5};

  optind = 1;
  EXPECT_THROW(NRG::Resample::Resample<double>(6, argv), std::invalid_argument);
}

TEST(resample, command_line_selects_valid_interpolation_methods) {
  const std::string input_filename = "resample_cli_methods_input.dat";
  const std::string grid_filename = "resample_cli_methods_grid.dat";
  const std::string linear_filename = "resample_cli_methods_linear.dat";
  const std::string cspline_filename = "resample_cli_methods_cspline.dat";
  const std::string steffen_filename = "resample_cli_methods_steffen.dat";
  {
    std::ofstream input(input_filename);
    input << "-3 4\n-1 -2\n0.5 1\n2 8\n5 -1\n9 6\n";
    std::ofstream grid(grid_filename);
    grid << "-2 0\n";
  }

  std::vector<std::string> linear_arguments{"resample", "-i", "linear", input_filename, grid_filename, linear_filename};
  std::vector<char *> linear_argv;
  for (auto &argument : linear_arguments) linear_argv.push_back(argument.data());
  optind = 1;
  NRG::Resample::Resample<double> linear(static_cast<int>(linear_argv.size()), linear_argv.data());
  linear.run();

  std::vector<std::string> cspline_arguments{"resample", "--interpolation=cspline", input_filename, grid_filename, cspline_filename};
  std::vector<char *> cspline_argv;
  for (auto &argument : cspline_arguments) cspline_argv.push_back(argument.data());
  optind = 1;
  NRG::Resample::Resample<double> cspline(static_cast<int>(cspline_argv.size()), cspline_argv.data());
  cspline.run();

  std::vector<std::string> steffen_arguments{"resample", "-i", "steffen", input_filename, grid_filename, steffen_filename};
  std::vector<char *> steffen_argv;
  for (auto &argument : steffen_arguments) steffen_argv.push_back(argument.data());
  optind = 1;
  NRG::Resample::Resample<double> steffen(static_cast<int>(steffen_argv.size()), steffen_argv.data());
  steffen.run();

  const auto linear_output = readtable<double, double>(linear_filename);
  const auto cspline_output = readtable<double, double>(cspline_filename);
  const auto steffen_output = readtable<double, double>(steffen_filename);
  ASSERT_EQ(linear_output.size(), 1);
  ASSERT_EQ(cspline_output.size(), 1);
  ASSERT_EQ(steffen_output.size(), 1);
  EXPECT_DOUBLE_EQ(linear_output[0].second, 1.0);
  EXPECT_NEAR(cspline_output[0].second, 0.11598527171152861, 1e-12);
  EXPECT_NEAR(steffen_output[0].second, 0.25, 1e-12);

  std::remove(input_filename.c_str());
  std::remove(grid_filename.c_str());
  std::remove(linear_filename.c_str());
  std::remove(cspline_filename.c_str());
  std::remove(steffen_filename.c_str());
}
