#include <gtest/gtest.h>

#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <common/tabulated.hpp>

namespace {

using Pairs = std::vector<std::pair<double, double>>;

std::string next_line(std::istringstream &input) {
  std::string line;
  std::getline(input, line);
  return line;
}

auto load_pairs(const std::string &text) {
  std::istringstream input{text};
  return NRG::Tools::load_pairs<Pairs>(input, next_line);
}

auto load_abs_pairs(const std::string &text) {
  std::istringstream input{text};
  return NRG::Tools::load_abs_pairs<Pairs>(
    input, true, [](const bool, const double x) { return x > 0.0; }, next_line);
}

auto load_values(const std::string &text) {
  std::istringstream input{text};
  std::vector<double> values;
  NRG::Tools::load_values(input, values, next_line);
  return values;
}

} // namespace

TEST(Tabulated, finite_values_and_surrounding_whitespace_are_accepted) { // NOLINT
  const auto pairs = load_pairs(" 0 +1.25 ignored\n-2e1 3e-2\n");
  ASSERT_EQ(pairs.size(), 2U);
  EXPECT_DOUBLE_EQ(pairs[0].first, 0.0);
  EXPECT_DOUBLE_EQ(pairs[0].second, 1.25);
  EXPECT_DOUBLE_EQ(pairs[1].first, -20.0);
  EXPECT_DOUBLE_EQ(pairs[1].second, 0.03);

  const auto abs_pairs = load_abs_pairs("-2 3\n2 4\n");
  ASSERT_EQ(abs_pairs.size(), 1U);
  EXPECT_DOUBLE_EQ(abs_pairs[0].first, 2.0);
  EXPECT_DOUBLE_EQ(abs_pairs[0].second, 4.0);

  const auto values = load_values(" 0 \t\n+1.25\n-2e1   \n");
  ASSERT_EQ(values.size(), 3U);
  EXPECT_DOUBLE_EQ(values[0], 0.0);
  EXPECT_DOUBLE_EQ(values[1], 1.25);
  EXPECT_DOUBLE_EQ(values[2], -20.0);
}

TEST(Tabulated, pair_loaders_reject_invalid_numeric_fields) { // NOLINT
  for (const std::string invalid : {"invalid", "1garbage", "nan", "inf", "-inf", "1e9999", "1e-9999"}) {
    SCOPED_TRACE(invalid);
    EXPECT_THROW(load_pairs(invalid + " 2\n"), std::runtime_error);
    EXPECT_THROW(load_pairs("1 " + invalid + "\n"), std::runtime_error);

    EXPECT_THROW(load_abs_pairs(invalid + " 2\n1 2\n"), std::runtime_error);
    EXPECT_THROW(load_abs_pairs("-1 " + invalid + "\n1 2\n"), std::runtime_error);
  }
}

TEST(Tabulated, single_value_loader_rejects_invalid_numeric_fields) { // NOLINT
  for (const std::string invalid : {"invalid", "1garbage", "nan", "inf", "-inf", "1e9999", "1e-9999", "1 trailing"}) {
    SCOPED_TRACE(invalid);
    EXPECT_THROW(load_values(invalid + "\n"), std::runtime_error);
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
