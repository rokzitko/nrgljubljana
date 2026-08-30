#include <gtest/gtest.h>

#include <array>
#include <iomanip>
#include <ios>
#include <limits>
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

auto read_strict_pairs(const std::string &text, const std::string &source = "table.dat") {
  std::istringstream input{text};
  return NRG::Tools::read_strict_pairs(input, source);
}

auto strict_read_error(std::istream &input, const std::string &source = "table.dat") {
  try {
    (void)NRG::Tools::read_strict_pairs(input, source);
  } catch (const std::runtime_error &error) {
    return std::string(error.what());
  }
  return std::string{};
}

auto strict_read_error(const std::string &text, const std::string &source = "table.dat") {
  std::istringstream input{text};
  return strict_read_error(input, source);
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

TEST(Tabulated, strict_pair_reader_accepts_comments_subnormals_and_empty_data) { // NOLINT
  std::ostringstream input;
  input << "\n\t \r\n# comment\n  \f\v# indented comment\n0 -0\n"
        << std::setprecision(std::numeric_limits<double>::max_digits10) << std::numeric_limits<double>::denorm_min()
        << " 2";
  const auto pairs = read_strict_pairs(input.str(), "points.dat");

  ASSERT_EQ(pairs.size(), 2U);
  EXPECT_DOUBLE_EQ(pairs[0].first, 0.0);
  EXPECT_DOUBLE_EQ(pairs[0].second, -0.0);
  EXPECT_DOUBLE_EQ(pairs[1].first, std::numeric_limits<double>::denorm_min());
  EXPECT_DOUBLE_EQ(pairs[1].second, 2.0);

  EXPECT_TRUE(read_strict_pairs("\n \t\r\n# comment\n  # another comment").empty());
}

TEST(Tabulated, strict_pair_reader_requires_exactly_two_fields) { // NOLINT
  EXPECT_EQ(strict_read_error("# heading\n\n1\n"),
            "table.dat:3: expected exactly 2 numeric fields; found 1.");
  EXPECT_EQ(strict_read_error("1 2 3\n", "extra.dat"),
            "extra.dat:1: expected exactly 2 numeric fields; found 3.");
  EXPECT_EQ(strict_read_error("1 2 # inline comment\n"),
            "table.dat:1: expected exactly 2 numeric fields; found 5.");
}

TEST(Tabulated, strict_pair_reader_rejects_invalid_numeric_fields_with_context) { // NOLINT
  for (const std::string invalid : {"invalid", "1garbage", "nan", "inf", "-inf", "1e9999", "1e-9999"}) {
    SCOPED_TRACE(invalid);
    const auto error = strict_read_error("# heading\n\n1 " + invalid + "\n", "values.dat");
    EXPECT_NE(error.find("values.dat:3: field 2: expected a finite representable number"), std::string::npos);
    EXPECT_NE(error.find("got '" + invalid + "'"), std::string::npos);
  }

  EXPECT_EQ(strict_read_error("\ninvalid 2", "first-field.dat"),
            "first-field.dat:2: field 1: expected a finite representable number; got 'invalid'.");
}

TEST(Tabulated, strict_pair_reader_handles_unterminated_lines_and_exception_masks) { // NOLINT
  const std::array masks{std::ios::goodbit, std::ios::eofbit, std::ios::failbit, std::ios::badbit,
                         std::ios::eofbit | std::ios::failbit | std::ios::badbit};
  for (const auto mask : masks) {
    for (const std::string text : {"# heading\n1 2", "# heading\n1 2\n"}) {
      SCOPED_TRACE(static_cast<int>(mask));
      SCOPED_TRACE(text);
      std::istringstream input{text};
      input.exceptions(mask);
      EXPECT_NO_THROW({
        const auto pairs = NRG::Tools::read_strict_pairs(input, "masked.dat");
        ASSERT_EQ(pairs.size(), 1U);
        EXPECT_EQ(pairs[0], (std::pair<double, double>{1.0, 2.0}));
      });
    }
  }
}

TEST(Tabulated, strict_pair_reader_distinguishes_io_failures) { // NOLINT
  std::istringstream broken;
  broken.setstate(std::ios::badbit);
  EXPECT_EQ(strict_read_error(broken, "broken.dat"), "broken.dat: I/O error after line 0.");

  std::istringstream exception_enabled;
  exception_enabled.exceptions(std::ios::badbit);
  EXPECT_THROW(exception_enabled.setstate(std::ios::badbit), std::ios_base::failure);
  const auto error = strict_read_error(exception_enabled, "masked-broken.dat");
  EXPECT_EQ(error.find("masked-broken.dat: I/O error after line 0"), 0U);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS(); // NOLINT
}
