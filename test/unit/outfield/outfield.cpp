#include <gtest/gtest.h>
#include <string>
using namespace std::string_literals;
#include <fstream>
#include <exception>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include <outfield.hpp>

constexpr int prec = 10;
constexpr int width = 20;

TEST(outfield, Outfield) {
  const auto desc = "Sz"s;
  auto o = NRG::Outfield(desc, prec, width);
  EXPECT_EQ(o.get_desc(), desc);
  const auto ref_h = fmt::format("{:>20} ", desc);
  std::ostringstream s_header;
  o.put_header(s_header);
  EXPECT_EQ(ref_h, s_header.str());

  const auto ref_v_nul = fmt::format("{:>20} ", ""s);
  std::ostringstream s_value_nul;
  o.put_value(s_value_nul);
  EXPECT_EQ(ref_v_nul, s_value_nul.str());

  constexpr auto val = 12.0;
  o.set_value(val);
  const auto ref_v = fmt::format("{:>20.10} ", val);
  std::ostringstream s_value;
  o.put_value(s_value);
  EXPECT_EQ(ref_v, s_value.str());
}

TEST(outfield, many) {
  const std::vector fields = { "A"s, "B"s, "C"s };
  NRG::Allfields all(fields, prec, width);

  std::ostringstream s_header;
  all.save_header(s_header);
  const auto ref_h = fmt::format("#{:>20} {:>20} {:>20} \n", "A"s, "B"s, "C"s);
  EXPECT_EQ(ref_h, s_header.str());

  all.set("A"s, 1.0);
  all.set("B"s, 2.0);
  all.set("C"s, 3.0);
  const auto ref_v = fmt::format(" {:>20.10} {:>20.10} {:>20.10} \n", 1.0, 2.0, 3.0);
  std::ostringstream s_values;
  all.save_values(s_values);
  EXPECT_EQ(ref_v, s_values.str());

  EXPECT_THROW(all.set("D"s, 4.0), std::runtime_error);
}

TEST(outfield, add1) {
  const std::vector fields = { "A"s, "C"s };
  NRG::Allfields all(fields, prec, width);
  all.add("B", 1);

  std::ostringstream s_header;
  all.save_header(s_header);
  const auto ref_h = fmt::format("#{:>20} {:>20} {:>20} \n", "A"s, "B"s, "C"s);
  EXPECT_EQ(ref_h, s_header.str());
}

TEST(outfield, add2) {
  const std::vector fields = { "A"s, "D"s };
  NRG::Allfields all(fields, prec, width);
  all.add("B", 1);
  all.add("C", 2);

  std::ostringstream s_header;
  all.save_header(s_header);
  const auto ref_h = fmt::format("#{:>20} {:>20} {:>20} {:>20} \n", "A"s, "B"s, "C"s, "D"s);
  EXPECT_EQ(ref_h, s_header.str());
}

TEST(outfield, merge1) {
  const std::vector fields1 = { "A"s, "D"s };
  NRG::Allfields all(fields1, prec, width);
  const std::vector fields2 = { "B"s, "C"s };
  NRG::Allfields more(fields2, prec, width);
  all.add(more);
  
  std::ostringstream s_header;
  all.save_header(s_header);
  const auto ref_h = fmt::format("#{:>20} {:>20} {:>20} {:>20} \n", "A"s, "D"s, "B"s, "C"s);
  EXPECT_EQ(ref_h, s_header.str());
}

TEST(outfield, merge2) {
  const std::vector fields1 = { "A"s, "D"s };
  NRG::Allfields all(fields1, prec, width);
  const std::vector fields2 = { "B"s, "C"s };
  NRG::Allfields more(fields2, prec, width);
  all.add(more, 1);
  
  std::ostringstream s_header;
  all.save_header(s_header);
  const auto ref_h = fmt::format("#{:>20} {:>20} {:>20} {:>20} \n", "A"s, "B"s, "C"s, "D"s);
  EXPECT_EQ(ref_h, s_header.str());
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
