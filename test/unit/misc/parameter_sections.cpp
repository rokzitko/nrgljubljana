#include <gtest/gtest.h>

#include <cstdio>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#include <misc.hpp>
#include <parameter_sections.hpp>
#include <common/parser.hpp>

using namespace NRG;
using namespace std::string_literals;

namespace {

// The error message of read_parameter_sections(text), or "" if it is accepted
std::string error_of(const std::string &text) {
  try {
    read_parameter_sections(text);
  } catch (const std::runtime_error &e) { return e.what(); }
  return {};
}

bool contains(const std::string &s, const std::string &part) { return s.find(part) != std::string::npos; }

// Each section of accepted text is what the solver and the tools read for that name
void expect_read_alike(const std::string &text) {
  const auto sections = read_parameter_sections(text);
  const auto filename = "parameter_sections.txt"s;
  {
    std::ofstream file(filename);
    file << text;
  }
  for (const auto &[name, section] : sections) {
    std::istringstream solver_input(text);
    ASSERT_TRUE(find_block(solver_input, name)) << name;
    EXPECT_EQ(parse_block(solver_input), section) << name;

    std::ifstream tools_input(filename);
    ASSERT_TRUE(Tools::find_block(tools_input, name)) << name;
    std::map<std::string, std::string> tools_section;
    Tools::parse_key_value_block(tools_input, tools_section);
    EXPECT_EQ(tools_section, section) << name;
  }
  std::remove(filename.c_str());
}

} // namespace

// ---- Line syntax ----

TEST(parameter_sections, line) {
  using Kind = ParameterLine::Kind;
  EXPECT_EQ(parameter_line("").kind, Kind::skip);
  EXPECT_EQ(parameter_line(" \t\r").kind, Kind::skip);
  EXPECT_EQ(parameter_line("  # x=1").kind, Kind::skip);
  EXPECT_EQ(parameter_line(" [a b] ").kind, Kind::header);
  EXPECT_EQ(parameter_line(" [ a ] ").name, " a ");
  EXPECT_EQ(parameter_line("[a").kind, Kind::malformed_header);
  EXPECT_EQ(parameter_line("[").kind, Kind::malformed_header);
  EXPECT_EQ(parameter_line("[a] x").kind, Kind::malformed_header);
  EXPECT_EQ(parameter_line("no equals").kind, Kind::other);
  const auto entry = parameter_line(" G imp[0,0] = d:A=B \r");
  EXPECT_EQ(entry.kind, Kind::entry);
  EXPECT_EQ(entry.name, "G imp[0,0]");
  EXPECT_EQ(entry.value, "d:A=B");
}

// ---- Accepted text ----

TEST(parameter_sections, sections_and_keys) {
  const auto text = "[a]\nx=1\ny=two words\n[b]\nx=3\n"s;
  const auto s = read_parameter_sections(text);
  ASSERT_EQ(s.size(), 2u);
  EXPECT_EQ(s.at("a").at("x"), "1");
  EXPECT_EQ(s.at("a").at("y"), "two words");
  EXPECT_EQ(s.at("b").at("x"), "3"); // the same key in another section is a different key
  expect_read_alike(text);
}

TEST(parameter_sections, whitespace_comments_and_blank_lines) {
  const auto text = "# leading comment\n\n[sec]\n\t key  =  value \t\n   # indented comment\n\nk2=v2\r\n[nekaj drugega]\n"s;
  const auto s = read_parameter_sections(text);
  ASSERT_EQ(s.size(), 2u);
  EXPECT_EQ(s.at("sec").size(), 2u);
  EXPECT_EQ(s.at("sec").at("key"), "value");
  EXPECT_EQ(s.at("sec").at("k2"), "v2"); // CRLF on a key line
  EXPECT_TRUE(s.at("nekaj drugega").empty());
  expect_read_alike(text);
}

TEST(parameter_sections, split_at_first_equals) {
  const auto text = "[c]\nG imp[0,0]=d:A_d-A_d\nexpr=a=b\n"s;
  const auto s = read_parameter_sections(text);
  EXPECT_EQ(s.at("c").at("G imp[0,0]"), "d:A_d-A_d"); // keys may contain spaces and brackets
  EXPECT_EQ(s.at("c").at("expr"), "a=b");
  expect_read_alike(text);
}

TEST(parameter_sections, empty_values_and_sections) {
  const auto text = "[a]\nempty=\nspaces=   \n[b]\n"s;
  const auto s = read_parameter_sections(text);
  EXPECT_EQ(s.at("a").at("empty"), "");
  EXPECT_EQ(s.at("a").at("spaces"), "");
  EXPECT_TRUE(s.at("b").empty());
  EXPECT_TRUE(read_parameter_sections("").empty());
  EXPECT_TRUE(read_parameter_sections("# only a comment\n\n").empty());
  expect_read_alike(text);
}

TEST(parameter_sections, hash_after_value_is_part_of_value) {
  const auto text = "[a]\nx=1 # not a comment\n"s;
  EXPECT_EQ(read_parameter_sections(text).at("a").at("x"), "1 # not a comment");
  expect_read_alike(text);
}

TEST(parameter_sections, no_final_newline) {
  const auto text = "[a]\nx=1"s;
  EXPECT_EQ(read_parameter_sections(text).at("a").at("x"), "1");
  expect_read_alike(text);
}

TEST(parameter_sections, param_file) {
  const auto text = "[param]\nsymtype=QS\nLambda=2.0\nops=n_d A_d\n\n[extra]\n# model\nU1=0.5\neps1=-0.25\n"s;
  const auto s = read_parameter_sections(text);
  EXPECT_EQ(s.at("param").at("ops"), "n_d A_d");
  EXPECT_EQ(s.at("extra").at("eps1"), "-0.25");
  expect_read_alike(text);
}

// ---- Rejected text: each error names the line ----

TEST(parameter_sections, line_without_equals) {
  const auto e = error_of("[a]\nx=1\ntypo line\n");
  EXPECT_TRUE(contains(e, "line 3: expected key=value")) << e;
  EXPECT_TRUE(contains(e, "\"typo line\"")) << e;
}

TEST(parameter_sections, malformed_section_header) {
  EXPECT_TRUE(contains(error_of("[a\nx=1\n"), "line 1: malformed section header"));
  EXPECT_TRUE(contains(error_of("[a] trailing\n"), "malformed section header"));
  EXPECT_TRUE(contains(error_of("[\n"), "malformed section header"));
  EXPECT_TRUE(contains(error_of("[]\n"), "empty section name"));
}

TEST(parameter_sections, whitespace_in_section_header) {
  for (const auto *text : {"[ sec ]\n", "[sec ]\n", "  [sec]\n", "[sec]  \n", "[sec]\r\nx=1\n", "\t[sec]\n", "[  ]\n"})
    EXPECT_TRUE(contains(error_of(text), "line 1: whitespace in a section header")) << text;
}

TEST(parameter_sections, duplicate_section) {
  const auto e = error_of("[a]\nx=1\n[b]\n[a]\ny=2\n");
  EXPECT_TRUE(contains(e, "line 4: section [a] given twice")) << e;
}

TEST(parameter_sections, key_before_section) {
  const auto e = error_of("# comment\nx=1\n[a]\n");
  EXPECT_TRUE(contains(e, "line 2: key before the first section")) << e;
}

TEST(parameter_sections, empty_key) {
  EXPECT_TRUE(contains(error_of("[a]\n=value\n"), "line 2: empty key"));
  EXPECT_TRUE(contains(error_of("[a]\n  = value\n"), "empty key"));
}

TEST(parameter_sections, duplicate_key) {
  const auto e = error_of("[a]\nx=1\n x = 2\n");
  EXPECT_TRUE(contains(e, "line 3: key \"x\" given twice")) << e;
}

// ---- The solver reads one section and skips what the strict reader rejects ----

TEST(parameter_sections, solver_reading) {
  std::istringstream F("x=0\n[ param ]\na=1\n[param]\nnot a key\nb=2\n[a\nc=3\n[param]\nd=4\n[other]\ne=1\ne=2\n");
  ASSERT_TRUE(find_block(F, "param")); // [ param ] is another section
  EXPECT_EQ(parse_block(F), (std::map<std::string, std::string>{{"b", "2"}})); // ends at [a; the second [param] unread
  EXPECT_TRUE(find_block(F, " param "));
  EXPECT_FALSE(find_block(F, "a"));
}
