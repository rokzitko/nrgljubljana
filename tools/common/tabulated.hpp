#ifndef _tools_common_tabulated_hpp_
#define _tools_common_tabulated_hpp_

#include <algorithm>
#include <cfloat>
#include <cctype>
#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace NRG::Tools {

inline auto split_fields(const std::string &s, const unsigned int atleast = 0) {
  const int len = s.length();
  int index = 0;
  while (index < len && std::isspace(static_cast<unsigned char>(s[index]))) { index++; }
  std::vector<std::string> substrings;
  while (index < len) {
    std::string substr;
    while (index < len && !std::isspace(static_cast<unsigned char>(s[index]))) {
      substr += s[index];
      index++;
    }
    substrings.push_back(substr);
    while (index < len && std::isspace(static_cast<unsigned char>(s[index]))) { index++; }
  }
  if (substrings.size() < atleast)
    throw std::runtime_error("At least " + std::to_string(atleast) + " columns expected.");
  return substrings;
}

template<typename Vec, typename Stream, typename NextLine>
auto load_pairs(Stream &input, NextLine next_line) {
  Vec vec;
  while (input) {
    const auto line = next_line(input);
    if (!input) break;
    const auto columns = split_fields(line, 2);
    vec.emplace_back(std::make_pair(std::atof(columns[0].c_str()), std::atof(columns[1].c_str())));
  }
  if (vec.empty())
    throw std::runtime_error("No data found.");
  return vec;
}

template<typename Vec>
void rescale_xy(Vec &vec, const double factorx, const double factory, std::ostream &out = std::cout) {
  for (auto &xy : vec) {
    xy.first *= factorx;
    xy.second *= factory;
  }
  out << "Rescaled to the interval [ " << vec.front().first << " : " << vec.back().first << " ]" << std::endl;
}

template<typename Vec>
void print_minmax(const Vec &vec, const std::string &name, std::ostream &out = std::cout) {
  auto miny = DBL_MAX;
  auto maxy = 0.0;
  for (const auto &xy : vec) {
    const auto y = xy.second;
    if (y > maxy) { maxy = y; }
    if (y < miny) { miny = y; }
  }
  out << "# min[" << name << "]=" << miny << " max[" << name << "]=" << maxy << std::endl;
}

template<typename Vec, typename Stream, typename Sign, typename AcceptSign, typename NextLine>
auto load_abs_pairs(Stream &input, const Sign sign, AcceptSign accept_sign, NextLine next_line) {
  Vec vec;
  while (input) {
    const auto line = next_line(input);
    if (!input) break;
    const auto columns = split_fields(line, 2);
    const auto x = std::atof(columns[0].c_str());
    const auto y = std::atof(columns[1].c_str());
    if (accept_sign(sign, x)) {
      if (y < 0.0)
        throw std::runtime_error("Negative y found.");
      vec.emplace_back(std::make_pair(std::abs(x), y));
    }
  }
  if (vec.empty())
    throw std::runtime_error("No data found.");
  std::sort(vec.begin(), vec.end());
  return vec;
}

template<typename Vec>
void print_interval(const std::string &filename, const std::string &label, const Vec &vec, std::ostream &out = std::cout) {
  out << "# " << filename << " - " << label << " - interval [ " << vec.front().first << " : " << vec.back().first << " ]" << std::endl;
}

template<typename Pair>
std::string pair_to_string(const Pair &p) {
  std::ostringstream str;
  str << p.first << " " << p.second;
  return str.str();
}

template<typename Vec>
void save_pairs(std::ostream &out, const Vec &vec) {
  std::transform(vec.begin(), vec.end(), std::ostream_iterator<std::string>(out, "\n"), pair_to_string<typename Vec::value_type>);
}

inline void save_values(std::ostream &out, const std::vector<double> &values) {
  out << std::setprecision(18);
  std::copy(values.begin(), values.end(), std::ostream_iterator<double>(out, "\n"));
}

template<typename Stream, typename NextLine>
void load_values(Stream &input, std::vector<double> &values, NextLine next_line) {
  values.clear();
  while (input) {
    const auto line = next_line(input);
    if (!input) break;
    values.push_back(std::atof(line.c_str()));
  }
}

} // namespace NRG::Tools

#endif
