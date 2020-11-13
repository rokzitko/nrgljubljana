// Discretization ODE solver for NRG
// ** Loading (and parsing) of tabulated data

#ifndef _adapt_load_hpp_
#define _adapt_load_hpp_

#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <sstream>
#include <stdexcept>
#include <cctype>
#include <cstdlib>

namespace NRG::Adapt {

enum class Sign { POS, NEG }; // positive vs. negative energies

// Split a string 's' into substrings. Leading spaces are ignored.
inline auto split_string(const std::string &s, unsigned int atleast = 0) {
  const int len = s.length();
  int index = 0;
  while (index < len && isspace(s[index])) { index++; }
  std::vector<std::string> substrings;
  while (index < len) {
    std::string substr = "";
    // Copy string until space or end of string
    while (index < len && !isspace(s[index])) {
      substr += s[index];
      index++;
    }
    substrings.push_back(substr);
    // Locate new substring
    while (index < len && isspace(s[index])) { index++; }
  }
  if (substrings.size() < atleast) 
    throw std::runtime_error("At least " + std::to_string(atleast) + " columns expected.");
  return substrings;
}

inline auto load_g(const std::string &filename) {
  std::ifstream F;
  safe_open(F, filename);
  Vec vecg;
  while (F) {
    const auto line = getnextline(F);
    if (!F) break;
    const auto columns = split_string(line, 2);
    vecg.emplace_back(std::make_pair(atof(columns[0]), atof(columns[1])));
  }
  if (vecg.size() == 0) 
    throw std::runtime_error("No data found.");
  return vecg;
}

inline void rescalevecxy(Vec &vec, const double factorx, const double factory) {
  for (int i = 0; i < vec.size(); i++) {
    vec[i].first  *= factorx;
    vec[i].second *= factory;
  }
  std::cout << "Rescaled to the interval [ " << vec.front().first << " : " << vec.back().first << " ]" << std::endl;
}

// Show minimal and maximal y in a table.
inline void minmaxvec(const Vec &vec, const std::string name) {
  auto miny   = DBL_MAX;
  auto maxy   = 0;
  for (int i = 0; i < vec.size(); i++) {
    const auto y = vec[i].second;
    if (y > maxy) { maxy = y; }
    if (y < miny) { miny = y; }
  }
  std::cout << "# min[" << name << "]=" << miny << " max[" << name << "]=" << maxy << std::endl;
}

// Load positive (sign=POS) or negative (sogn=NEG) part of the hybridisation function into a vector.
inline Vec load_rho(const std::string &filename, const Sign sign) {
  std::ifstream F;
  safe_open(F, filename);
  Vec vecrho;
  while (F) {
    const auto line = getnextline(F);
    if (!F) break;
    const auto columns = split_string(line, 2);
    const auto  x = atof(columns[0]);
    const auto  y = atof(columns[1]);
    if ((sign == Sign::POS && x > 0) || (sign == Sign::NEG && x < 0)) {
      // y must be positive (or zero)
      if (y < 0.0) 
        throw std::runtime_error("Negative y found.");
      // Disregard sign of x !
      vecrho.emplace_back(std::make_pair(abs(x), y));
    }
  }
  if (vecrho.size() == 0) 
    throw std::runtime_error("No data found.");
  std::sort(vecrho.begin(), vecrho.end());
  std::cout << "# " << filename << " ";
  std::cout << "- " << (sign == Sign::POS ? "POS" : "NEG") << " ";
  std::cout << "- interval [ " << vecrho.front().first << " : ";
  std::cout << vecrho.back().first << " ]" << std::endl;
  return vecrho;
}

inline std::string tostring(const Pair &p) {
  std::ostringstream str;
  str << p.first << " " << p.second;
  return str.str();
}

inline void save(const std::string &fn, const Vec &v) {
  std::ofstream F(fn.c_str());
  if (!F) 
    throw std::runtime_error("Failed to open " + fn + " for writing.");
  std::transform(v.begin(), v.end(), std::ostream_iterator<std::string>(F, "\n"), tostring);
}

inline void save(const std::string &fn, const std::vector<double> &v) {
  std::ofstream F(fn.c_str());
  if (!F) 
    throw std::runtime_error("Failed to open " + fn + " for writing.");
  F << std::setprecision(18);
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(F, "\n"));
}

inline void load(const std::string &fn, std::vector<double> &v) {
  std::ifstream F;
  safe_open(F, fn);
  v.clear();
  while (F) {
    const auto line = getnextline(F);
    if (!F) break;
    const auto x = atof(line.c_str());
    v.push_back(x);
  }
}

} // namespace

#endif
