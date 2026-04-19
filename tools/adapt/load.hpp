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

#include "../common/tabulated.hpp"

namespace NRG::Adapt {

enum class Sign { POS, NEG }; // positive vs. negative energies

// Split a string 's' into substrings. Leading spaces are ignored.
inline auto split_string(const std::string &s, unsigned int atleast = 0) {
  return NRG::Tools::split_fields(s, atleast);
}

inline auto load_g(const std::string &filename) {
  std::ifstream F;
  safe_open(F, filename);
  return NRG::Tools::load_pairs<Vec>(F, [](std::ifstream &input) { return getnextline(input); });
}

inline void rescalevecxy(Vec &vec, const double factorx, const double factory) {
  NRG::Tools::rescale_xy(vec, factorx, factory);
}

// Show minimal and maximal y in a table.
inline void minmaxvec(const Vec &vec, const std::string name) {
  NRG::Tools::print_minmax(vec, name);
}

// Load positive (sign=POS) or negative (sogn=NEG) part of the hybridisation function into a vector.
inline Vec load_rho(const std::string &filename, const Sign sign) {
  std::ifstream F;
  safe_open(F, filename);
  auto vecrho = NRG::Tools::load_abs_pairs<Vec>(
    F,
    sign,
    [](const Sign s, const double x) { return (s == Sign::POS && x > 0) || (s == Sign::NEG && x < 0); },
    [](std::ifstream &input) { return getnextline(input); });
  NRG::Tools::print_interval(filename, sign == Sign::POS ? "POS" : "NEG", vecrho);
  return vecrho;
}

inline void save(const std::string &fn, const Vec &v) {
  std::ofstream F(fn.c_str());
  if (!F) 
    throw std::runtime_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_pairs(F, v);
}

inline void save(const std::string &fn, const std::vector<double> &v) {
  std::ofstream F(fn.c_str());
  if (!F) 
    throw std::runtime_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_values(F, v);
}

inline void load(const std::string &fn, std::vector<double> &v) {
  std::ifstream F;
  safe_open(F, fn);
  NRG::Tools::load_values(F, v, [](std::ifstream &input) { return getnextline(input); });
}

} // namespace

#endif
