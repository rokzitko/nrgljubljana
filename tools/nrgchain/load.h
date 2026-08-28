// Discretization ODE solver for NRG
//
// ** Loading (and parsing) of tabulated data

#include <algorithm>
#include <exception>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../common/tabulated.hpp"

[[noreturn]] inline void fail_with_error(const std::string &message) {
  throw std::runtime_error(message);
}

template<typename F>
auto run_or_exit(F &&fn) -> decltype(fn()) {
  try {
    return fn();
  } catch (const std::exception &e) {
    fail_with_error(e.what());
  }
}

// Split a string 's' into substrings. Leading spaces are ignored.
std::vector<std::string> split_string(const std::string &s, unsigned int atleast = 0) {
  return run_or_exit([&] { return NRG::Tools::split_fields(s, atleast); });
}

Vec load_g(const std::string &filename) {
  std::ifstream F;
  safe_open(F, filename);
  return run_or_exit([&] { return NRG::Tools::load_pairs<Vec>(F, [](std::ifstream &input) { return getnextline(input); }); });
}

void rescalevecxy(Vec &vec, double factorx, double factory) {
  NRG::Tools::rescale_xy(vec, factorx, factory);
}

// Show minimal and maximal y in a table.
void minmaxvec(Vec &vec, std::string name) {
  NRG::Tools::print_minmax(vec, name);
}

enum SIGN { POS, NEG }; // positive vs. negative energies

// Load positive (sign=POS) or negative (sogn=NEG) part of the
// hybridisation function into a vector.
Vec load_rho(const std::string &filename, SIGN sign) {
  std::ifstream F;
  safe_open(F, filename);
  auto vecrho = run_or_exit([&] {
    return NRG::Tools::load_abs_pairs<Vec>(
      F,
      sign,
      [](const SIGN s, const double x) { return (s == POS && x > 0) || (s == NEG && x < 0); },
      [](std::ifstream &input) { return getnextline(input); });
  });
  NRG::Tools::print_interval(filename, sign == POS ? "POS" : "NEG", vecrho);
  return vecrho;
}

void save(const std::string &fn, const Vec &v) {
  std::ofstream F(fn.c_str());
  if (!F) fail_with_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_pairs(F, v);
  F.close();
  if (!F) fail_with_error("Failed writing " + fn + ".");
}

void save(const std::string &fn, const std::vector<double> &v) {
  std::ofstream F(fn.c_str());
  if (!F) fail_with_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_values(F, v);
  F.close();
  if (!F) fail_with_error("Failed writing " + fn + ".");
}

void load(const std::string &fn, std::vector<double> &v) {
  std::ifstream F;
  safe_open(F, fn);
  run_or_exit([&] {
    NRG::Tools::load_values(F, v, [](std::ifstream &input) { return getnextline(input); });
    return 0;
  });
}
