// Discretization ODE solver for NRG
//
// ** Loading (and parsing) of tabulated data

#include <algorithm>

#include "../common/tabulated.hpp"

[[noreturn]] inline void fail_with_error(const std::string &message) {
  cerr << "ERROR: " << message << endl;
  exit(1);
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
vector<string> split_string(const string &s, unsigned int atleast = 0) {
  return run_or_exit([&] { return NRG::Tools::split_fields(s, atleast); });
}

Vec load_g(const string &filename) {
  ifstream F;
  safe_open(F, filename);
  return run_or_exit([&] { return NRG::Tools::load_pairs<Vec>(F, [](ifstream &input) { return getnextline(input); }); });
}

void rescalevecxy(Vec &vec, double factorx, double factory) {
  NRG::Tools::rescale_xy(vec, factorx, factory);
}

// Show minimal and maximal y in a table.
void minmaxvec(Vec &vec, string name) {
  NRG::Tools::print_minmax(vec, name);
}

enum SIGN { POS, NEG }; // positive vs. negative energies

// Load positive (sign=POS) or negative (sogn=NEG) part of the
// hybridisation function into a vector.
Vec load_rho(const string &filename, SIGN sign) {
  ifstream F;
  safe_open(F, filename);
  auto vecrho = run_or_exit([&] {
    return NRG::Tools::load_abs_pairs<Vec>(
      F,
      sign,
      [](const SIGN s, const double x) { return (s == POS && x > 0) || (s == NEG && x < 0); },
      [](ifstream &input) { return getnextline(input); });
  });
  NRG::Tools::print_interval(filename, sign == POS ? "POS" : "NEG", vecrho);
  return vecrho;
}

void save(const string &fn, const Vec &v) {
  ofstream F(fn.c_str());
  if (!F) fail_with_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_pairs(F, v);
}

void save(const string &fn, const vector<double> &v) {
  ofstream F(fn.c_str());
  if (!F) fail_with_error("Failed to open " + fn + " for writing.");
  NRG::Tools::save_values(F, v);
}

void load(const string &fn, vector<double> &v) {
  ifstream F;
  safe_open(F, fn);
  run_or_exit([&] {
    NRG::Tools::load_values(F, v, [](ifstream &input) { return getnextline(input); });
    return 0;
  });
}
