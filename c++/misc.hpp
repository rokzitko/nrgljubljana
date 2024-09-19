// misc.h - Miscelaneous functions
// Copyright (C) 2005-2022 Rok Zitko

#ifndef _misc_hpp_
#define _misc_hpp_

#include <stdexcept>
#include <string>
#include <iterator>
#include <optional>
#include <deque>
#include <vector>
#include <set>
#include <map>
#include <fstream>
#include <iostream>
#include <sstream>
#include <cstring> // stdcasecmp
#include <exception>
#include <cstdio> // stdout
#include <unistd.h> // isatyy

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>

#include <fmt/format.h>

#include "basicio.hpp"
#include "portabil.hpp"

#include <Eigen/Dense>

namespace NRG {

inline bool is_stdout_redirected() {
  return !isatty(fileno(stdout));
}

inline auto file_exists(const char *fileName) {
  std::ifstream file(fileName);
  return file.good();
}

inline auto contains(const std::string &str, const char c) {
  return str.find(c) != std::string::npos;
}

// x raised to the power of n
template<typename T, typename N>
constexpr inline auto intpow(const T x, const N n) {
  assert(n >= 0);
  T res = 1;
  for (N i = 1; i <= n; i++) res *= x;
  return res;
}

template<typename T> auto get_back(T &d) { // usually T is list or deque
  if (d.empty()) throw std::runtime_error("Error: List empty! File: "s + (std::string)__FILE__ + " Line: "s +  std::to_string(__LINE__));
  auto i = d.back();
  d.pop_back();
  return i;
}

template<typename T> auto get_front(T &d) {
  if (d.empty()) throw std::runtime_error("Error: List empty! File: "s + (std::string)__FILE__ + " Line: "s +  std::to_string(__LINE__));
  auto i = d.front();
  d.pop_front();
  return i;
}

// switch statement with three cases
template <typename T, typename T1>
inline T switch3(const T1 x0, const T1 x1, const T y1, const T1 x2, const T y2, const T1 x3, const T y3) {
  if (x0 == x1) return y1;
  if (x0 == x2) return y2;
  if (x0 == x3) return y3;
  throw std::runtime_error("Error: No match! File: "s + (std::string)__FILE__ + " Line: "s +  std::to_string(__LINE__));
}

// Get next line from stream F, skipping empty lines and comments.
inline std::optional<std::string> nextline(std::istream &F) {
  std::string line;
  while (F) {
    std::getline(F, line);
    if (!F) return std::nullopt;       // bail out
    if (line.length() == 0) continue;  // skip empty lines
    if (line[0] == '#') continue;      // skip comment lines
    return line;
  }
  return std::nullopt;                 // error
}

inline std::string strip_trailing_whitespace(const std::string &in) {
  auto s(in);
  auto it = s.rbegin();
  while (it != s.rend() && std::isspace(*it)) {
    s.erase(--it.base());
    it = s.rbegin();
  }
  return s;
}

// Parse a block of "keyword=value" lines.
inline auto parse_block(std::istream &F) {
  std::map<std::string, std::string> parsed_params;
  while (F) {
    if (const auto l = nextline(F)) {
      const auto line = l.value();
      if (line[0] == '[') // new block, we're done!
        break;
      const auto pos_eq = line.find_first_of('=');
      if (pos_eq == std::string::npos) // not found
        continue;
      const auto keyword = line.substr(0, pos_eq);
      // Important: Strip trailing whitespace to avoid hard-to-detect problems!
      const auto value   = strip_trailing_whitespace(line.substr(pos_eq+1));
      if (parsed_params.count(keyword))
        throw std::runtime_error("Duplicate keyword: " + keyword);
      parsed_params[keyword] = value;
    }
  }
  return parsed_params;
}

// Locate block [name] in a file stream. Returns true if succeessful.
inline bool find_block(std::istream &F, const std::string &s) {
  std::string target = "[" + s + "]";
  F.clear();
  F.seekg(0, std::ios::beg);
  while (F) {
    if (auto l = nextline(F))
      if(target.compare(l.value()) == 0) { return true; }
  }
  return false;
}

// Parse the [param] block of an input file.
inline auto parser(const std::string &filename, const std::string &block) {
  auto F = safe_open_for_reading(filename);
  if (!find_block(F, block))
    throw std::runtime_error(fmt::format("Block {} not found in input file {}.", block, filename));
  return parse_block(F);
}

// Simple tokenizer class
class string_token {
 private:
   std::istringstream iss;
   std::set<std::string> l;
 public:
   explicit string_token(const std::string &s) : iss(s),
     l(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()) {};
   [[nodiscard]] auto find(const std::string &x) const { return l.count(x) != 0; }
};

// Skip comment lines in the input stream 'f'.
inline void skip_comments(std::istream &f, const bool output = false, std::ostream &OUT = std::cout) {
  while (f) {
    const auto ch = f.peek();
    // skip white space and line breaks
    if (ch == ' ' || ch == '\t' || ch == '\n') {
      f.ignore();
      continue;
    }
    // break the loop if a non-comment like found
    if (ch != '#') break;
    f.ignore(); // ignore '#'
    std::string line;
    std::getline(f, line);
    if (output) OUT << ">> " << line << std::endl;
  }
}

// Sort according to the first component of the pair. Second component is ignored (unlike in the default sort
// function).
struct sortfirst {
  template <typename T1, typename T2> constexpr bool operator()(const std::pair<T1, T2> &xy1, const std::pair<T1, T2> &xy2) { return xy1.first < xy2.first; }
};

template<typename T> auto range0(const T b) { return boost::irange(T{0}, b); }
template<typename T> auto range1(const T b) { return boost::irange(T{1}, b+1); }

// Returns true if the data file contains complex values
inline bool complex_data(const std::string &filename = "data") {
  std::ifstream F(filename);
  if (!F) throw std::runtime_error("Can't load initial data.");
  std::string l;
  std::getline(F, l);
  std::getline(F, l);
  std::getline(F, l); // third line
  const auto pos = l.find("COMPLEX");
  return pos != std::string::npos;
}

template<typename K, typename V>
auto vector_of_keys(const std::map<K,V> &container)
{
  std::vector<K> keys;
  for (const auto &k: container | boost::adaptors::map_keys)
    keys.push_back(k);
  return keys;
}

inline double atof(const std::string &s) { return std::atof(s.c_str()); }
inline int atoi(const std::string &s) { return std::atoi(s.c_str()); }

// Read data from stream F.
template <typename T1, typename T2>
std::vector<std::pair<T1, T2>> readtable(const std::string &filename, const bool verbose = false)
{
  auto F = safe_open_for_reading(filename);
  std::vector<std::pair<T1, T2>> v;
  while (F)
  {
      skip_comments(F);
      const auto x = read_one<T1>(F);
      const auto y = read_one<T2>(F);
      if (F.fail()) break;
      assert(std::isfinite(x) && std::isfinite(y));
      v.push_back(std::make_pair(x, y));
  }
  if (verbose) std::cout << v.size() << " lines read." << std::endl;
  return v;
}

template <typename T1, typename T2>
void writetable(const std::vector<std::pair<T1, T2>> &re, std::string filename, const int output_precision = 16) {
  auto F = safe_open(filename);
  F << std::setprecision(output_precision);
  for (const auto & [x, y] : re) F << x << " " << y << std::endl;
}

inline std::vector<std::string> split (const std::string &s, char delim) {
    std::vector<std::string> result;
    std::stringstream ss (s);
    std::string item;

    while (getline (ss, item, delim)) {
        result.push_back (item);
    }

    return result;
}

} // namespace

#endif
