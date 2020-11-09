// misc.h - Miscelaneous functions
// Copyright (C) 2005-2020 Rok Zitko

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
#include <cstring> // stdcasecmp

#include <boost/range/irange.hpp>
#include <boost/range/adaptor/map.hpp>

#define FMT_HEADER_ONLY
#include <fmt/format.h>

#include "basicio.hpp"
#include "portabil.hpp"

namespace NRG {

template<typename T> auto get_back(T &d) { // usually T is list or deque
  my_assert(!d.empty());
  auto i = d.back();
  d.pop_back();
  return i;
}

template<typename T> auto get_front(T &d) {
  my_assert(!d.empty());
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
  my_assert_not_reached();
}

// Get next line from stream F, skipping empty lines and comments.
inline std::optional<std::string> nextline(std::ifstream &F) {
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
inline auto parse_block(std::ifstream &F) {
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
inline bool find_block(std::ifstream &F, const std::string &s) {
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
   const std::set<std::string> l;
 public:
   explicit string_token(std::string s) : iss(s), 
     l(std::istream_iterator<std::string>(iss), std::istream_iterator<std::string>()) {};
   [[nodiscard]] auto find(const std::string &x) const { return l.count(x) != 0; }
};

// Skip comment lines in the input stream 'f'.
inline void skip_comments(std::istream &f, const bool output = false, std::ostream &OUT = std::cout) {
  while (f) {
    char ch = f.peek();
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
  template <typename T1, typename T2> bool operator()(const std::pair<T1, T2> &xy1, const std::pair<T1, T2> &xy2) { return xy1.first < xy2.first; }
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

} // namespace

#endif
