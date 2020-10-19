// misc.h - Miscelaneous functions
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _misc_h_
#define _misc_h_

// Conversion functions
template <class T> T fromstring(const std::string &str) {
  T result;
  try {
    result = boost::lexical_cast<T>(str);
  } catch (boost::bad_lexical_cast &) { throw std::runtime_error(fmt::format("Lexical cast [{}] failed.", str)); }
  return result;
}
template <> bool fromstring(const std::string &str) { return (strcasecmp(str.c_str(), "true") == 0 ? true : false); }
// for T=int, std::to_string is used
template <class T> std::string to_string(const T val) { return boost::lexical_cast<std::string>(val); }

// switch statement with three cases
template <typename T, typename T1> T switch3(const T1 x0, const T1 x1, const T y1, const T1 x2, const T y2, const T1 x3, const T y3) {
  if (x0 == x1) return y1;
  if (x0 == x2) return y2;
  if (x0 == x3) return y3;
  my_assert_not_reached();
}

// Get next line from stream F, skipping empty lines and comments.
std::string getnextline(std::ifstream &F) {
  std::string line;
  while (F) {
    std::getline(F, line);
    if (!F) // bail out
      break;
    if (line.length() == 0) // skip empty lines
      continue;
    if (line[0] == '#') // skip comment lines
      continue;
    return line;
  }
  return ""; // error
}

void strip_trailing_whitespace(std::string &s) {
  std::string::reverse_iterator it = s.rbegin();
  while (it != s.rend() && isspace(*it)) {
    s.erase(--it.base());
    it = s.rbegin();
  }
}

// Parse a block of "keyword=value" lines.
auto parse_block(std::ifstream &F) {
  std::map<std::string, std::string> parsed_params; 
  while (F) {
    std::string line = getnextline(F);
    if (!F) break;
    if (line[0] == '[') // new block, we're done!
      break;
    std::string::size_type pos_eq = line.find_first_of('=');
    if (pos_eq == std::string::npos) // not found
      continue;
    const std::string keyword = line.substr(0, pos_eq);
    std::string value         = line.substr(pos_eq + 1);
    // Important: Strip trailing whitespace to avoid hard-to-detect problems!
    // (Note: INI parsers do this by convention!)
    strip_trailing_whitespace(value);
    if (parsed_params.count(keyword))
      throw std::runtime_error("Duplicate keyword: " + keyword);
    parsed_params[keyword] = value;
  }
  return parsed_params;
}

// Locate block [name] in a file stream. Returns true if succeessful.
bool find_block(std::ifstream &F, const std::string &s) {
  std::string target = "[" + s + "]";
  F.clear();
  F.seekg(0, std::ios::beg);
  while (F) {
    std::string line;
    std::getline(F, line);
    if (F && target.compare(line) == 0) { break; }
  }
  return !F.fail(); // True if found.
}

// Parse the [param] block of an input file.
auto parser(const std::string &filename, const std::string &block) {
  auto F = safe_open_for_reading(filename);
  if (!find_block(F, block))
    throw std::runtime_error(fmt::format("Block {} not found in input file {}.", block, filename));
  return parse_block(F);
}

// Input/output
template <typename T1, typename T2> std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &p) { return os << p.first << ' ' << p.second; }
template <typename T> std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
  for (const auto &x : vec) os << x << " ";
  return os;
}
template <typename T> std::ostream &operator<<(std::ostream &os, const ublas::vector<T> &vec) {
  for (const auto &x : vec) os << x << " ";
  return os;
}
template <typename T> std::ostream &operator<<(std::ostream &os, const ublas::matrix<T> &m) {
  for (auto r1 = 0; r1 < m.size1(); r1++) {
    for (auto r2 = 0; r2 < m.size2(); r2++)
      os << m(r1, r2) << ' ';
    os << std::endl;
  }
  return os;
}

// Simple tokenizer class
class string_token {
 private:
   std::string s;
   std::list<std::string> l;
 public:
   explicit string_token(string _s) : s(std::move(_s)) {
     string::size_type pos = 0;
     string::size_type first, last;
     while ((first = s.find_first_not_of(" ", pos)) != std::string::npos) {
       last              = s.find_first_of(" ", first);
       std::string token = std::string(s, first, last - first);
       l.push_back(token);
       if (last == std::string::npos)
         break;
       else
         pos = last + 1;
     }
   }
   bool find(std::string x) const { return std::find(l.begin(), l.end(), x) != l.end(); }
};

// Skip comment lines in the input stream 'f'.
void skip_comments(std::istream &f, const bool output = false, std::ostream &OUT = std::cout) {
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

// Sort according to the first component of the pair. Second
// component is ignored (unlike in the default sort function).
struct sortfirst {
  template <typename T1, typename T2> bool operator()(const pair<T1, T2> &xy1, const pair<T1, T2> &xy2) { return xy1.first < xy2.first; }
};

#endif
