// misc.h - Miscelaneous functions
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _misc_h_
#define _misc_h_

// Conversion functions
template <class T> T fromstring(const string &str) {
  T result;
  try {
    result = boost::lexical_cast<T>(str);
  } catch (boost::bad_lexical_cast &) { throw std::runtime_error(fmt::format("Lexical cast [{}] failed.", str)); }
  return result;
}
template <> bool fromstring(const string &str) { return (strcasecmp(str.c_str(), "true") == 0 ? true : false); }
// for T=int, std::to_string is used
template <class T> string to_string(const T val) { return boost::lexical_cast<string>(val); }

// switch statement with three cases
template <typename T, typename T1> T switch3(T1 x0, T1 x1, T y1, T1 x2, T y2, T1 x3, T y3) {
  if (x0 == x1) return y1;
  if (x0 == x2) return y2;
  if (x0 == x3) return y3;
  my_assert_not_reached();
}

// Get next line from stream F, skipping empty lines and comments.
string getnextline(ifstream &F) {
  string line;
  while (F) {
    getline(F, line);
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

void strip_trailing_whitespace(string &s) {
  string::reverse_iterator it = rbegin(s);
  while (it != rend(s) && isspace(*it)) {
    s.erase(--it.base());
    it = rbegin(s);
  }
}

// Parse a block of "keyword=value" lines.
auto parse_block(ifstream &F) {
  map<string, string> parsed_params; 
  while (F) {
    string line = getnextline(F);
    if (!F) break;
    if (line[0] == '[') // new block, we're done!
      break;
    string::size_type pos_eq = line.find_first_of('=');
    if (pos_eq == string::npos) // not found
      continue;
    const string keyword = line.substr(0, pos_eq);
    string value         = line.substr(pos_eq + 1);
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
bool find_block(ifstream &F, const string &s) {
  string target = "[" + s + "]";
  F.clear();
  F.seekg(0, ios::beg);
  while (F) {
    string line;
    getline(F, line);
    if (F && target.compare(line) == 0) { break; }
  }
  return !F.fail(); // True if found.
}

// Parse the [param] block of an input file.
map<string, string> parser(const string &filename, const string &block) {
  ifstream F = safe_open_for_reading(filename);
  if (find_block(F, block))
    return parse_block(F);
  return {};
}

// Input/output
template <typename T1, typename T2> ostream &operator<<(ostream &os, const pair<T1, T2> &p) { return os << p.first << ' ' << p.second; }
template <typename T> ostream &operator<<(ostream &os, const std::vector<T> &vec) {
  for (const auto &x : vec) os << x << " ";
  return os;
}
template <typename T> ostream &operator<<(ostream &os, const ublas::vector<T> &vec) {
  for (const auto &x : vec) os << x << " ";
  return os;
}
template <typename T> ostream &operator<<(ostream &os, const ublas::matrix<T> &m) {
  for (auto r1 = 0; r1 < m.size1(); r1++) {
    for (auto r2 = 0; r2 < m.size2(); r2++)
      os << m(r1, r2) << ' ';
    os << endl;
  }
  return os;
}

// Simple tokenizer class
class string_token {
  private:
  string s;
  list<string> l;
  public:
  explicit string_token(string _s) : s(std::move(_s)) {
    string::size_type pos = 0;
    string::size_type first, last;
    while ((first = s.find_first_not_of(" ", pos)) != string::npos) {
      last         = s.find_first_of(" ", first);
      string token = string(s, first, last - first);
      l.push_back(token);
      if (last == string::npos)
        break;
      else
        pos = last + 1;
    }
  }
  bool find(string x) const { return std::find(begin(l), end(l), x) != end(l); }
};

// Skip comment lines in the input stream 'f'.
void skip_comments(istream &f, bool output = false, ostream &OUT = cout) {
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
    string line;
    getline(f, line);
    if (output) OUT << ">> " << line << endl;
  }
}

// Sort according to the first component of the pair. Second
// component is ignored (unlike in the default sort function).
struct sortfirst {
  template <typename T1, typename T2> bool operator()(const pair<T1, T2> &xy1, const pair<T1, T2> &xy2) { return xy1.first < xy2.first; }
};

#endif
