// misc.h - Miscelaneous functions
// Copyright (C) 2005-2009 Rok Zitko

#ifndef _misc_h_
#define _misc_h_

// Conversion functions

template <class T> T fromstring(const string &str)
{
   T result;
   try {
      result = boost::lexical_cast<T>(str);
   }
   catch (boost::bad_lexical_cast &) {
      my_error("Lexical cast [%s] failed.", str.c_str());
   }
   return result;
}

template<> bool fromstring(const string &str)
{ 
   return (strcasecmp(str.c_str(), "true") == 0 ? true : false); 
}

template <class T> string tostring(const T val)
{
  return boost::lexical_cast<string>(val);
}

template<typename T, typename T1>
T switch3(T1 x0, T1 x1, T y1, T1 x2, T y2, T1 x3, T y3)
{
     if (x0 == x1)
         return y1;
     if (x0 == x2)
         return y2;
     if (x0 == x3)
         return y3;
     my_assert_not_reached();
     return 0.0; // avoid warnings
}

// Get next line from stream F, skipping empty lines and comments.
string getnextline(ifstream &F)
{
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

void strip_trailing_whitespace(string &s)
{
   string::reverse_iterator it = rbegin(s);
   while (it != rend(s) && isspace(*it)) {
      s.erase(--it.base());
      it = rbegin(s);
   }
}

// Parse a block of "keyword=value" lines.
void parse_block(map<string, string> &parsed_params, ifstream &F)
{
  while (F) {
    string line = getnextline(F);
    if (!F)
      break;
    if (line[0] == '[') // new block, we're done!
      break;
    string::size_type pos_eq = line.find_first_of('=');
    if (pos_eq == string::npos) // not found
      continue;
    const string keyword = line.substr(0, pos_eq);
    string value = line.substr(pos_eq+1);
    // Important: Strip trailing whitespace to avoid hard-to-detect problems!
    // (Note: INI parsers do this by convention!)
    strip_trailing_whitespace(value);
    if (parsed_params.count(keyword)) {
       cout << "Duplicate keyword: " << keyword << endl;
       exit(1);
    }
    parsed_params[keyword] = value;
  }
}

void safe_open(ifstream &F, const string &filename)
{
   F.open(filename.c_str());
   if (!F) 
      exit1("Can't open " << filename << " for reading.");
}

// Locate block [name] in a file stream. Returns true if succeessful.
bool find_block(ifstream &F, const string &s)
{
  string target = "[" + s + "]";
  F.clear();
  F.seekg(0, ios::beg);
  while (F) {
    string line;
    getline(F, line);
    if (F && target.compare(line) == 0) {
      break;
    }
  }
  return !F.fail(); // True if found.
}

void parser(map<string, string> &parsed_params, const string &filename)
{
  ifstream F;
  safe_open(F, filename);
  if (find_block(F, "param"))
    parse_block(parsed_params, F);
}

// Input/output

template <typename T1, typename T2>
ostream & operator<<(ostream &os, const pair<T1, T2> &p)
{
  return os << p.first << ' ' << p.second;
}

template <typename T>
ostream & operator<<(ostream &os, const std::vector<T> &vec)
{
   for(const auto &x : vec)
      os << x << " ";
   return os;
}

class string_token
{
 private:
  string s;
  list<string> l;

 public:
  string_token(const string &_s) : s(_s) {
    string::size_type pos = 0;
    string::size_type first, last;
    while ((first = s.find_first_not_of(" ", pos)) != string::npos) {
      last = s.find_first_of(" ", first);
      string token = string(s, first, last-first);
      l.push_back(token);
      if (last == string::npos)
        break;
      else
        pos = last+1;
    }
  }
  bool find(string x) const {
    return std::find(begin(l), end(l), x) != end(l);
  }
};

// Skip comment lines in the input stream 'f'.
void skip_comments(istream &f, bool output = false, ostream & OUT = cout)
{
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
struct sortfirst
{
   template <typename T1, typename T2>
      bool operator()(const pair<T1, T2> &xy1, const pair<T1, T2> &xy2) {
	 return xy1.first < xy2.first;
      }
};

#endif
