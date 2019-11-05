// Discretization ODE solver for NRG
//
// ** Loading (and parsing) of tabulated data

#include <algorithm>

// Split a string 's' into substrings. Leading spaces are ignored.
vector<string> split_string(const string &s, unsigned int atleast = 0)
{
  int index = 0;
  int len = s.length();

  while (index < len && isspace(s[index])) {
    index++;
  }

  vector<string> substrings;

  while (index < len) {
    string substr = "";

    // Copy string until space or end of string
    while (index < len && !isspace(s[index])) {
      substr += s[index];
      index++;
    }
    substrings.push_back(substr);

    // Locate new substring
    while (index < len && isspace(s[index])) {
      index++;
    }
  }

  if (substrings.size() < atleast) {
    cerr << "ERROR: At least " << atleast << " columns expected." << endl;
    exit(1);
  }

  return substrings;
}

Vec load_g(const string &filename)
{
  ifstream F;
  safe_open(F, filename);

  Vec vecg;
  while (F) {
    string line = getnextline(F);
    if (!F)
      break;

    vector<string> columns = split_string(line, 2);

    double x = atof(columns[0]);
    double y = atof(columns[1]);

    vecg.push_back(make_pair(x, y));
  }
  if (vecg.size() == 0) {
    cerr << "ERROR: No data found." << endl;
    exit(1);
  }

  return vecg;
}

void rescalevecxy(Vec &vec, double factorx, double factory)
{
   const int len = vec.size();
   for (int i = 0; i < len; i++) {
      vec[i].first *= factorx;
      vec[i].second *= factory;
   }
   
   cout << "Rescaled to the interval [ " << vec.front().first 
        << " : " << vec.back().second << " ]" << endl;
}

// Show minimal and maximal y in a table.
void minmaxvec(Vec &vec, string name)
{
  double miny = DBL_MAX;
  double maxy = 0;
  const int len = vec.size();
  for (int i = 0; i < len; i++) {
    double y = vec[i].second;
    if (y > maxy) {
      maxy = y;
    }
    if (y < miny) {
      miny = y;
    }
  }

  cout << "# min[" << name << "]=" << miny;
  cout << " max[" << name << "]=" << maxy << endl;
}

enum SIGN { POS, NEG }; // positive vs. negative energies

// Load positive (sign=POS) or negative (sogn=NEG) part of the
// hybridisation function into a vector.
Vec load_rho(const string &filename, SIGN sign)
{
  ifstream F;
  safe_open(F, filename);

  Vec vecrho;
  while (F) {
    string line = getnextline(F);
    if (!F)
      break;

    vector<string> columns = split_string(line, 2);

    double x = atof(columns[0]);
    double y = atof(columns[1]);

    if ((sign == POS && x > 0) || (sign == NEG && x < 0)) {
      // y must be positive (or zero)
      if (y < 0.0) {
        cerr << "ERROR: Negative y found." << endl;
        exit(1);
      }

      // Disregard sign of x !
      vecrho.push_back(make_pair(abs(x), y));
    }
  }

  if (vecrho.size() == 0) {
    cerr << "ERROR: No data found." << endl;
    exit(1);
  }

  sort(vecrho.begin(), vecrho.end());
  cout << "# " << filename << " ";
  cout << "- " << (sign == POS ? "POS" : "NEG") << " ";
  cout << "- interval [ " << vecrho.front().first << " : ";
  cout << vecrho.back().first << " ]" << endl;

  return vecrho;
}

#include <algorithm>
#include <iterator>
#include <sstream>

string tostring(const Pair &p)
{
   ostringstream str;
   str << p.first << " " << p.second;
   return str.str();
}

void save(const string &fn, const Vec &v)
{
   ofstream F(fn.c_str());
   if (!F) {
      cerr << "Failed to open " << fn << " for writing." << endl;
      exit(1);
   }

   transform(v.begin(),
	     v.end(),
	     ostream_iterator<string>(F, "\n"),
	     tostring);
}

void save(const string &fn, const vector<double> &v)
{
   ofstream F(fn.c_str());
   if (!F) {
      cerr << "Failed to open " << fn << " for writing." << endl;
      exit(1);
   }
   
   F << setprecision(18);

   copy(v.begin(),
	v.end(),
	ostream_iterator<double>(F, "\n"));
}

void load(const string &fn, vector<double> &v)
{
   ifstream F;
   safe_open(F, fn);
   
   v.clear();
   
   while (F) {
      string line = getnextline(F);
      if (!F)
	break;

      const double x = atof(line.c_str());
      
      v.push_back(x);
   }
}
