// tdavg - Averaging with interpolation for thermodynamics
// Part of "NRG Ljubljana", Rok Zitko, rok.zitko@ijs.si, Aug 2009

// CHANGE LOG
// 28.8.2009 - first version
// 21.4.2010 - copy comment line

#define PROGRAM "tdavg"
#define DESCRIPTION "thermodynamics averaging tool"
#define VERSION "0.2"
#define USAGE "[input] [reference]"
#define AUTHOR "Rok Zitko, rok.zitko@ijs.si, 2009"

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <cstring>
#include <sys/stat.h>

#include <unistd.h>
#include <getopt.h>

using namespace std;

bool verbose     = false; // output verbosity level
bool veryverbose = false; // horribly detailed output

string name;     // input file
string ref_name; // reference file (to be subtracted element by element)

bool copycomments      = false;
string lastcommentline = "";

using dvec = vector<double>;

typedef pair<double, dvec> Line;
using multiVec = vector<Line>;

typedef pair<double, double> dpair;
using Vec = vector<dpair>;

using DataType = vector<multiVec>;
DataType input;     // input data
DataType reference; // reference data (to be subtracted element by element)
dvec mesh;          // output mesh

unsigned int Nz;          // number of data sets (twist parameters z)
unsigned int columns = 0; // number of columns (excluding the first column!)

const int OUTPUT_PRECISION = 16;

void usage(ostream &F = cout) { F << "Usage: " << PROGRAM << " " << USAGE << endl; }

// Check for the existance of a regular file.
bool file_exists(string filename) {
  struct stat s{};
  int result = stat(filename.c_str(), &s);
  if (result == 0 && S_ISREG(s.st_mode)) return true;
  return false;
}

void cmd_line(int argc, char *argv[]) {
  char c;
  while ((c = getopt(argc, argv, "vVc")) != -1) {
    switch (c) {
      case 'v': verbose = true; break;
      case 'V': veryverbose = true; break;
      case 'c': copycomments = true; break;
      default: abort();
    }
  }
  int remaining = argc - optind; // arguments left
  if (remaining > 2) {
    usage();
    exit(1);
  }
  if (remaining >= 1) {
    name = string(argv[optind]); // Input filename
  } else {
    name = string("td.dat");
  }
  if (remaining >= 2) {
    ref_name = string(argv[optind]); // Reference filename
  } else {
    const string ref_name_default = "td-ref.dat";
    if (file_exists("td-ref.dat")) {
      ref_name = ref_name_default;
    } else {
      ref_name = string("");
    }
  }
}

// Split a string containing double floating point number into a vector.
// All empty spaces are stripped.
dvec split(const string &s) {
  dvec elements;
  string::const_iterator i = s.begin();
  // Skip to first non-space
  while (i != s.end() && isspace(*i)) i++;
  while (i != s.end()) {
    string substring;
    while (i != s.end() && !isspace(*i)) {
      substring += *i;
      i++;
    }
    elements.push_back(atof(substring.c_str()));
    // Skip to next non-space
    while (i != s.end() && isspace(*i)) i++;
  }
  return elements;
}

ostream &operator<<(ostream &os, const dvec &v) {
  for (double i : v) { os << i << " "; }
  return os;
}

ostream &operator<<(ostream &os, const Vec &v) {
  for (const auto & i : v) { os << i.first << " " << i.second << endl; }
  return os;
}

DataType load(const string &filename) {
  ifstream f(filename.c_str());
  if (!f.good() || f.eof() || !f.is_open()) {
    cerr << "Error opening file " << filename << endl;
    exit(1);
  }
  cout << "Reading " << filename << endl;
  DataType result;
  // skip comment lines and locate the first data block
  while (f.good() && f.peek() == '#') {
    string line;
    getline(f, line);
    if (veryverbose) cout << "##" << line << endl;
    if (copycomments) lastcommentline = line;
  }
  do {
    multiVec data; // contains data for one block
    while (f.good() && f.peek() != '#') {
      string line;
      getline(f, line);
      if (veryverbose) cout << line << endl;
      if (f.fail()) break;
      dvec elements = split(line);
      if (elements.size() == 0) {
        cerr << "Error parsing line: " << line << endl;
        exit(1);
      }
      if (columns == 0) {
        columns = elements.size() - 1;
      } else {
        if (elements.size() - 1 != columns) {
          cerr << "Unequal number of columns in line: " << endl << line << endl;
          exit(1);
        }
      }
      double T    = elements[0]; // temperature
      dvec values = dvec(elements.begin() + 1, elements.end());
      if (veryverbose) cout << T << " " << values << endl;
      data.push_back(make_pair(T, values));
    }
    if (data.size() != 0) {
      if (verbose) cout << "Block length " << data.size() << endl;
      result.push_back(data); // accumulate results for all blocks
    }
    // skip comments to next block
    while (f.good() && f.peek() == '#') {
      string line;
      getline(f, line);
      if (veryverbose) { cout << "#" << line << endl; }
    }
  } while (f.good() && !f.eof());
  return result;
}

// Load all the input data.
void read_files() {
  input = load(name);
  Nz    = input.size();
  if (Nz == 0) {
    cerr << "Failed loading the input data!" << endl;
    exit(1);
  }
  if (ref_name != "") {
    reference           = load(ref_name);
    unsigned int Nz_ref = reference.size();
    if (Nz != Nz_ref) {
      cerr << "Reference data sets do not match the input!" << endl;
      exit(1);
    }
  }
  if (verbose) cout << "Nz=" << Nz << " columns=" << columns << endl;
}

// Linear interpolation class
class LinInt {
  protected:
  Vec vec;               // tabulated data
  int len{};               // length of vec
  int index{};             // index of the interval where last x was found
  double x0{}, x1{};         // last x was in [x0:x1]
  double f0{}, f1{};         // f(x0), f(x1)
  double deriv{};          // (f1-f0)/(x1-x0)
  bool newintegral_flag{}; // set to true when we switch to a new interval
  double xmin{}, xmax{};     // lowest and highest x contained in vec
  double fxmin{}, fxmax{};   // f(xmin), f(xmax)
  public:
  LinInt()= default;
  LinInt(Vec &in_vec) : vec(in_vec) {
    len              = vec.size();
    index            = -1;
    newintegral_flag = false;
    xmin             = vec.front().first;
    fxmin            = vec.front().second;
    xmax             = vec.back().first;
    fxmax            = vec.back().second;
  };
  void findindex(double x);
  double operator()(double x);
};

// Serach for the interval so that x is contained in [x0:x1].
void LinInt::findindex(double x) {
  if (index == -1) {
    // When interpolation class object is constructed, index is
    // initialized to -1. We have no prior knowledge, thus we search in
    // the full interval.
    for (int i = 0; i < len - 1; i++) {
      x0 = vec[i].first;
      x1 = vec[i + 1].first;
      if (x0 <= x && x <= x1) {
        index = i;
        break;
      }
    }
  } else {
    if (x >= x1) {
      for (int i = index + 1; i < len - 1; i++) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else {
      for (int i = index - 1; i >= 0; i--) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    }
  }
  if (!(0 <= index && index < len && x0 <= x && x <= x1)) {
    cerr << "findindex() error."
         << " x=" << x << " x0=" << x0 << " x1=" << x1 << " index=" << index << " len=" << len << endl;
    exit(1);
  }
  f0             = vec[index].second;     // f(x0)
  f1             = vec[index + 1].second; // f(x1)
  double Delta_y = f1 - f0;
  double Delta_x = x1 - x0;
  deriv          = Delta_y / Delta_x;
}

// Return y(x) using linear interpolation between the tabulated values.
double LinInt::operator()(double x) {
  // Extrapolate if necessary
  if (x <= xmin) return fxmin;
  if (x >= xmax) return fxmax;
  if (index == -1 || !(x0 <= x && x < x1)) {
    newintegral_flag = true;
    findindex(x);
  }
  double dx = x - x0;
  return f0 + deriv * dx;
}

const double EPS = 1e-6;

inline bool eq_approx(double a, double b) { return abs((a - b) / a) < EPS; }

dvec merge_meshes() {
  dvec mesh;
  for (unsigned int i = 0; i < Nz; i++) {
    const unsigned int len = input[i].size();
    for (unsigned int j = 0; j < len; j++) { mesh.push_back(input[i][j].first); }
  }
  sort(mesh.begin(), mesh.end());
  auto new_end = unique(mesh.begin(), mesh.end(), eq_approx);
  mesh.erase(new_end, mesh.end());
  return mesh;
}

using LinIntVector = vector<LinInt>;
using IntType = vector<LinIntVector>;
IntType f, reff; // interpolation objects

IntType interpolate(const DataType &data) {
  IntType result;
  for (unsigned int i = 0; i < Nz; i++) {
    LinIntVector f0;
    for (unsigned int j = 0; j < columns; j++) {
      unsigned int len = data[i].size();
      Vec v(len);
      for (unsigned int k = 0; k < len; k++) v[k] = make_pair(data[i][k].first, data[i][k].second[j]);
      sort(v.begin(), v.end()); // important!
      f0.push_back(LinInt(v));
    }
    result.push_back(f0);
  }
  return result;
}

multiVec avg() {
  multiVec result;
  const unsigned int len = mesh.size();
  for (unsigned int j = 0; j < len; j++) {
    const double T = mesh[j]; // temperature
    dvec sum(columns, 0.0);   // vector of column sums
    for (unsigned int k = 0; k < columns; k++) {
      for (unsigned int i = 0; i < Nz; i++) {
        sum[k] += f[i][k](T) / Nz;
        if (ref_name != "") sum[k] -= reff[i][k](T) / Nz;
      }
    }
    result.push_back(make_pair(T, sum));
  }

  return result;
}

void save(const multiVec &result, ostream &f) {
  const unsigned int len = result.size();
  for (unsigned int j = 0; j < len; j++) {
    f << result[j].first << " ";
    assert(result[j].second.size() == columns);
    for (unsigned int k = 0; k < columns; k++) f << result[j].second[k] << " ";
    f << endl;
  }
}

void hello() {
  cout << PROGRAM << " - " << DESCRIPTION << " - " << VERSION << endl;
  cout << AUTHOR << endl;
}

int main(int argc, char *argv[]) {
  hello();
  cout << setprecision(OUTPUT_PRECISION);
  cmd_line(argc, argv);
  const string output_name = (name == "td.dat" ? "td-avg.dat" : name + ".avg.dat");
  ofstream OF(output_name.c_str());
  if (!OF) {
    cerr << "Failed opening " << output_name << endl;
    exit(1);
  }
  OF << setprecision(OUTPUT_PRECISION);
  read_files();
  mesh = merge_meshes();
  f    = interpolate(input);
  if (ref_name != "") reff = interpolate(reference);
  multiVec result = avg();
  if (copycomments && lastcommentline != "") OF << lastcommentline << endl;
  save(result, OF);
}
