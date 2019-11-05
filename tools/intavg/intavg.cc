// intavg - Averaging with interpolation
// Part of "NRG Ljubljana", Rok Zitko, rok.zitko@ijs.si, Aug 2009

// CHANGE LOG
// 25.8.2009 - first version
// 12.1.2010 - missing header added

#define VERSION "0.1"

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

#include <unistd.h>
#include <getopt.h>

using namespace std;

bool verbose = false; // output verbosity level
bool veryverbose = false; // horribly detailed output
string name; // filename of binary files containing the raw data
int Nz; // Number of spectra (1..Nz)

typedef pair<double, double> Pair;
typedef vector<Pair> Vec;
typedef vector<double> dvec;

vector<Vec> input; // input data
dvec mesh; // output mesh

void usage(ostream &F = cout)
{
   F << "Usage: intavg <name> <Nz>\n";
}

void cmd_line(int argc, char *argv[])
{
   char c;
   
   while ((c = getopt(argc, argv, "vV")) != -1) {
      switch (c) {
       case 'v':
	 verbose = true;
	 break;
	 
       case 'V':
	 veryverbose = true;
	 break;
      
       default:
	 abort();
      }
   }
   
   int remaining = argc-optind; // arguments left
   
   if (remaining != 2) {
      usage();
      exit(1);
   }
   name = string(argv[optind]); // Name of spectral density files
   Nz = atoi(argv[optind+1]); // Number of z-values
   assert(Nz >= 1);
}

string tostring(int i)
{
   ostringstream S;
   S << i;
   return S.str();
}

Vec load(int i)
{
   // i-th z-value defines the name of the directory where the results of
   // the NRG calculation are contained.
   const string filename = tostring(i) + "/" + name;

   ifstream f(filename.c_str());
   if (!f.good() || f.eof() || !f.is_open()) {
      cerr << "Error opening file " << filename << endl;
      exit(1);
   }
   if (verbose) {
      cout << "Reading " << filename << endl;
   }

   Vec data;
   
   while (f) {
      if (f.peek() == '#') { // skip comment lines
	 string line;
	 getline(f, line);
      } else {
	 double x, y;
	 f >> x >> y;
	 
	 if (f.fail())
	   break;
	 
	 assert(isfinite(x) && isfinite(y));
	 data.push_back(make_pair(x, y));
      }
   }
   
   return data;
}

// Load all the input data.
void read_files()
{
   for (int i = 1; i <= Nz; i++) {
      Vec v = load(i);
      input.push_back(v);
   }
}

// Linear interpolation class
class LinInt
{
protected:
  Vec vec; // tabulated data
  int len; // length of vec
  int index; // index of the interval where last x was found
  double x0, x1; // last x was in [x0:x1]
  double f0, f1; // f(x0), f(x1)
  double deriv; // (f1-f0)/(x1-x0)
  bool newintegral_flag; // set to true when we switch to a new interval
  double xmin, xmax; // lowest and highest x contained in vec
  double fxmin, fxmax; // f(xmin), f(xmax)

public:
  LinInt() {};
  LinInt(Vec &in_vec) : vec(in_vec) {
    len = vec.size();
    index = -1;
    newintegral_flag = false;
    xmin = vec.front().first;
    fxmin = vec.front().second;
    xmax = vec.back().first;
    fxmax = vec.back().second;
  };
  void findindex(double x);
  double operator()(double x);
};

// Serach for the interval so that x is contained in [x0:x1].
void LinInt::findindex(double x)
{
  if (index == -1) {
    // When interpolation class object is constructed, index is
    // initialized to -1. We have no prior knowledge, thus we search in
    // the full interval.
    for (int i = 0; i < len-1; i++) {
      x0 = vec[i].first;
      x1 = vec[i+1].first;
      if (x0 <= x && x <= x1) {
        index = i;
        break;
      }
    }
  } else {
    if (x >= x1) {
      for (int i = index+1; i < len-1; i++) {
        x0 = vec[i].first;
        x1 = vec[i+1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else {
      for (int i = index-1; i >= 0; i--) {
        x0 = vec[i].first;
        x1 = vec[i+1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    }
  }

  if (!(0 <= index && index < len && x0 <= x && x <= x1)) {
    cerr << "findindex() error."
         << " x=" << x
         << " x0=" << x0
         << " x1=" << x1
         << " index=" << index
         << " len=" << len << endl;
    exit(1);
  }

  f0 = vec[index].second; // f(x0)
  f1 = vec[index+1].second; // f(x1)
  double Delta_y = f1-f0;
  double Delta_x = x1-x0;
  deriv = Delta_y/Delta_x;
}

// Return y(x) using linear interpolation between the tabulated values.
double LinInt::operator()(double x)
{
  // Extrapolate if necessary
  if (x <= xmin) {
    return fxmin;
  }
  if (x >= xmax) {
    return fxmax;
  }

  if (index == -1 || !(x0 <= x && x < x1)) {
    newintegral_flag = true;
    findindex(x);
  }

  double dx = x-x0;
  return f0 + deriv * dx;
}

const double EPS = 1e-6;

inline bool eq_approx(double a, double b)
{
   return abs((a-b)/a) < EPS;
}

dvec merge_meshes()
{
   dvec mesh;
   
   for (int i = 0; i < Nz; i++) {
      const int len = input[i].size();
      for (int j = 0; j < len; j++) {
	   mesh.push_back(input[i][j].first);
      }
   }
   
   sort(mesh.begin(), mesh.end());
   dvec::iterator new_end = unique(mesh.begin(), mesh.end(), eq_approx);
   mesh.erase(new_end, mesh.end());
   return mesh;
}

vector<LinInt> f; // interpolation objects

void interpolate()
{
   for (int i = 0; i < Nz; i++) {
      f.push_back(LinInt(input[i]));
   }
}

Vec avg()
{
   Vec result;
   
   const int len = mesh.size();
   for (int j = 0; j < len; j++) {
      const double x = mesh[j];
      double sum = 0.0;
      for (int i = 0; i < Nz; i++) {
	 sum += f[i](x);
      }
      const double avg = sum/Nz;
      result.push_back(make_pair(x, avg));
   }
   
   return result;
}

void save(const Vec &result, string filename)
{
   ofstream f(filename.c_str());
   if (!f) {
      cerr << "Failed opening " << filename << endl;
      exit(1);
   }
   f << setprecision(16);
   
   const int len = result.size();
   for (int j = 0; j < len; j++) {
      f << result[j].first << " " << result[j].second << endl;
   }
}

int main(int argc, char *argv[])
{
   cout << "intavg - Averaging tool - " << VERSION << endl;
   cout << "Rok Zitko, rok.zitko@ijs.si, 2009" << endl;
   cout << setprecision(16);

   cmd_line(argc, argv);
   read_files();
   mesh = merge_meshes();
   interpolate();
   Vec result = avg();
   save(result, name);
}
