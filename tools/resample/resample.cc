// Resample a tabulated function on a new X grid using smooth interpolation functions (GSL).
// Rok Zitko, rok.zitko@ijs.si, May 2014

// The input file must consist of a table of space-separated (energy,
// value) pairs. Gauss-Kronrod quadrature rules are used.

// CHANGE LOG
// 22.5.2014 - first version

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <getopt.h>

using namespace std;

typedef pair<double, double> XYPOINT;
typedef vector<XYPOINT> XYFUNC;
typedef vector<double> DVEC;

// number of digits of precision in the output
#define OUTPUT_PRECISION 16

// Dump additional information to stdout?
bool verbose = false; // enable with -v

string inputfn; // Filename for input data
string gridfn; // Filename for a new X grid
string outputfn; // Filename for resampled data. May be the same as gridfn.

void about()
{
   cout << "Resampling tool" << endl;
#ifdef __TIMESTAMP__
   cout << "Timestamp: " << __TIMESTAMP__ << endl;
#endif
   cout << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
}

void usage()
{
   cout << "\nUsage: resample [-v] <input> <grid> <output>" << endl;
   cout << "-v: toggle verbose messages (now=" << verbose << ")" << endl;
}

void cmd_line(int argc, char *argv[])
{
   char c;
   while ((c = getopt(argc, argv, "v")) != -1) {
      switch (c) {
       case 'v':
	 verbose = true;
	 break;
       default:
	 abort();
      }
   }
   int remaining = argc-optind;
   if (remaining != 3) {
      about();
      usage();
      exit(1);
   }
   inputfn = string(argv[optind++]);
   gridfn = string(argv[optind++]);
   outputfn = string(argv[optind++]);
}

// Read data from stream F.
void readtable(istream &F, XYFUNC &v)
{
  while (F) {
    if (F.peek() == '#') { // skip comment lines
       string line;
       getline(F, line);
    } else {
       double x, y;
       F >> x >> y;
       if (F.fail())
	 break;    
       assert (isfinite(x) && isfinite(y));
       v.push_back(make_pair(x, y));
    }
  }
  if (verbose)
     cout << v.size() << " lines read." << endl;
}

int len; // number of data points
double Xmin, Xmax; // the interval boundaries 
  
DVEC Xpts, Ypts;
 
gsl_interp_accel *acc;
gsl_spline *spline;

const size_t limit = 1000;
gsl_integration_workspace *w;

const double EPSABS = 1e-12; // numeric integration epsilon (absolute)
const double EPSREL = 1e-8; // numeric integration epsilon (relative)
 
void init(XYFUNC &im)
{
   sort(im.begin(), im.end());
   len = im.size();
   Xmin = im[0].first;
   Xmax = im[len-1].first;
   if (verbose)
     cout << "Range: [" << Xmin << " ; " << Xmax << "]" << endl;
   // Xpts are increasing
   Xpts = DVEC(len);
   Ypts = DVEC(len);
   for (int i = 0; i < len; i++) {
      Xpts[i] = im[i].first;
      Ypts[i] = im[i].second;
   }
   acc = gsl_interp_accel_alloc();
   //const gsl_interp_type * Interp_type = gsl_interp_linear;
   //const gsl_interp_type * Interp_type = gsl_interp_cspline;
   const gsl_interp_type * Interp_type = gsl_interp_akima;
   spline = gsl_spline_alloc(Interp_type, len);
   gsl_spline_init(spline, &Xpts[0], &Ypts[0], len);
   w = gsl_integration_workspace_alloc(limit);
   gsl_set_error_handler_off();
}

void writetable(XYFUNC &re, ostream &F)
{
  F << setprecision(OUTPUT_PRECISION);
  for (XYFUNC::iterator i = re.begin(); i != re.end(); i++)
    F << i->first << " " << i->second << endl;
}

void resample(XYFUNC &grid)
{
   for (XYFUNC::iterator i = grid.begin(); i != grid.end(); i++)
     i->second = gsl_spline_eval(spline, i->first, acc);
}

void done()
{
   gsl_spline_free(spline);
   gsl_interp_accel_free(acc);
}

int main(int argv, char *argc[])
{
   cmd_line(argv, argc);
   if (verbose)
     about();
   ifstream Fin;
   Fin.open(inputfn.c_str());
   if (!Fin) {
      cerr << "Can't open " << inputfn << " for reading." << endl;
      exit(2);
   }
   XYFUNC f;
   readtable(Fin, f);
   init(f);
   ifstream Fgrid;
   Fgrid.open(gridfn.c_str());
   if (!Fgrid) {
      cerr << "Can't open " << gridfn << " for reading." << endl;
      exit(2);
   }
   XYFUNC grid;
   readtable(Fgrid, grid);
   resample(grid);
   ofstream Fout;
   Fout.open(outputfn.c_str());
   if (!Fout) {
      cerr << "Can't open " << outputfn << " for writing." << endl;
      exit(2);
   }
   writetable(grid, Fout);
   done();
}
