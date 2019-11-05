// Numerical integration tool: compute accurate spectral weights and moments
// of tabulated functions using smooth interpolation functions (GSL).
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si, May 2014

// The input file must consist of a table of space-separated (energy,
// value) pairs. Gauss-Kronrod quadrature rules are used.

// CHANGE LOG
// 21.5.2014 - first version

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
bool veryverbose = false; // enable with -V
bool showwarnings = true; // disable with -s

double T = 1e-99; // Temperature. Default is (essentially) 0.
string inputfn; // Filename for input data

double sum; // integral over input function on [Xmin:Xmax] interval
double total, totalabs, pos, neg, fermi; // Results
string out = "total"; // What to output to STDOUT

void about()
{
   cout << "Integration tool - GSL" << endl;
#ifdef __TIMESTAMP__
   cout << "Timestamp: " << __TIMESTAMP__ << endl;
#endif
   cout << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
}

void usage()
{
   cout << "\nUsage: integ [-v] [-V] [-w] [-T temp] <input> [-p|n|a|f]" << endl;
   cout << "-v: toggle verbose messages (now=" << verbose << ")" << endl;
   cout << "-V: toggle very verbose messages (now=" << veryverbose << ")" << endl;
   cout << "-w: toggle warnings (now=" << showwarnings << ")" << endl;
   cout << "-T: temperature T" << endl;
   cout << "-p: integral over positive X range" << endl;
   cout << "-n: integral over negative X range" << endl;
   cout << "-a: integral over |f|" << endl;
   cout << "-f: integral weighted with Fermi-Dirac function for temperature T" << endl;
}

void cmd_line(int argc, char *argv[])
{
   char c;
   
   while ((c = getopt(argc, argv, "vVwt:T:pnaf")) != -1) {
      switch (c) {
       case 'v':
	 verbose = true;
	 break;
	          
       case 'V':
	 veryverbose = true;
	 break;
	 
       case 'w':
         showwarnings = false;
         break;
	          
       case 't': // case insensitive!
       case 'T': 
	 T = atof(optarg);
	 if (verbose) {
	    cout << "T=" << T << endl;
	 }
	 break;
	 
       case 'p':
	 out = "pos";
	 break;
	 
       case 'n':
	 out = "neg";
	 break;
	 
       case 'a':
	 out = "abs";
	 break;
	 
       case 'f':
	 out = "fermi";
	 break;
	 
       default:
	 abort();
      }
   }
   
   int remaining = argc-optind;
   
   if (remaining != 1) {
      about();
      usage();
      exit(1);
   }
   
   inputfn = string(argv[optind]);
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

   sum = gsl_spline_eval_integ(spline, Xmin, Xmax, acc);
   
   if (!isfinite(sum)) {
      cerr << "Error: Integral is not a finite number." << endl;
      exit(1);
   }
   
   if (verbose)
     cout << "Sum=" << sum << endl;

   w = gsl_integration_workspace_alloc(limit);
   
   gsl_set_error_handler_off();
}

inline double f_neg(double X, void * params)
{
   return (X < 0 ? gsl_spline_eval(spline, X, acc) : 0);
}

inline double f_pos(double X, void * params)
{
   return (X > 0 ? gsl_spline_eval(spline, X, acc) : 0);
}

inline double f_total(double X, void * params)
{
   return gsl_spline_eval(spline, X, acc);
}

inline double f_abs(double X, void * params)
{ 
  return fabs(gsl_spline_eval(spline, X, acc));
}

inline double f_fermi(double X, void * params)
{
   double fd = 1.0/(1.0+exp(X/T));
   return gsl_spline_eval(spline, X, acc) * fd;
}

void handle_qag(int status)
{
   if (status && showwarnings) {
      cerr << "WARNING - qag error: " << status
	<< " -- " << gsl_strerror(status) << endl;
   }
}

// NOTE about Gauss-Kronrod: The higher-order rules give better accuracy
// for smooth functions, while lower-order rules save time when the
// function contains local difficulties, such as discontinuities. [GSL manual]

// On each iteration the adaptive integration strategy bisects the
// interval with the largest error estimate. The subintervals and their
// results are stored in the memory provided by workspace. The maximum
// number of subintervals is given by limit, which may not exceed the
// allocated size of the workspace. [GSL manual]

double calc(double(*fnc)(double, void *))
{
   gsl_function F;
   F.function = fnc;
//   F.params = &Z;
   
   double integral;
   double integration_error;
   int status = gsl_integration_qag(&F, // integrand function
				    Xmin,  // lower integration boundary
				    Xmax, // upper integration boundary
				    EPSABS,
				    EPSREL,
				    limit, // size of workspace w
				    GSL_INTEG_GAUSS15, // Gauss-Kronrod rule
				    w, // integration workspace
				    &integral, // final approximation
				    &integration_error); // estimate of absolute error

   handle_qag(status);

   if (veryverbose) {
     cout << scientific;
     cout << "Result=" << integral << endl;
     cout << "Int. error=" << integration_error << endl;
     cout.unsetf(ios_base::scientific);
   }

   return integral;
}

void done()
{
   gsl_spline_free(spline);
   gsl_interp_accel_free(acc);
   gsl_integration_workspace_free(w);
}

int main(int argv, char *argc[])
{
   cmd_line(argv, argc);
   if (verbose)
     about();

   if (verbose)
     cout << "T=" << T << endl;
   
   ifstream Fin;

   Fin.open(inputfn.c_str());
   if (!Fin) {
      cerr << "Can't open " << inputfn << " for reading." << endl;
      exit(2);
   }

   XYFUNC f;
   readtable(Fin, f);
   init(f);
   total = calc(f_total);
   pos = calc(f_pos);
   neg = calc(f_neg);
   totalabs = calc(f_abs);
   fermi = calc(f_fermi);
   done();
   
   if (verbose) {
      cout << "Total=" << total << endl;
      cout << "Positive=" << pos << endl;
      cout << "Negative=" << neg << endl;
      cout << "Total|f|=" << totalabs << endl;
      cout << "Fermi-Dirac weighted=" << fermi << endl;
   }
   

   cout << setprecision(OUTPUT_PRECISION);
   if (out == "total")
     cout << total << endl;
   if (out == "pos")
     cout << pos << endl;
   if (out == "neg")
     cout << neg << endl;
   if (out == "abs")
     cout << totalabs << endl;
   if (out == "fermi") 
     cout << fermi << endl;
}
