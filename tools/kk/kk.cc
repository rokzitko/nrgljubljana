// Kramers-Kronig transformation tool
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si, 2007-2017

// The input file must consist of a table of space-separated (energy,
// value) pairs. The energy grid must be symmetric with respect to zero and
// the file must contain an even number of lines. Gauss-Kronrod quadrature
// rules are used. At singularity points the derivative is computed using
// GSL interpolation routines (cubic splines).

// CHANGE LOG
// 12.10.2007 - first version
// 17.10.2007 - combination of gsl_integ and naive trapezoid rule
// 18.10.2007 - extrapolation beyond the [Xmin:Xmax] limits -> no good!
//            - correction term -> better!
// 22. 2.2008 - code cleanup
//            - correction term that is actually correct -> much better!
//  7. 7.2008 - support for comment lines commencing with '#'
//  15.5.2009 - GSL integration -> significant improvement in the accuracy!
//            - code simplification and cleanup
//  12.1.2010 - missing headers added

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

using namespace std;

typedef pair<double, double> XYPOINT;
using XYFUNC = vector<XYPOINT>;
using DVEC = vector<double>;

// number of digits of precision in the generated output file
#define OUTPUT_PRECISION 16

// Output width for verbose mode
#define WIDTH 24

// Dump additional information to stdout?
bool verbose     = false;
bool veryverbose = false;

// mode==FILES: we read from a file and write to a file
// mode==STD: we read from stdin and write to stdout. No other output
// is sent to stdout (but errors are reported to stderr).
enum MODE { UNDEF, FILES, STD } mode;

// Read data from stream F.
void readim(istream &F, XYFUNC &v) {
  while (F) {
    if (F.peek() == '#') { // skip comment lines
      string line;
      getline(F, line);
    } else {
      double x, y;
      F >> x >> y;

      if (F.fail()) break;

      assert(isfinite(x) && isfinite(y));
      v.push_back(make_pair(x, y));
    }
  }

  if (mode == FILES) cout << v.size() << " lines read." << endl;
}

int len;           // number of data points
int nr;            // =len/2
double Xmin, Xmax; // the interval boundaries
double sum;        // integral over input function on [Xmin:Xmax] interval

DVEC Xpts, Ypts;
DVEC Xpos; // Only positive X points [grid]

gsl_interp_accel *acc;
gsl_spline *spline;

const size_t limit = 1000;
gsl_integration_workspace *w;

const double EPSABS = 1e-12; // numeric integration epsilon (absolute)
const double EPSREL = 1e-8;  // numeric integration epsilon (relative)

// Initialize the KK transformer
void init(XYFUNC &im) {
  sort(im.begin(), im.end());
  len = im.size();
  assert(len % 2 == 0);
  nr   = len / 2;
  Xmin = im[0].first;
  Xmax = im[len - 1].first;
  if (mode == FILES) cout << "Range: [" << Xmin << " ; " << Xmax << "]" << endl;
  if (gsl_fcmp(-Xmin, Xmax, 1.e-8) != 0) {
    cerr << "ERROR (KK::KK) - Only symmetric intervals are supported!\n";
    exit(1);
  }
  // Xpts are increasing
  Xpts = DVEC(len);
  Ypts = DVEC(len);
  for (int i = 0; i < len; i++) {
    Xpts[i] = im[i].first;
    Ypts[i] = im[i].second;
  }
  // NOTE: With akime splines there are problems with the loss of the
  // floating point precision in the numeric integration step. We use
  // cubic splines instead where no such difficulties appear.
  // NOTE (Dec 2017): monotonic cubic interpolation would be a better choice!
  acc = gsl_interp_accel_alloc();
  //const gsl_interp_type * Interp_type = gsl_interp_linear;
  //const gsl_interp_type * Interp_type = gsl_interp_cspline;
  const gsl_interp_type *Interp_type = gsl_interp_akima;
  spline                             = gsl_spline_alloc(Interp_type, len);
  gsl_spline_init(spline, &Xpts[0], &Ypts[0], len);
  sum = gsl_spline_eval_integ(spline, Xmin, Xmax, acc);
  if (!isfinite(sum)) {
    cerr << "Error: Integral is not a finite number." << endl;
    exit(1);
  }
  if (mode == FILES) cout << "Sum=" << sum << endl;
  for (int i = nr; i < len; i++) {
    // Check for the symmetry of the grid
    assert(gsl_fcmp(Xpts[len - i - 1], -Xpts[i], 1e-8) == 0);
  }
  // Xpos are positive and increasing!
  Xpos = DVEC(nr);
  copy(Xpts.begin() + nr, Xpts.end(), Xpos.begin());
  w = gsl_integration_workspace_alloc(limit);
  gsl_set_error_handler_off();
}

// Integrand. Method: we take a sum of the contributions for positive and
// negative x and return their sum. This is helpful for even integrands,
// since it leads to possible cancellations and better accuracy of the
// final result, in particular for small Z.
inline double f(double X, double Z) {
  double a, b;
  // [ f(x) - f(z) ] / (x-z)
  if (X != Z)
    a = (gsl_spline_eval(spline, X, acc) - gsl_spline_eval(spline, Z, acc)) / (X - Z);
  else
    a = gsl_spline_eval_deriv(spline, X, acc);
  // [ f(-x) - f(z) ] / (-x-z)
  if (-X != Z)
    b = (gsl_spline_eval(spline, -X, acc) - gsl_spline_eval(spline, Z, acc)) / (-X - Z);
  else
    b = gsl_spline_eval_deriv(spline, -X, acc);
  return a + b;
}

// Wrapper function for integration using GSL
inline double f_gsl(double X, void *params) {
  const double Z = *(double *)params;
  return f(X, Z);
}

void handle_qag(int status) {
  if (status && veryverbose) cerr << "WARNING - qag error: " << status << " -- " << gsl_strerror(status) << endl;
}

// Perform the calculation for one point. Note: this is the critical part
// of the code, both for computational requirements as for the accuracy of
// the results. Optimise wisely!
double calc(double Z) {
  // NOTE about Gauss-Kronrod: The higher-order rules give better accuracy
  // for smooth functions, while lower-order rules save time when the
  // function contains local difficulties, such as discontinuities. [GSL manual]
  // On each iteration the adaptive integration strategy bisects the
  // interval with the largest error estimate. The subintervals and their
  // results are stored in the memory provided by workspace. The maximum
  // number of subintervals is given by limit, which may not exceed the
  // allocated size of the workspace. [GSL manual]
  gsl_function F;
  F.function = &f_gsl;
  F.params   = &Z;
  double integral;
  double integration_error;
  int status = gsl_integration_qag(&F,   // integrand function
                                   0,    // lower integration boundary
                                   Xmax, // upper integration boundary
                                   EPSABS, EPSREL,
                                   limit,               // size of workspace w
                                   GSL_INTEG_GAUSS15,   // Gauss-Kronrod rule
                                   w,                   // integration workspace
                                   &integral,           // final approximation
                                   &integration_error); // estimate of absolute error
  handle_qag(status);
  // Add an approximation of the (-inf,-Xmax] and [Xmax,+inf) intervals.
  double correction = 0.0;
  if (fabs(Z) != Xmax) correction = -gsl_spline_eval(spline, Z, acc) * 2. * gsl_atanh(Z / Xmax);
  double sum = integral + correction;
  if (mode == FILES && verbose) {
    cout << scientific;
    cout << setw(WIDTH) << Z << " t=" << setw(WIDTH) << integral / M_PI << " c=" << setw(WIDTH) << correction / M_PI << " ratio=" << setw(WIDTH)
         << correction / integral << endl;
    cout.unsetf(ios_base::scientific);
  }
  sum /= M_PI; // Don't forget to divide by pi!
  return sum;
}

// Perform the calculations for all points.
void calc(XYFUNC &re) {
  re.clear();
  re.reserve(len);
  for (int i = 0; i < len; i++) {
    double Z = Xpts[i];
    re.push_back(make_pair(Z, calc(Z)));
  }
}

void writere(XYFUNC &re, ostream &F) {
  F << setprecision(OUTPUT_PRECISION);
  for (auto & i : re) F << i.first << " " << i.second << endl;
}

void done() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(w);
}

void about() {
  cout << "Kramers-Kronig transformation tool - GSL" << endl;
#ifdef __TIMESTAMP__
  cout << "Timestamp: " << __TIMESTAMP__ << endl;
#endif
  cout << "Compiled on " << __DATE__ << " at " << __TIME__ << endl;
}

void usage() {
  cout << "\nUsage: kk <input> <output>\n";
  cout << "\nAlternative usage: kk -\n";
  cout << "\nIn this mode, kk reads from STDIN and outputs to STDOUT." << endl;
}

int main(int argv, char *argc[]) {
  mode = UNDEF;
  if (argv == 3) mode = FILES;
  if (argv == 2 && strcmp(argc[1], "-") == 0) mode = STD;
  if (mode != STD) about();
  if (mode == UNDEF) {
    usage();
    exit(1);
  }
  ifstream Fin;
  ofstream Fout;
  if (mode == FILES) {
    char *inputfn  = argc[1];
    char *outputfn = argc[2];
    cout << inputfn << " --> " << outputfn << endl;
    Fin.open(inputfn);
    if (!Fin) {
      cerr << "Can't open " << inputfn << " for reading." << endl;
      exit(2);
    }
    Fout.open(outputfn);
    if (!Fout) {
      cerr << "Can't open " << outputfn << " for writing." << endl;
      exit(2);
    }
  }
  cout << setprecision(OUTPUT_PRECISION);
  XYFUNC im, re;
  readim((mode == FILES ? Fin : cin), im);
  init(im);
  calc(re);
  writere(re, (mode == FILES ? Fout : cout));
  done();
  cout << "KK done!" << endl;
}
