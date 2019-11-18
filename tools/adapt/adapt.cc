// Discretization ODE solver for NRG
// Adaptable discretization mesh code
//
// Linear interpolation of rho(omega), integration using trapezoidal
// method. 4th-order Runge-Kutta ODE solver. Secant method for
// refinement of parameter A.

/*
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

// CHANGE LOG
// 18.6.2009 - bug fix for Gamma(1)=0 case
// 5.3.2010 - parameters 'hardgap' and 'boundary' for specifying an excluded
//            interval around omega=0
// 30.10.2011 - parameter 'max_abs' for maximal value of |f(x)|; if exceeded,
//              the program aborts. Default is 100.
// 19.9.2016 - comments added, some code reorganisation in timestep()
//           - max_subdiv default decreased to 10
// 23.9.2016 - bandrescale support

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <ctime>

using namespace std;

#include "lambda.h"
#include "linint.h"
#include "io.h"
#include "parser.h"
#include "load.h"
#include "calc.h"

string param_fn = "param"; // file with input parameters
SIGN sign       = POS;
LAMBDA Lambda; // discretization parameter

Vec vecrho; // GLOBAL: rho in tabulated form
// GLOBAL: interpolation objects
LinInt rho;
IntLinInt intrho1, intrho2;
// NOTE: There are two interpolation objects for the integral of rho: one
// for the lower boundary and one for the higher boundary. This leads to a
// significant performance improvement.
LinInt g;

// GLOBAL integration variables.
double x;                                          // GLOBAL: running x
double y;                                          // GLOBAL: running y(x)
double max_error;                                  // GLOBAL: maximum error in Delta(y)
double allowed_error;                              // GLOBAL: maximum error allowed
int max_subdiv;                                    // GLOBAL: maximum number of interval subdivisions
double xmax, xfine, output_step, dx_fine, dx_fast; // GLOBAL, described in set_parameters()
double convergence_eps = 1e-4;                     // Convergence epsilon for secant method
double factor0         = 1.0 + 1e-7;               // Shift of A in secant method
int max_iter           = 10;                       // Maximum number of iterations in the secant method
double max_abs         = 100.0;                    // Maximum value of |f(x)|.
double bandrescale     = 1.0;                      // Rescale the input data by this scale factor

bool adapt; // If adapt=false --> g(x)=1.

bool hardgap;
double boundary;

double intA; // intA=int_0^1 rho(w) dw (trapezoidal method).
double A;    // parameter in the shooting method. Initially A=intA.
             // Equal to intA if adapt=false.

// Right-hand-side of the differential equation. y=g !
double rhs_G(double x, double y) {
  double powL = Lambda.power(2.0 - x);
  return Lambda.logL() * (y - A / rho(y * powL));
}

void save(ostream &OUT) { OUT << x << " " << y << endl; }

bool check2 = false; // XXX

// Perform a step of length (at most) try_dx. Returns actual dx used.
double timestep(double try_dx, double (*rhs)(double, double)) {
  double dx = try_dx; // Current step length
  double dy, error;

  int subdivision = 0;

  // Approach: we attempt Runge-Kutta steps, reducing dx if necessary so
  // as to avoid crossing the tabulation data interval boundaries (with
  // too large steps) and to avoid too large errors.
  int cnt           = 0;
  const int max_cnt = 2 * max_subdiv;
  while (true) {
    assert(isfinite(dx));

    rho.clear_flag();

    double k1, k2, k3, k4; // Runge-Kutta intermediate results
    k1 = dx * rhs(x, y);
    k2 = dx * rhs(x + dx / 2.0, y + k1 / 2.0);
    k3 = dx * rhs(x + dx / 2.0, y + k2 / 2.0);
    k4 = dx * rhs(x + dx, y + k3);

    // Have we crossed the interval boundary ?
    if (rho.flag()) {
      if (subdivision < max_subdiv) {
        subdivision++;
        dx /= 2.0;
        continue;
      }
    }

    dy = (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;

    // Calculate a second-order result using Heun's formula and estimate
    // the error in y.
    double k2_heun = dx * rhs(x + dx, y + k1);
    double dy_heun = (k1 + k2_heun) / 2.0;
    error          = abs(dy_heun - dy);
    max_error      = max(max_error, error); // update max_error

    // Check the flags for k2_heun evaluation.
    if (rho.flag()) {
      if (subdivision < max_subdiv) {
        subdivision++;
        dx /= 2.0;
        continue;
      }
    }

    // Error too large? Then try again!
    if (error > allowed_error) {
      if (subdivision < max_subdiv) {
        subdivision++;
        dx /= 2.0;
        continue;
      }
      // If we cannot subdivide, we need to give up and break from the loop!
      break;
    }

    cnt++;
    if (cnt > max_cnt) {
      // Paranoid: we should never see this!
      cerr << "WARNING: cnt=" << cnt << " in timestep()" << endl;
    }

    break; // exit loop here: the RK step is appropriate
  }

  if (check2) {
    // Sanity checks
    const double dEdx      = dy / dx - y * Lambda.logL();
    const double tolerance = 1e-5;
    if (dEdx >= 0.0) { // E(x) must be monotonously decreasing!
      cerr << "WARNING: dE/dx is not negative." << endl;
      cerr << " x=" << x << " y=" << y << " dEdx=" << dEdx << " dy/dx=" << dy / dx << " y*log(Lambda)=" << y * Lambda.logL() << endl;
    }
    if (dEdx > tolerance) {
      cerr << "Aborting!" << endl;
      exit(1);
    }
  }

  x += dx;
  y += dy;

  return dx;
}

// Integrate with step dx up to xf.
void int_with_to(double dx, double xf, double (*rhs)(double, double)) {
  // Maximum discrepancies in x when approaching mesh points.
  const double max_diff_x = 1e-10;

  int cnt           = 0;
  const int max_cnt = 1000000;
  do {
    // If dx is too long, steplength should be decreased to xf-x.
    const double steplength = (x + dx < xf ? dx : xf - x);

    timestep(steplength, rhs);

    cnt++;
    if (cnt > max_cnt) {
      // Paranoid: we should never see this!
      cerr << "WARNING: cnt=" << cnt << " in int_with_to()" << endl;
      cerr << "Try decreasing max_subdiv!" << endl;
    }
  } while (abs(x - xf) > max_diff_x); // x is a global variable, updated in timestep()
}

// returns the g(xmax)/[A/rho(0)] ratio
double shoot_g(Vec &vecg) {
  string filename_g = "GSOL" + string(sign == POS ? "" : "NEG") + ".dat";
  ofstream OUTG;
  safe_open(OUTG, filename_g);

  cout << "#  A=" << A << endl;
  vecg.clear();

  rho(1); // NECESSARY!

  double dx = dx_fine; // initial step size

  // Initial conditions
  x         = 2.0;
  y         = 1.0; // y=g here!
  max_error = 0.0;

  save(OUTG); // save x=2 data point
  vecg.push_back(make_pair(x, y));

  double x_st = x; // Target x for next output line
  do {
    x_st += output_step;
    int_with_to(dx, x_st, rhs_G); // rhs_G !!
    save(OUTG);
    vecg.push_back(make_pair(x, y));
    if (x > xfine) { dx = dx_fast; }
  } while (x < xmax);

  double factor = A / rho(0);
  double ratio  = y / factor;
  cout << "#  x_last=" << x << " "
       << "g_last/factor=" << ratio << endl;
  cout << "#  eps_last=" << y * Lambda.power(2 - x) << " max_error=" << max_error << endl;
  return ratio;
}

void init_A() {
  intA = integrate_ab(vecrho, 0.0, 1.0);
  cout << "# intA=" << intA << endl;

  A = intA;
}

Vec calc_g(const Vec &vecrho) {
  Vec vecg;

  double ratio1 = shoot_g(vecg);
  if (abs(ratio1 - 1.0) < convergence_eps) { // We're done
    return vecg;
  }

  // Otherwise more effort is required...

  if (ratio1 > 1.0) { // intA too small
    A = A * factor0;
  } else { // intA too large
    A = A / factor0;
  }

  double ratio2 = shoot_g(vecg);
  if (abs(ratio2 - 1.0) < convergence_eps) { // We're done
    return vecg;
  }

  // Refine using secant method
  double x0 = intA;
  double x1 = A;
  double y0 = ratio1 - 1.0;
  double y1 = ratio2 - 1.0;
  int iter  = 0;
  double ynew;
  do {
    double xnew     = x1 - (x1 - x0) / (y1 - y0) * y1;
    A               = xnew;
    double rationew = shoot_g(vecg);
    ynew            = rationew - 1.0;
    // Shift
    x0 = x1;
    y0 = y1;
    x1 = xnew;
    y1 = ynew;
    iter++;
  } while (abs(ynew) > convergence_eps && iter < max_iter);

  if (abs(ynew) > convergence_eps) {
    cerr << "Secant method failed to converge in " << max_iter << " steps." << endl;
    exit(1);
  } else {
    cout << "# Converged!" << endl;
  }

  return vecg;
}

// Rescaling of for excluding finite intervals around omega=0. The
// new accumulation point is determined by the global variable
// 'boundary'.
double rescale(double omega) { return (1.0 - boundary) * omega + boundary; }

// eps(x) = D g(x) Lambda^(2-x) for x>2.
double eps(double x) {
  const double gx = (adapt ? g(x) : 1.0);
  double epsilon  = (x <= 2.0 ? 1.0 : gx * Lambda.power(2.0 - x));
  if (hardgap) { epsilon = rescale(epsilon); }
  return epsilon;
}

// Eps(x) = D f(x) Lambda^(2-x)
inline double Eps(double x, double f) {
  assert(x >= 1 && f > 0);
  return f * Lambda.power(2.0 - x);
}

// Right-hand-side of the differential equation. y=f !
double rhs_F(double x, double y) {
  assert(isfinite(x));
  assert(isfinite(y));

  const double term1 = Lambda.logL() * y;

  const double integral = intrho2(eps(x)) - intrho1(eps(x + 1));
  const double powL     = Lambda.power(2.0 - x);
  const double denom    = powL * rho(y * powL);

  double term2 = integral / denom;

  if (denom == 0.0) {
    cout << "# Warning: denom=0 with integral=" << integral << endl;
    cout << "# (denom=0 may arise from zeros in the hybridisation function)" << endl;
    term2 = 0.0;
  }

  assert(isfinite(term2));

  const double result = term1 - term2;

  return result;
}

void about(ostream &F = cout) {
  F << "# Discretization ODE solver" << endl;
  F << "# Rok Zitko, rok.zitko@ijs.si, 2008-2019" << endl;
}

void cmd_line(int argc, char *argv[]) {
  if (argc >= 2) {
    char first = toupper(argv[1][0]);
    switch (first) {
      case 'P': sign = POS; break;
      case 'N': sign = NEG; break;
      default: cerr << "Usage: adapt [P|N] [param_filename]" << endl; exit(1);
    }
  }

  if (argc == 3) { param_fn = string(argv[2]); }

  cout << "# ++ " << (sign == POS ? "POSITIVE" : "NEGATIVE") << endl;
}

void add_zero_point (Vec &vecrho)
{
  double x0 = vecrho.front().first;
  double y0 = vecrho.front().second;
  const double SMALL = 1e-99;
  if (x0 > SMALL)
    vecrho.push_back(make_pair(SMALL, y0));
  sort(begin(vecrho), end(vecrho));
}

void load_init_rho() {
  string rhofn = Pstr("dos", "Delta.dat");
  vecrho       = load_rho(rhofn, sign);
  add_zero_point(vecrho);
      
  rescalevecxy(vecrho, 1.0 / bandrescale, bandrescale);
  minmaxvec(vecrho, "rho");
  rho = LinInt(vecrho);
  cout << "# rho(0)=" << rho(0) << " rho(1)=" << rho(1) << endl;

  Vec vecintrho(vecrho);
  integrate(vecintrho);

  intrho1 = IntLinInt(vecrho, vecintrho);
  intrho2 = intrho1;
}

void set_parameters() {
  cout << setprecision(PREC);

  Lambda = LAMBDA(P("Lambda", 2.0));
  assert(Lambda > 1.0);

  adapt = Pbool("adapt", false); // Enable adaptable g(x)? Default is false!!

  hardgap  = Pbool("hardgap", false); // Exclude an interval around omega=0 ?
  boundary = P("boundary", 0.0);      // The boundary of the exclusion interval.

  bandrescale = P("bandrescale", 1.0); // band rescaling parameter

  xmax = P("xmax", 30); // Integrate over [1..xmax]
  assert(xmax > 1);
  xfine = P("xfine", 5); // Fine stepsize integral [1..xfine]
  assert(xfine > 1);
  output_step = P("outputstep", 1.0 / 64.0); // Stepsize for output file
  assert(output_step <= 1.0);
  dx_fine       = P("dx_fine", 1e-5);        // Integration stepsize in [1..xfine]
  dx_fast       = P("dx_fast", 1e-4);        // Integration stepsize in [xfine..xmax]
  allowed_error = P("allowed_error", 1e-10); // error control for adaptable stepsize
  max_subdiv    = Pint("max_subdiv", 10);    // maximum nr of integ. step subdivisions
  assert(dx_fine * pow(0.5, max_subdiv) > DBL_EPSILON);
  max_abs = P("max_abs", 100.0); // Maximal |f(x)|
  assert(max_abs > 0.0);

  convergence_eps = P("secant_eps", 1e-4);
  factor0         = 1.0 + P("secant_factor", 1e-7);
  max_iter        = Pint("secant_max_iter", 10);

  cout << "# ++ " << (adapt ? "ADAPTIVE" : "FIXED-GRID") << endl;
  cout << "# Lambda=" << Lambda;
  cout << " bandrescale=" << bandrescale;
  cout << " xmax=" << xmax;
  cout << " xfine=" << xfine;
  cout << " output_step=" << output_step;
  cout << " dx_fine=" << dx_fine << " dx_fast=" << dx_fast;
  cout << endl;
  cout << "# allowed_error=" << allowed_error;
  cout << " max_subdiv=" << max_subdiv;
  cout << " max_abs=" << max_abs;
  cout << endl;
}

void load_or_calc_g() {
  Vec vecg;
  const bool loadg = Pbool("loadg", false);
  if (loadg) {
    string gfn = "GSOL" + string(sign == POS ? "" : "NEG") + ".dat";
    vecg       = load_g(gfn);
  } else {
    vecg = calc_g(vecrho);
  }

  minmaxvec(vecg, "g");
  g = LinInt(vecg);
}

void calc_f() {
  string filename_f = "FSOL" + string(sign == POS ? "" : "NEG") + ".dat";
  ofstream OUTF;
  safe_open(OUTF, filename_f);

  check2 = true; // Ensure that Eps(x) is strictly decreasing

  rho(1); // NECESSARY!

  double dx = dx_fine; // initial step size

  // Initial conditions
  x         = 1.0;
  y         = 1.0 / Lambda; // y=f here!
  max_error = 0.0;

  save(OUTF); // save x=1 data point

  double x_st = x; // Target x for next output line
  do {
    x_st += output_step;
    int_with_to(dx, x_st, rhs_F); // rhs_F !!
    save(OUTF);
    if (x > xfine) { dx = dx_fast; }
    if (abs(y) > max_abs) {
      cout << "***** y=" << y << " |y|>max_abs=" << max_abs << endl;
      cout << "***** Terminating!" << endl;
    }
  } while (x < xmax && abs(y) <= max_abs);

  double factor = Lambda.factor() * A / rho(0);
  cout << "# x_last=" << x << endl;
  cout << "# f_last=" << y << " [f_last/factor=" << y / factor << "]" << endl;
  cout << "# eps_last=" << eps(x) << " (smallest energy point considered [input])" << endl;
  cout << "# Eps_last=" << y * Lambda.power(2 - x) << " (smallest energy scale obtained [output])" << endl;
  cout << "# max_error=" << max_error << " (maximum integration error)" << endl;
}

int main(int argc, char *argv[]) {
  clock_t start_clock = clock();

  about();
  cmd_line(argc, argv);
  parser(param_fn);
  set_parameters();

  load_init_rho();
  init_A();

  if (adapt) { load_or_calc_g(); }

  calc_f();

  clock_t end_clock = clock();
  cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << endl;
}
