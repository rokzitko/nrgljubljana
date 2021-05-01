#ifndef _adapt_adapt_hpp_
#define _adapt_adapt_hpp_

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
#include <functional>

using namespace std;
using namespace std::string_literals;

#include "lambda.hpp"
#include "linint.hpp"
#include "io.hpp"
#include "parser.hpp"
#include "load.hpp"
#include "calc.hpp"

namespace NRG::Adapt {

inline void add_zero_point(Vec &v, const double small = 1e-99)
{
  const auto [x0,y0] = v.front();
  if (x0 > small)
    v.emplace_back(std::make_pair(small, y0));
  std::sort(begin(v), end(v));
}

class Adapt {
 public:
   Params P;
   Sign sign;     // positive or negative frequencies
   LAMBDA Lambda; // discretization parameter
   Vec vecrho;    // Density of states (rho) in tabulated form
   LinInt rho;    // Linear interpolation object
   IntLinInt intrho1, intrho2; // Integrated interpolation function objects
   LinInt g;
   double x;                                          // running x
   double y;                                          // running y(x)
   double max_error;                                  // maximum error in Delta(y)
   double allowed_error;                              // maximum error allowed
   int max_subdiv;                                    // maximum number of interval subdivisions
   double xmax, xfine, output_step, dx_fine, dx_fast; // described in set_parameters()
   double convergence_eps = 1e-4;                     // Convergence epsilon for secant method
   double factor0         = 1.0 + 1e-7;               // Shift of A in secant method
   int max_iter           = 10;                       // Maximum number of iterations in the secant method
   double max_abs         = 100.0;                    // Maximum value of |f(x)|.
   double bandrescale     = 1.0;                      // Rescale the input data by this scale factor
   bool adapt; // If adapt=false --> g(x)=1.
   bool hardgap;
   double boundary;
   double intA; // intA=int_0^1 rho(w) dw (trapezoidal method).
   double A;    // parameter in the shooting method. Initially A=intA. Equal to intA if adapt=false.
   // Right-hand-side of the differential equation. y=g !
   auto rhs_G(const double x, const double y) {
     const auto powL = Lambda.power(2.0 - x);
     return Lambda.logL() * (y - A / rho(y * powL));
   }
   void save(std::ostream &OUT) { OUT << x << " " << y << std::endl; }
   template<typename FNC> auto rk_step(const double dx, FNC rhs)
   {
     const auto k1 = dx * rhs(x, y);
     const auto k2 = dx * rhs(x + dx / 2.0, y + k1 / 2.0);
     const auto k3 = dx * rhs(x + dx / 2.0, y + k2 / 2.0);
     const auto k4 = dx * rhs(x + dx, y + k3);
     return (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
   }
   template<typename FNC> auto heun_step(const double dx, FNC rhs)
   {
     const auto k1 = dx * rhs(x, y);
     const auto k2 = dx * rhs(x + dx, y + k1);
     return (k1 + k2) / 2.0;
   }
   // Perform a step of length (at most) try_dx. Returns actual dx used.
   template<typename FNC> double timestep(const double try_dx, FNC rhs, bool check_f) {
     double dx = try_dx; // Current step length
     double dy, error;
     int subdivision = 0;
     // Approach: we attempt Runge-Kutta steps, reducing dx if necessary so as to avoid crossing the tabulation data
     // interval boundaries (with too large steps) and to avoid too large errors.
     int cnt           = 0;
     const int max_cnt = 2 * max_subdiv;
     while (true) {
       assert(std::isfinite(dx));
       rho.clear_flag();
       dy = rk_step(dx, rhs);
       // Have we crossed the interval boundary ?
       if (rho.flag()) {
         if (subdivision < max_subdiv) {
           subdivision++;
           dx /= 2.0;
           continue;
         }
       }
       // Calculate a second-order result using Heun's formula and estimate the error in y.
       const auto dy_heun = heun_step(dx, rhs);
       error          = abs(dy_heun - dy);
       max_error      = std::max(max_error, error); // update max_error
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
       if (cnt > max_cnt) 
         throw std::runtime_error("Number of subdivisions exceeded.");
       break; // exit loop here: the RK step is appropriate
     }
     if (check_f) {
       // Sanity checks
       const double dEdx      = dy/dx - y*Lambda.logL();
       const double tolerance = 1e-5;
       if (dEdx >= 0.0) { // E(x) must be monotonously decreasing!
         std::cerr << "WARNING: dE/dx is not negative." << std::endl;
         std::cerr << " x=" << x << " y=" << y << " dEdx=" << dEdx << " dy/dx=" << dy/dx << " y*log(Lambda)=" << y*Lambda.logL() << std::endl;
       }
       if (dEdx > tolerance)
         throw std::runtime_error("Tolerance criterium not satisfied.");
     }
     x += dx;
     y += dy;
     return dx;
   }
   // Integrate with step dx up to xf. max_diff_x = maximum discrepancy in x when approaching mesh points.
   template<typename FNC>
   void int_with_to(const double dx, const double xf, FNC rhs, const bool check_f, const double max_diff_x = 1e-10, const int max_cnt = 1000000) {
     int cnt = 0;
     do {
       const double steplength = x+dx < xf ? dx : xf-x; // If dx is too long, steplength should be decreased to xf-x.
       timestep(steplength, rhs, check_f);
       cnt++;
       if (cnt > max_cnt) 
         std::runtime_error("Number of subdivisions exceeded. max_subdiv should be increased.");
     } while (abs(x-xf) > max_diff_x); // x is a global variable, updated in timestep()
   }

   // returns the g(xmax)/[A/rho(0)] ratio and vecg
   auto shoot_g() {
     const std::string filename_g = "GSOL" + std::string(sign == Sign::POS ? "" : "NEG") + ".dat";
     std::ofstream OUTG;
     safe_open(OUTG, filename_g);
     std::cout << "#  A=" << A << std::endl;
     Vec vecg;
     rho(1); // NECESSARY!
     double dx = dx_fine; // initial step size
     // Initial conditions
     x         = 2.0;
     y         = 1.0; // y=g here!
     max_error = 0.0;
     save(OUTG); // save x=2 data point
     vecg.emplace_back(std::make_pair(x, y));
     double x_st = x; // Target x for next output line
     do {
       x_st += output_step;
       int_with_to(dx, x_st, [this](const auto x, const auto y) { return rhs_G(x, y); }, false); // rhs_G !!
       save(OUTG);
       vecg.emplace_back(std::make_pair(x, y));
       if (x > xfine) { dx = dx_fast; }
     } while (x < xmax);
     const auto factor = A / rho(0);
     const auto ratio  = y / factor;
     std::cout << "#  x_last=" << x << " " << "g_last/factor=" << ratio << std::endl;
     std::cout << "#  eps_last=" << y * Lambda.power(2 - x) << " max_error=" << max_error << std::endl;
     return std::make_pair(ratio, vecg);
   }
   void init_A() {
     intA = integrate_ab(vecrho, 0.0, 1.0);
     std::cout << "# intA=" << intA << std::endl;
     A = intA;
   }
   Vec calc_g(const Vec &vecrho) {
     const auto [ratio1, vecg1] = shoot_g();
     if (abs(ratio1 - 1.0) < convergence_eps)  // We're done
       return vecg1;
     // Otherwise more effort is required...
     A = ratio1 > 1.0 ? A * factor0 : A / factor0;
     const auto [ratio2, vecg2] = shoot_g();
     if (abs(ratio2 - 1.0) < convergence_eps)  // We're done
       return vecg2;
     // Refine using secant method
     double x0 = intA;
     double x1 = A;
     double y0 = ratio1 - 1.0;
     double y1 = ratio2 - 1.0;
     int iter  = 0;
     double ynew;
     do {
       const auto xnew            = x1 - (x1 - x0) / (y1 - y0) * y1;
       A                          = xnew;
       const auto [rationew,vecg] = shoot_g();
       ynew                       = rationew - 1.0;
       if (abs(rationew-1) < convergence_eps)
         return vecg;
       // Shift
       x0 = x1;
       y0 = y1;
       x1 = xnew;
       y1 = ynew;
       iter++;
     } while (iter < max_iter);
     throw std::runtime_error("Secant method failed to converge in " + std::to_string(max_iter) + " steps.");
   }
   // Rescaling of omega for excluding finite intervals around omega=0. The new accumulation point is determined by
   // the variable 'boundary'.
   auto rescale(const double omega) { return (1.0 - boundary) * omega + boundary; }
   // eps(x) = D g(x) Lambda^(2-x) for x>2.
   auto eps(const double x) {
     const auto gx = adapt ? g(x) : 1.0;
     double epsilon  = x <= 2.0 ? 1.0 : gx * Lambda.power(2.0-x);
     if (hardgap) { epsilon = rescale(epsilon); }
     return epsilon;
   }
   // Eps(x) = D f(x) Lambda^(2-x)
   inline auto Eps(const double x, const double f) {
     assert(x >= 1 && f > 0);
     return f * Lambda.power(2.0-x);
   }
   // Right-hand-side of the differential equation. y=f !
   auto rhs_F(const double x, const double y) {
     assert(std::isfinite(x));
     assert(std::isfinite(y));
     const double term1 = Lambda.logL() * y;
     const double integral = intrho2(eps(x)) - intrho1(eps(x + 1));
     const double powL     = Lambda.power(2.0 - x);
     const double denom    = powL * rho(y * powL);
     double term2 = integral / denom;
     if (denom == 0.0) {
       std::cout << "# Warning: denom=0 with integral=" << integral << " at x=" << x << std::endl;
       std::cout << "# (denom=0 may arise from zeros in the hybridisation function)" << std::endl;
       term2 = 0.0;
     }
     assert(std::isfinite(term2));
     return term1 - term2;
   }
   void load_init_rho() {
     std::string rhofn = P.Pstr("dos", "Delta.dat");
     vecrho       = load_rho(rhofn, sign);
     add_zero_point(vecrho);
     rescalevecxy(vecrho, 1.0/bandrescale, bandrescale);
     minmaxvec(vecrho, "rho");
     rho = LinInt(vecrho);
     std::cout << "# rho(0)=" << rho(0) << " rho(1)=" << rho(1) << std::endl;
     Vec vecintrho(vecrho);
     integrate(vecintrho);
     intrho1 = IntLinInt(vecrho, vecintrho);
     intrho2 = intrho1; // trick: two objects, one at upper, another at lower boundary. significant performance improvement!

   }
   void report_parameters() {
     std::cout << "# ++ " << (adapt ? "ADAPTIVE" : "FIXED-GRID") << std::endl;
     std::cout << "# Lambda=" << Lambda;
     std::cout << " bandrescale=" << bandrescale;
     std::cout << " xmax=" << xmax;
     std::cout << " xfine=" << xfine;
     std::cout << " output_step=" << output_step;
     std::cout << " dx_fine=" << dx_fine << " dx_fast=" << dx_fast;
     std::cout << std::endl;
     std::cout << "# allowed_error=" << allowed_error;
     std::cout << " max_subdiv=" << max_subdiv;
     std::cout << " max_abs=" << max_abs;
     std::cout << std::endl;
   }
   void set_parameters(const int PREC = 16) {
     std::cout << std::setprecision(PREC);
     Lambda = LAMBDA(P.P("Lambda", 2.0));
     assert(Lambda > 1.0);
     adapt = P.Pbool("adapt", false); // Enable adaptable g(x)? Default is false!!
     hardgap  = P.Pbool("hardgap", false); // Exclude an interval around omega=0 ?
     boundary = P.P("boundary", 0.0);      // The boundary of the exclusion interval.
     bandrescale = P.P("bandrescale", 1.0); // band rescaling parameter
     xmax = P.P("xmax", 30); // Integrate over [1..xmax]
     assert(xmax > 1);
     xfine = P.P("xfine", 5); // Fine stepsize integral [1..xfine]
     assert(xfine > 1);
     output_step = P.P("outputstep", 1.0 / 64.0); // Stepsize for output file
     assert(output_step <= 1.0);
     dx_fine       = P.P("dx_fine", 1e-5);        // Integration stepsize in [1..xfine]
     dx_fast       = P.P("dx_fast", 1e-4);        // Integration stepsize in [xfine..xmax]
     allowed_error = P.P("allowed_error", 1e-10); // error control for adaptable stepsize
     max_subdiv    = P.Pint("max_subdiv", 10);    // maximum nr of integ. step subdivisions
     assert(dx_fine * pow(0.5, max_subdiv) > DBL_EPSILON);
     max_abs = P.P("max_abs", 100.0); // Maximal |f(x)|
     assert(max_abs > 0.0);
     convergence_eps = P.P("secant_eps", 1e-4);
     factor0         = 1.0 + P.P("secant_factor", 1e-7);
     max_iter        = P.Pint("secant_max_iter", 10);
   }
   Adapt(const Params &P, const Sign &sign) : P(P), sign(sign) {
     set_parameters();
     report_parameters();
   }
   auto g_fn(const Sign &sign) { return "GSOL" + (sign == Sign::POS ? ""s : "NEG"s) + ".dat"; }
   auto f_fn(const Sign &sign) { return "FSOL" + (sign == Sign::POS ? ""s : "NEG"s) + ".dat"; }
   void load_or_calc_g() {
     const auto vecg = P.Pbool("loadg", false) ? load_g(g_fn(sign)) : calc_g(vecrho);
     minmaxvec(vecg, "g");
     g = LinInt(vecg);
   }
   void report() {
     const double factor = Lambda.factor() * A / rho(0);
     std::cout << "# x_last=" << x << std::endl;
     std::cout << "# f_last=" << y << " [f_last/factor=" << y / factor << "]" << std::endl;
     std::cout << "# eps_last=" << eps(x) << " (smallest energy point considered [input])" << std::endl;
     std::cout << "# Eps_last=" << y * Lambda.power(2 - x) << " (smallest energy scale obtained [output])" << std::endl;
     std::cout << "# max_error=" << max_error << " (maximum integration error)" << std::endl;
   }
   void calc_f() {
     std::ofstream OUTF;
     safe_open(OUTF, f_fn(sign));
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
       int_with_to(dx, x_st, [this](const auto x, const auto y){ return rhs_F(x, y); }, true); // rhs_F !!
       save(OUTF);
       if (x > xfine) { dx = dx_fast; }
       if (abs(y) > max_abs) {
         std::cout<< "***** y=" << y << " |y|>max_abs=" << max_abs << std::endl;
         std::cout<< "***** Terminating!" << std::endl;
       }
     } while (x < xmax && abs(y) <= max_abs);
   }
   void run() {
     load_init_rho();
     init_A();
     if (adapt) { load_or_calc_g(); }
     calc_f();
     report();
   }
};

} // namespace

#endif
