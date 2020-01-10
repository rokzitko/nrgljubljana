// Calculation of NRG chain coefficients
// Rok Zitko, rok.zitko@ijs.si, 2009-2020

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
#include <limits>

#include <gmp.h>

using namespace std;

#include "lambda.h"
#include "linint.h"
#include "io.h"
#include "parser.h"
#include "load.h"
#include "calc.h"

string param_fn = "param";  // file with input parameters
LAMBDA Lambda;              // discretization parameter
double z;                   // twist parameter
double xmax;                // higher boundary of the x=j+z interval, where the ODE
                            // was solved numerically
unsigned int mMAX;          // the number of coefficients computed (max index)
unsigned int Nmax;          // the length of the Wilson chain (max index)
double bandrescale = 1.0;   // band rescaling factor
bool rescalexi     = false; // rescale coefficients xi

unsigned int preccpp; // precision for GMP

Vec vecrho_pos, vecrho_neg; // rho, for positive and negative energies
LinInt rho_pos, rho_neg;
IntLinInt intrho_pos, intrho_neg; // integrated rho
LinInt g_pos, g_neg;              // g(x)
LinInt f_pos, f_neg;

using Table = vector<double>;

// Input to the tridiagonalisation.
Table de_pos, de_neg, du_pos, du_neg;

bool adapt; // If adapt=false --> g(x)=1.

string band; // If band="flat", we use an analytical expression for f,
             // otherwise we load "FSOL.dat" and "FSOLNEG.dat"

bool nrgchain_tables_save; // If nrg_tables_save=true, coefficient tables are
                           // written to files.
bool nrgchain_tables_load; // If nrg_tables_load=true, coefficient tables are
                           // read from files.
bool nrgchain_tridiag;     // If nrgchain_tridiag=true, tridiagonalisation is
                           // performed.

// eps(x) = D g(x) Lambda^(2-x) for x>2.
// This is only an auxiliary quantity which defines the discretization
// mesh.
double eps_pos(double x) {
  const double gx = (adapt ? g_pos(x) : 1.0);
  return (x <= 2.0 ? 1.0 : gx * Lambda.power(2.0 - x));
}

double eps_neg(double x) {
  const double gx = (adapt ? g_neg(x) : 1.0);
  return (x <= 2.0 ? 1.0 : gx * Lambda.power(2.0 - x));
}

// Analytical expression for Epsilon(x) in the case of a flat band.
// Cf. PRB 79, 085106 (2009), Eqs. (25) & (36).
inline double Eps_flat(double x) {
  assert(x >= 1.0);
  const int j    = floor(x);
  const double z = x - j;

  if (j == 1) {
    return (1.0 - Lambda.power(-z)) / Lambda.logL() + 1.0 - z;
  } else {
    return (1.0 - Lambda.power(-1.0)) / Lambda.logL() * Lambda.power(2.0 - j - z);
  }
}

// Eps(x) = D f(x) Lambda^(2-x)
// This are the "representative energies" of the grid.
inline double Eps_pos(double x) {
  if (band == "flat") return Eps_flat(x);

  assert(x >= 1.0);
  const double f = f_pos(x);
  return f * Lambda.power(2.0 - x);
}

inline double Eps_neg(double x) {
  if (band == "flat") return Eps_flat(x);

  assert(x >= 1.0);
  const double f = f_neg(x);
  return f * Lambda.power(2.0 - x);
}

void about(ostream &F = cout) {
  F << "# Calculation of NRG coefficients" << endl;
  F << "# Rok Zitko, rok.zitko@ijs.si, 2009" << endl;
}

// Called before the parser.
void cmd_line(int argc, char *argv[]) {}

// Called after the parser.
void cmd_line_post(int argc, char *argv[]) {
  // Argument 's': save tables, do not tridiagonalise
  if (argc == 2 && string(argv[1]) == "s") {
    nrgchain_tables_load = false;
    nrgchain_tables_save = true;
    nrgchain_tridiag     = false;
  }

  // Argument 'l': load tables, tridiagonalise
  if (argc == 2 && string(argv[1]) == "l") {
    nrgchain_tables_load = true;
    nrgchain_tables_save = false;
    nrgchain_tridiag     = true;
  }
}

void set_parameters() {
  cout << setprecision(PREC);

  Lambda = LAMBDA(P("Lambda", 2.0));
  assert(Lambda > 1.0);

  z = P("z", 1.0);
  assert(0 < z && z <= 1.0);

  adapt = Pbool("adapt", false); // Enable adaptable g(x)? Default is false!!

  bandrescale = P("bandrescale", 1.0);
  rescalexi   = Pbool("rescalexi", false);

  xmax = P("xmax", 30); // Interval [1..xmax]
  assert(xmax >= 1.0);

  Nmax = Pint("Nmax", 0); // Maximal site index in the Wilson chain

  mMAX = Pint("mMAX", 2 * Nmax); // Maximal index of coefficients (e,f)

  assert(mMAX > 0);

  preccpp = Pint("preccpp", 2000); // Precision for GMP
  assert(preccpp > 10);

  band = Pstr("band", "adapt"); // Default: load FSOL*.dat

  nrgchain_tables_save = Pbool("nrgchain_tables_save", false);
  nrgchain_tables_load = Pbool("nrgchain_tables_load", false);
  nrgchain_tridiag     = Pbool("nrgchains_tridiag", true);

  cout << "# Lambda=" << Lambda;
  cout << " bandrescale=" << bandrescale;
  cout << " z=" << z << endl;
  cout << "# xmax=" << xmax;
  cout << " mMAX=" << mMAX;
  cout << " Nmax=" << Nmax << endl;
  cout << "# band=" << band << endl;
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

void load_rho() {
  const string rhofn = Pstr("dos", "Delta.dat");

  vecrho_pos         = load_rho(rhofn, POS);
  rescalevecxy(vecrho_pos, 1.0 / bandrescale, bandrescale);
  add_zero_point(vecrho_pos);

  vecrho_neg = load_rho(rhofn, NEG);
  rescalevecxy(vecrho_neg, 1.0 / bandrescale, bandrescale);
  add_zero_point(vecrho_neg);
}

void init_rho() {
  rho_pos = LinInt(vecrho_pos);
  rho_neg = LinInt(vecrho_neg);

  Vec vecintrho_pos(vecrho_pos);
  integrate(vecintrho_pos);
  intrho_pos = IntLinInt(vecrho_pos, vecintrho_pos);

  Vec vecintrho_neg(vecrho_neg);
  integrate(vecintrho_neg);
  intrho_neg = IntLinInt(vecrho_neg, vecintrho_neg);
}

void load_g() {
  const string gfn_pos = "GSOL.dat";
  Vec vecg_pos         = load_g(gfn_pos);
  g_pos                = LinInt(vecg_pos);

  const string gfn_neg = "GSOLNEG.dat";
  Vec vecg_neg         = load_g(gfn_neg);
  g_neg                = LinInt(vecg_neg);
}

void load_f() {
  const string ffn_pos = "FSOL.dat";
  Vec vecf_pos         = load_g(ffn_pos); // same load_g() function as for g
  f_pos                = LinInt(vecf_pos);

  const string ffn_neg = "FSOLNEG.dat";
  Vec vecf_neg         = load_g(ffn_neg);
  f_neg                = LinInt(vecf_neg);
}

// The factor that multiplies eigenvalues of the Wilson chain Hamiltonian
// in order to obtain the eigenvalues of the true Hamiltonian (at scale D).
double SCALE(int N) { return (1.0 - 1. / Lambda) / log(Lambda) * pow(Lambda, -(N - 1.0) / 2.0 + 1.0 - z); }

inline double sqr(double x) { return x * x; }

void tables() {
  const double int_pos1 = integrate_ab(vecrho_pos, 0.0, 1.0);
  const double int_neg1 = integrate_ab(vecrho_neg, 0.0, 1.0);
  const double theta1   = int_pos1 + int_neg1;
  cout << "# int_pos1=" << int_pos1 << " int_neg1=" << int_neg1 << " theta1=" << theta1 << endl;
  const double int_pos2 = intrho_pos(eps_pos(z+1)) - intrho_pos(eps_pos(z + mMAX + 2));
  const double int_neg2 = intrho_neg(eps_neg(z+1)) - intrho_neg(eps_neg(z + mMAX + 2));
  const double theta2 = int_pos2 + int_neg2;
  cout << "# int_pos2=" << int_pos2 << " int_neg2=" << int_neg2 << " theta2=" << theta2 << endl;

  // For consistency with df_pos & df_neg, we use set 2
  const double int_pos = int_pos2;
  const double int_neg = int_neg2;
  const double theta = theta2;
  
  ofstream THETA;
  safe_open(THETA, "theta.dat"); // theta (hybridisation fnc. weight)
  THETA << theta << endl;
  THETA.close();

  Table df_pos(mMAX + 1), df_neg(mMAX + 1);
  Table du0_neg(mMAX + 1), du0_pos(mMAX + 1);

  de_pos.resize(mMAX + 1);
  de_neg.resize(mMAX + 1);

  for (unsigned int m = 0; m <= mMAX; m++) {
    df_pos[m] = intrho_pos(eps_pos(z + m + 1)) - intrho_pos(eps_pos(z + m + 2));
    df_neg[m] = intrho_neg(eps_neg(z + m + 1)) - intrho_neg(eps_neg(z + m + 2));

    du0_pos[m] = sqrt(df_pos[m]) / sqrt(theta);
    du0_neg[m] = sqrt(df_neg[m]) / sqrt(theta);

    de_pos[m] = Eps_pos(z + m + 1);
    de_neg[m] = Eps_neg(z + m + 1);
  }

  double checksum = 0.0;
  for (unsigned int m = 0; m <= mMAX; m++) checksum += sqr(du0_pos[m]) + sqr(du0_neg[m]);

  cout << "# 1-checksum=" << 1 - checksum << endl;

  // A large deviation probably indicates a serious problem!
  const double CHECKSUM_LIMIT = 1e-10;
  if (abs(1 - checksum) > CHECKSUM_LIMIT) {
    cerr << "Checksum test failed." << endl;
    exit(1);
  }

  du_pos.resize(mMAX + 1);
  du_neg.resize(mMAX + 1);

  for (unsigned int m = 0; m <= mMAX; m++) {
    du_pos[m] = du0_pos[m] / sqrt(checksum);
    du_neg[m] = du0_neg[m] / sqrt(checksum);
  }

  for (unsigned int m = 0; m <= mMAX; m++) {
    cout << "# " << m << " " << du_pos[m] << " " << du_neg[m] << " " << de_pos[m] << " " << de_neg[m] << endl;
  }
}

void save_tables() {
  save("de_pos.dat", de_pos);
  save("de_neg.dat", de_neg);
  save("du_pos.dat", du_pos);
  save("du_neg.dat", du_neg);
}

void load_tables() {
  load("de_pos.dat", de_pos);
  load("de_neg.dat", de_neg);
  load("du_pos.dat", du_pos);
  load("du_neg.dat", du_neg);
}

class my_mpf {
  private:
  mpf_t val{};

  public:
  my_mpf() { mpf_init(val); }
  // Copy constructor is mendatory!
  my_mpf(const my_mpf &x) {
    mpf_init(val);
    mpf_set(val, x.val);
  }
  ~my_mpf() { mpf_clear(val); }
  inline operator mpf_t &() { return val; }
};

using vmpf = std::vector<my_mpf>;

// Fix normalization of u_{n,m}, v_{n,m} to 1. IMPORTANT: pass by
// reference!
void fix_norm(vmpf &up, vmpf &um, unsigned int mMAX) {
  // Constants
  my_mpf mpZERO, mpONE;
  mpf_set_str(mpONE, "1.e0", 10);

  my_mpf sum, temp, tempsq;

  mpf_set(sum, mpZERO);
  for (unsigned int m = 0; m <= mMAX; m++) {
    mpf_mul(tempsq, up[m], up[m]);
    mpf_add(sum, sum, tempsq);
    mpf_mul(tempsq, um[m], um[m]);
    mpf_add(sum, sum, tempsq);
  }
  mpf_sqrt(temp, sum);

  for (unsigned int m = 0; m <= mMAX; m++) {
    mpf_div(up[m], up[m], temp);
    mpf_div(um[m], um[m], temp);
  }
}

#define HIGHPREC(val) setw(30) << setprecision(16) << (val) << setprecision(PREC)

// Triagonalisation by iteration.
//
// INPUT: tables du_pos, du_neg, de_pos, de_neg
// OUTPUT: written to files "xi.dat" and "zeta.dat"

void tridiag() {
  ofstream XI, ZETA;
  safe_open(XI, "xi.dat");     // hopping constants
  safe_open(ZETA, "zeta.dat"); // on-site energies

  mpf_set_default_prec(preccpp);
  cout << "Using precision of " << preccpp << " digits." << endl;

  // Constants
  my_mpf mpZERO;

  // Temporary MP variables
  my_mpf temp, tempsq, sum;

  my_mpf mpxi;   // xi
  my_mpf xi2;    // xi^2
  my_mpf mpzeta; // zeta

  my_mpf xi_prev, xi2_prev; // values in previous iteration

  vmpf up(mMAX + 1);
  vmpf up_prev(mMAX + 1);
  vmpf up_prev2(mMAX + 1);
  vmpf um(mMAX + 1);
  vmpf um_prev(mMAX + 1);
  vmpf um_prev2(mMAX + 1);
  vmpf ep1(mMAX + 1);
  vmpf em1(mMAX + 1);
  vmpf ep2(mMAX + 1);
  vmpf em2(mMAX + 1);
  for (unsigned int m = 0; m <= mMAX; m++) {
    mpf_set_d(up_prev[m], du_pos[m]);
    mpf_set_d(um_prev[m], du_neg[m]);
    mpf_set_d(ep1[m], de_pos[m]);
    mpf_set_d(em1[m], de_neg[m]);
    mpf_mul(ep2[m], ep1[m], ep1[m]);
    mpf_mul(em2[m], em1[m], em1[m]);
  }

  fix_norm(up_prev, um_prev, mMAX);

  for (unsigned int n = 0; n <= Nmax; n++) {
    // Calculate zeta_n, xi2_n and xi_n
    mpf_set(mpzeta, mpZERO);
    mpf_set(xi2, mpZERO);
    for (unsigned int m = 0; m <= mMAX; m++) {
      // up_prev = u^+_{n,m}
      mpf_mul(tempsq, up_prev[m], up_prev[m]);
      mpf_mul(temp, tempsq, ep2[m]);
      mpf_add(xi2, xi2, temp);
      mpf_mul(temp, tempsq, ep1[m]);
      mpf_add(mpzeta, mpzeta, temp);

      // um_prev = u^-_{n,m}
      mpf_mul(tempsq, um_prev[m], um_prev[m]);
      mpf_mul(temp, tempsq, em2[m]);
      mpf_add(xi2, xi2, temp);
      mpf_mul(temp, tempsq, em1[m]);
      mpf_sub(mpzeta, mpzeta, temp);
    }

    // subtract xi^2_{n-1}
    mpf_sub(xi2, xi2, xi2_prev);

    // subtract zeta^2_n
    mpf_mul(tempsq, mpzeta, mpzeta);
    mpf_sub(xi2, xi2, tempsq);

    if (!mpf_cmp_d(xi2, 0.0)) {
      cerr << "xi2 negative, aborting." << endl;
      exit(1);
    }

    mpf_sqrt(mpxi, xi2);

    // compute u_{n+1,m}, v_{n+1,m}
    for (unsigned int m = 0; m <= mMAX; m++) {
      // zeta=zeta_n
      mpf_sub(temp, ep1[m], mpzeta);
      // up_prev[m]=u_{n,m}
      mpf_mul(up[m], temp, up_prev[m]);

      // xi_prev=xi_{n-1}, up_prev2=u_{n-1,m}
      mpf_mul(temp, xi_prev, up_prev2[m]);
      mpf_sub(up[m], up[m], temp);

      // xi=xi_n
      mpf_div(up[m], up[m], mpxi);

      mpf_neg(temp, em1[m]);
      mpf_sub(temp, temp, mpzeta);
      mpf_mul(um[m], temp, um_prev[m]);
      mpf_mul(temp, xi_prev, um_prev2[m]);
      mpf_sub(um[m], um[m], temp);
      mpf_div(um[m], um[m], mpxi);
    }

    fix_norm(up, um, mMAX);

    // Recalculate xi, xi2
    mpf_set(sum, mpZERO);
    for (unsigned int m = 0; m <= mMAX; m++) {
      mpf_mul(temp, up[m], up_prev[m]);
      mpf_mul(temp, temp, ep1[m]);
      mpf_add(sum, sum, temp);
      mpf_mul(temp, um[m], um_prev[m]);
      mpf_mul(temp, temp, em1[m]);
      mpf_sub(sum, sum, temp);
    }
    mpf_set(mpxi, sum);
    mpf_mul(xi2, mpxi, mpxi);

    // Save results
    double dxi       = mpf_get_d(mpxi);
    double dzeta     = mpf_get_d(mpzeta);
    double coef_xi   = dxi / (rescalexi == true ? SCALE(n + 1) : 1.0);
    double coef_zeta = dzeta; // NEVER RESCALED!!!
    
    coef_xi = coef_xi * bandrescale; // by analogy with initial.m
    coef_zeta = coef_zeta * bandrescale;

    XI << coef_xi << endl;
    ZETA << coef_zeta << endl;

    cout << "  xi(" << n << ")=" << HIGHPREC(dxi) << " --> " << HIGHPREC(coef_xi) << endl;
    cout << "zeta(" << n << ")=" << HIGHPREC(dzeta) << endl;

    // Store results from previous iteration
    mpf_set(xi_prev, mpxi);
    mpf_set(xi2_prev, xi2);

    for (unsigned int m = 0; m <= mMAX; m++) {
      mpf_set(um_prev2[m], um_prev[m]);
      mpf_set(up_prev2[m], up_prev[m]);
      mpf_set(um_prev[m], um[m]);
      mpf_set(up_prev[m], up[m]);
    }
  }
}

void calc_tables() {
  if (band != "flat") {
    load_rho();
    if (adapt) { load_g(); }
    load_f();
  } else {
    // Flat band
    const double mindbl = numeric_limits<double>::min();
    Vec v;
    v.push_back(make_pair(mindbl, 0.5));
    v.push_back(make_pair(1.0, 0.5));
    vecrho_pos = v;
    vecrho_neg = v;
  }

  init_rho();

  tables();

  if (nrgchain_tables_save) { save_tables(); }
}

int main(int argc, char *argv[]) {
  clock_t start_clock = clock();

  about();
  cmd_line(argc, argv);
  parser(param_fn);
  set_parameters();
  cmd_line_post(argc, argv);

  if (nrgchain_tables_load) {
    load_tables();
  } else {
    calc_tables();
  }

  if (nrgchain_tridiag) tridiag();

  clock_t end_clock = clock();
  cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << endl;
}
