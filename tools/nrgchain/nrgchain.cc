// Calculation of NRG chain coefficients
// Rok Zitko, rok.zitko@ijs.si

#include <exception>
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
#include <iterator>
#include <ostream>
#include <stdexcept>

#include <gmp.h>

#include "../common/tabulated_density.hpp"

#ifndef NRGCHAIN_NO_MAIN
#include <common/version.hpp>

#include "../common/diagnostics.hpp"
#endif

using namespace std;

#include "lambda.h"
#include "linint.h"
#include "io.h"
#include "parser.h"
#include "load.h"
#include "nrgchain.hpp"

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
NRG::Tools::TabulatedDensity rho_pos, rho_neg;
LinInt g_pos, g_neg;              // g(x)
LinInt f_pos, f_neg;

using Table = vector<double>;

Table result_xi, result_zeta;
double result_theta = 0.0;

// Input to the tridiagonalisation.
Table de_pos, de_neg, du_pos, du_neg;

bool adapt; // If adapt=false --> g(x)=1.
NRG::Tools::InterpolationMethod density_interpolation = NRG::Tools::InterpolationMethod::linear;

string band; // If band="flat", we use an analytical expression for f,
             // otherwise we load "FSOL.dat" and "FSOLNEG.dat"

bool nrgchain_tables_save; // If nrg_tables_save=true, coefficient tables are
                           // written to files.
bool nrgchain_tables_load; // If nrg_tables_load=true, coefficient tables are
                           // read from files.
bool nrgchain_tridiag;     // If nrgchain_tridiag=true, tridiagonalisation is
                           // performed.
filesystem::path nrgchain_output_dir = ".";

string output_path(const string &filename) { return (nrgchain_output_dir / filename).string(); }

void close_output_checked(ofstream &output, const string &filename) {
  output.close();
  if (!output) throw runtime_error("Error writing " + filename + ".");
}

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
  const double z_ = x - j;

  if (j == 1) {
    return (1.0 - Lambda.power(-z_)) / Lambda.logL() + 1.0 - z_;
  } else {
    return (1.0 - Lambda.power(-1.0)) / Lambda.logL() * Lambda.power(2.0 - j - z_);
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
}

constexpr auto usage_text = "Usage: nrgchain [options] [s|l] [parameter_file]";

void usage(ostream &F = cout) {
  F << usage_text << endl;
  F << " -h, --help -- show help" << endl;
  F << " -v -- show resolved configuration and verbose diagnostics on stderr" << endl;
  F << " -vv -- increase verbosity further" << endl;
  F << " -V, --version -- show project version and exit" << endl;
  F << " s -- calculate and save tables without tridiagonalising" << endl;
  F << " l -- load tables and tridiagonalise" << endl;
}

struct CommandLineOptions {
  string param_filename = "param";
  NRG::Tools::NrgChain::TableMode mode = NRG::Tools::NrgChain::TableMode::Calculate;
  int verbosity = 0;
};

#ifndef NRGCHAIN_NO_MAIN
int command_line_verbosity = 0;
string command_line_param_filename = "param";
#endif

CommandLineOptions cmd_line(int argc, char *argv[]) {
  CommandLineOptions options;
  bool mode_set  = false;
  bool param_set = false;

  for (int i = 1; i < argc; i++) {
    const string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage();
      exit(EXIT_SUCCESS);
    }
    if (arg.size() >= 2 && arg[0] == '-' && all_of(arg.begin() + 1, arg.end(), [](const char ch) { return ch == 'v'; })) {
      options.verbosity += static_cast<int>(arg.size() - 1);
      continue;
    }
    if (arg == "s" || arg == "l") {
      if (mode_set) throw invalid_argument("Mode specified more than once.\n" + string(usage_text));
      options.mode = arg == "s" ? NRG::Tools::NrgChain::TableMode::SaveOnly
                                : NRG::Tools::NrgChain::TableMode::LoadAndTridiagonalize;
      mode_set = true;
      continue;
    }
    if (!arg.empty() && arg[0] == '-') throw invalid_argument("Unknown option: " + arg + "\n" + usage_text);
    if (param_set) throw invalid_argument("Unexpected argument: " + arg + "\n" + usage_text);
    options.param_filename = arg;
    param_set              = true;
  }

  return options;
}

void set_parameters() {
  cout << setprecision(PREC);

  Lambda = LAMBDA(P("Lambda", 2.0));
  if (!(Lambda > 1.0)) throw std::invalid_argument("Lambda must be greater than 1.");

  z = P("z", 1.0);
  if (!(0 < z && z <= 1.0)) throw std::invalid_argument("z must satisfy 0 < z <= 1.");

  adapt = Pbool("adapt", false); // Enable adaptable g(x)? Default is false!!

  bandrescale = P("bandrescale", 1.0);
  if (!(std::isfinite(bandrescale) && bandrescale > 0.0))
    throw std::invalid_argument("bandrescale must be positive and finite.");
  density_interpolation = NRG::Tools::parse_density_interpolation_method(
    Pstr("density_interpolation", "linear"));
  rescalexi   = Pbool("rescalexi", false);

  xmax = P("xmax", 30); // Interval [1..xmax]
  if (!(std::isfinite(xmax) && xmax >= 1.0))
    throw std::invalid_argument("xmax must be finite and greater than or equal to 1.");

  const auto nmax_value = Pint("Nmax", 0); // Maximal site index in the Wilson chain
  if (nmax_value < 0) throw std::invalid_argument("Nmax must be nonnegative.");
  Nmax = static_cast<unsigned int>(nmax_value);

  if (Nmax > static_cast<unsigned int>(std::numeric_limits<int>::max() / 2))
    throw std::invalid_argument("Nmax is too large to derive mMAX.");
  const auto mmax_value = Pint("mMAX", static_cast<int>(2 * Nmax)); // Maximal index of coefficients (e,f)
  if (mmax_value <= 0) throw std::invalid_argument("mMAX must be greater than 0.");
  mMAX = static_cast<unsigned int>(mmax_value);

  const auto precision_value = Pint("preccpp", 2000); // Precision for GMP
  if (precision_value <= 10) throw std::invalid_argument("preccpp must be greater than 10.");
  preccpp = static_cast<unsigned int>(precision_value);

  band = Pstr("band", "adapt"); // Default: load FSOL*.dat

  nrgchain_tables_save = Pbool("nrgchain_tables_save", false);
  nrgchain_tables_load = Pbool("nrgchain_tables_load", false);
  // Keep the historical misspelling as a fallback for existing parameter files.
  nrgchain_tridiag = params.contains("nrgchain_tridiag") ? Pbool("nrgchain_tridiag", true)
                                                          : Pbool("nrgchains_tridiag", true);

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
  rho_pos = NRG::Tools::TabulatedDensity(vecrho_pos, density_interpolation);
  rho_neg = NRG::Tools::TabulatedDensity(vecrho_neg, density_interpolation);
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
  const double int_pos1 = rho_pos.integral(0.0, 1.0);
  const double int_neg1 = rho_neg.integral(0.0, 1.0);
  const double theta1   = int_pos1 + int_neg1;
  cout << "# int_pos1=" << int_pos1 << " int_neg1=" << int_neg1 << " theta1=" << theta1 << endl;
  const double int_pos2 = rho_pos.integral(eps_pos(z + mMAX + 2), eps_pos(z + 1));
  const double int_neg2 = rho_neg.integral(eps_neg(z + mMAX + 2), eps_neg(z + 1));
  const double theta2 = int_pos2 + int_neg2;
  cout << "# int_pos2=" << int_pos2 << " int_neg2=" << int_neg2 << " theta2=" << theta2 << endl;

  // For consistency with df_pos & df_neg, we use set 2
  const double theta = theta2;
  if (!(std::isfinite(theta) && theta > 0.0))
    throw runtime_error("Hybridisation weight theta must be positive and finite.");
  result_theta = theta;
  
  ofstream THETA;
  const auto theta_filename = output_path("theta.dat");
  safe_open(THETA, theta_filename); // theta (hybridisation fnc. weight)
  THETA << setprecision(18) << theta << endl;
  close_output_checked(THETA, theta_filename);

  Table df_pos(mMAX + 1), df_neg(mMAX + 1);
  Table du0_neg(mMAX + 1), du0_pos(mMAX + 1);

  de_pos.resize(mMAX + 1);
  de_neg.resize(mMAX + 1);

  for (unsigned int m = 0; m <= mMAX; m++) {
    df_pos[m] = rho_pos.integral(eps_pos(z + m + 2), eps_pos(z + m + 1));
    df_neg[m] = rho_neg.integral(eps_neg(z + m + 2), eps_neg(z + m + 1));

    du0_pos[m] = sqrt(df_pos[m]) / sqrt(theta);
    du0_neg[m] = sqrt(df_neg[m]) / sqrt(theta);

    de_pos[m] = Eps_pos(z + m + 1);
    de_neg[m] = Eps_neg(z + m + 1);
    if (!(std::isfinite(de_pos[m]) && de_pos[m] > 0.0 && std::isfinite(de_neg[m]) && de_neg[m] > 0.0))
      throw runtime_error("Representative energies must be positive and finite.");
  }

  double checksum = 0.0;
  for (unsigned int m = 0; m <= mMAX; m++) checksum += sqr(du0_pos[m]) + sqr(du0_neg[m]);
  if (!(std::isfinite(checksum) && checksum > 0.0))
    throw runtime_error("Wilson-state normalization is not positive and finite.");

  cout << "# 1-checksum=" << 1 - checksum << endl;

  // A large deviation probably indicates a serious problem!
  const double CHECKSUM_LIMIT = 1e-10;
  if (std::abs(1 - checksum) > CHECKSUM_LIMIT) {
    throw runtime_error("Checksum test failed.");
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
  save(output_path("de_pos.dat"), de_pos);
  save(output_path("de_neg.dat"), de_neg);
  save(output_path("du_pos.dat"), du_pos);
  save(output_path("du_neg.dat"), du_neg);
}

void load_tables() {
  load("de_pos.dat", de_pos);
  load("de_neg.dat", de_neg);
  load("du_pos.dat", du_pos);
  load("du_neg.dat", du_neg);

  Table theta_values;
  load("theta.dat", theta_values);
  if (theta_values.size() != 1 || !(std::isfinite(theta_values.front()) && theta_values.front() > 0.0))
    throw runtime_error("theta.dat must contain one positive finite value.");
  result_theta = theta_values.front();

  const auto expected_size = static_cast<std::size_t>(mMAX) + 1;
  if (de_pos.size() != expected_size || de_neg.size() != expected_size || du_pos.size() != expected_size
      || du_neg.size() != expected_size)
    throw runtime_error("Loaded Wilson tables must each contain mMAX+1 values.");

  double checksum = 0.0;
  for (std::size_t index = 0; index < expected_size; ++index) {
    if (!(std::isfinite(de_pos[index]) && de_pos[index] > 0.0 && std::isfinite(de_neg[index])
          && de_neg[index] > 0.0))
      throw runtime_error("Loaded representative energies must be positive and finite.");
    if (!(std::isfinite(du_pos[index]) && du_pos[index] >= 0.0 && std::isfinite(du_neg[index])
          && du_neg[index] >= 0.0))
      throw runtime_error("Loaded Wilson amplitudes must be nonnegative and finite.");
    checksum += sqr(du_pos[index]) + sqr(du_neg[index]);
  }
  constexpr double checksum_limit = 1e-10;
  if (!std::isfinite(checksum) || std::abs(1.0 - checksum) > checksum_limit)
    throw runtime_error("Loaded Wilson-table normalization failed.");
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
void fix_norm(vmpf &up, vmpf &um, unsigned int mMAX_) {
  // Constants
  my_mpf mpZERO, mpONE;
  mpf_set_str(mpONE, "1.e0", 10);

  my_mpf sum, temp, tempsq;

  mpf_set(sum, mpZERO);
  for (unsigned int m = 0; m <= mMAX_; m++) {
    mpf_mul(tempsq, up[m], up[m]);
    mpf_add(sum, sum, tempsq);
    mpf_mul(tempsq, um[m], um[m]);
    mpf_add(sum, sum, tempsq);
  }
  mpf_sqrt(temp, sum);

  for (unsigned int m = 0; m <= mMAX_; m++) {
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
  const auto xi_filename = output_path("xi.dat");
  const auto zeta_filename = output_path("zeta.dat");
  safe_open(XI, xi_filename);     // hopping constants
  safe_open(ZETA, zeta_filename); // on-site energies
  result_xi.clear();
  result_zeta.clear();

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

    if (mpf_cmp_d(xi2, 0.0) <= 0) {
      throw runtime_error("xi2 non-positive, aborting.");
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
    result_xi.push_back(coef_xi);
    result_zeta.push_back(coef_zeta);

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
  close_output_checked(XI, xi_filename);
  close_output_checked(ZETA, zeta_filename);
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
    if (density_interpolation == NRG::Tools::InterpolationMethod::steffen)
      v.push_back(make_pair(0.5, 0.5));
    v.push_back(make_pair(1.0, 0.5));
    vecrho_pos = v;
    vecrho_neg = v;
  }

  init_rho();

  tables();

  if (nrgchain_tables_save) { save_tables(); }
}

void reset_calculation_state() {
  params.clear();
  Lambda = LAMBDA();
  z = 0.0;
  xmax = 0.0;
  mMAX = 0;
  Nmax = 0;
  bandrescale = 1.0;
  rescalexi = false;
  preccpp = 0;
  vecrho_pos.clear();
  vecrho_neg.clear();
  rho_pos = NRG::Tools::TabulatedDensity();
  rho_neg = NRG::Tools::TabulatedDensity();
  g_pos = LinInt();
  g_neg = LinInt();
  f_pos = LinInt();
  f_neg = LinInt();
  de_pos.clear();
  de_neg.clear();
  du_pos.clear();
  du_neg.clear();
  adapt = false;
  density_interpolation = NRG::Tools::InterpolationMethod::linear;
  band.clear();
  nrgchain_tables_save = false;
  nrgchain_tables_load = false;
  nrgchain_tridiag = true;
  result_xi.clear();
  result_zeta.clear();
  result_theta = 0.0;
  nrgchain_output_dir = ".";
}

namespace NRG::Tools::NrgChain {

namespace {

void apply_mode(const TableMode mode) {
  if (mode == TableMode::SaveOnly) {
    nrgchain_tables_load = false;
    nrgchain_tables_save = true;
    nrgchain_tridiag = false;
  }
  if (mode == TableMode::LoadAndTridiagonalize) {
    nrgchain_tables_load = true;
    nrgchain_tables_save = false;
    nrgchain_tridiag = true;
  }
}

#ifndef NRGCHAIN_NO_MAIN
const char *mode_name(const TableMode mode) {
  switch (mode) {
    case TableMode::Calculate: return "parameter-file";
    case TableMode::SaveOnly: return "save-only";
    case TableMode::LoadAndTridiagonalize: return "load-and-tridiagonalize";
  }
  return "unknown";
}

void report_configuration(const TableMode mode) {
  if (command_line_verbosity == 0) return;
  NRG::Tools::ConfigurationReport report("nrgchain");
  report.value("verbosity", command_line_verbosity);
  report.value("parameter_file", command_line_param_filename);
  report.value("output_directory", nrgchain_output_dir);
  report.value("mode", mode_name(mode));
  report.value("Lambda", Lambda);
  report.value("z", z);
  report.value("xmax", xmax);
  report.value("Nmax", Nmax);
  if (params.contains("mMAX"))
    report.value("mMAX", mMAX);
  else
    report.resolved("mMAX", mMAX, "2*Nmax");
  report.value("band", band);
  report.value("bandrescale", bandrescale);
  report.value("adapt", adapt);
  report.value("rescalexi", rescalexi);
  report.value("gmp_precision", preccpp);
  report.value("output_precision", PREC);
  const auto density_method = NRG::Tools::interpolation_method_name(density_interpolation);
  if (nrgchain_tables_load) {
    report.resolved("density_source", "inactive", "tables_load=true");
    report.resolved("density_interpolation", "inactive", "tables_load=true");
    report.resolved("density_integration", "inactive", "tables_load=true");
    report.resolved("dos", "inactive", "tables_load=true");
    report.resolved("g_positive", "inactive", "tables_load=true");
    report.resolved("g_negative", "inactive", "tables_load=true");
    report.resolved("f_positive", "inactive", "tables_load=true");
    report.resolved("f_negative", "inactive", "tables_load=true");
  } else if (band == "flat") {
    report.value("density_source", "flat");
    report.resolved("density_interpolation", density_method, "flat density is interpolation-independent");
    report.value("density_integration", "exact interpolant primitive");
    report.resolved("dos", "inactive", "flat band");
    report.resolved("g_positive", "inactive", "flat band");
    report.resolved("g_negative", "inactive", "flat band");
    report.resolved("f_positive", "inactive", "flat band");
    report.resolved("f_negative", "inactive", "flat band");
  } else {
    report.value("density_source", "file");
    if (params.contains("density_interpolation"))
      report.value("density_interpolation", density_method);
    else
      report.resolved("density_interpolation", density_method, "parameter default");
    report.value("density_integration", "exact interpolant primitive");
    report.value("dos", Pstr("dos", "Delta.dat"));
    if (adapt) {
      report.value("g_positive", "GSOL.dat");
      report.value("g_negative", "GSOLNEG.dat");
    } else {
      report.resolved("g_positive", "inactive", "adapt=false");
      report.resolved("g_negative", "inactive", "adapt=false");
    }
    report.value("f_positive", "FSOL.dat");
    report.value("f_negative", "FSOLNEG.dat");
  }
  if (mode == TableMode::Calculate) {
    report.value("tables_save", nrgchain_tables_save);
    report.value("tables_load", nrgchain_tables_load);
    report.value("tridiagonalize", nrgchain_tridiag);
  } else {
    report.resolved("tables_save", nrgchain_tables_save, "command mode override");
    report.resolved("tables_load", nrgchain_tables_load, "command mode override");
    report.resolved("tridiagonalize", nrgchain_tridiag, "command mode override");
  }

  if (nrgchain_tables_load)
    report.value("theta_input", "theta.dat");
  else
    report.value("theta_output", output_path("theta.dat"));

  if (nrgchain_tridiag) {
    report.value("xi_output", output_path("xi.dat"));
    report.value("zeta_output", output_path("zeta.dat"));
  } else {
    report.resolved("xi_output", "inactive", "tridiagonalize=false");
    report.resolved("zeta_output", "inactive", "tridiagonalize=false");
  }

  const auto report_table = [&report](const string &key, const string &filename) {
    if (nrgchain_tables_load)
      report.value(key + "_input", filename);
    else if (nrgchain_tables_save)
      report.value(key + "_output", output_path(filename));
    else
      report.resolved(key, "not persisted", "tables_save=false");
  };
  report_table("de_positive_table", "de_pos.dat");
  report_table("de_negative_table", "de_neg.dat");
  report_table("du_positive_table", "du_pos.dat");
  report_table("du_negative_table", "du_neg.dat");
  report.write(cerr);
}
#endif

WilsonData make_result() {
  WilsonData data;
  data.channels.push_back(WilsonChannel{result_theta, result_xi, result_zeta});
  return data;
}

WilsonData run_calculation(const TableMode mode, const filesystem::path &output_dir) {
  nrgchain_output_dir = output_dir;
  set_parameters();
  apply_mode(mode);
  if (nrgchain_tables_load && nrgchain_tables_save)
    throw invalid_argument("nrgchain_tables_load and nrgchain_tables_save cannot both be true.");
#ifndef NRGCHAIN_NO_MAIN
  report_configuration(mode);
#endif

  if (nrgchain_tables_load) {
    load_tables();
  } else {
    calc_tables();
  }

  if (nrgchain_tridiag) tridiag();

  return make_result();
}

} // namespace

WilsonData calculate_from_params(const std::map<std::string, std::string> &param_values, const TableMode mode,
                                 const filesystem::path &output_dir) {
  reset_calculation_state();
  params = param_values;
  return run_calculation(mode, output_dir);
}

WilsonData calculate_from_file(const std::string &param_filename, const TableMode mode,
                               const filesystem::path &output_dir) {
  reset_calculation_state();
  parser(param_filename);
  return run_calculation(mode, output_dir);
}

} // namespace NRG::Tools::NrgChain

#ifndef NRGCHAIN_NO_MAIN
int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "nrgchain")) return EXIT_SUCCESS;
  try {
    clock_t start_clock = clock();

    about();
    const auto options = cmd_line(argc, argv);
    command_line_verbosity = options.verbosity;
    command_line_param_filename = options.param_filename;
    NRG::Tools::NrgChain::calculate_from_file(options.param_filename, options.mode);

    clock_t end_clock = clock();
    cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << endl;
  } catch (const std::exception &e) {
    cerr << "nrgchain: error: " << e.what() << endl;
    return EXIT_FAILURE;
  } catch (...) {
    cerr << "nrgchain: error: unknown exception" << endl;
    return EXIT_FAILURE;
  }
}
#endif
