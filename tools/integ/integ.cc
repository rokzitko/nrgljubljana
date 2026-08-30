// Numerical integration tool: compute accurate spectral weights and moments
// of tabulated functions using smooth interpolation functions (GSL).
// Part of "NRG Ljubljana"
// Rok Zitko, rok.zitko@ijs.si

// The input file must consist of a table of space-separated (energy,
// value) pairs. Gauss-Kronrod quadrature rules are used.

// CHANGE LOG
// 21.5.2014 - first version

#include <cstddef>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ios>
#include <istream>
#include <limits>
#include <ostream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <getopt.h>

#include <common/version.hpp>

#include "../common/diagnostics.hpp"
#include "../common/gsl_config.hpp"

using namespace std;

using NRG::Tools::GslErrorPolicy;
using NRG::Tools::InterpolationMethod;
using NRG::Tools::QagRule;

typedef pair<double, double> XYPOINT;
using XYFUNC = vector<XYPOINT>;
using DVEC = vector<double>;

// number of digits of precision in the output
#define OUTPUT_PRECISION 16

int verbosity = 0;

double T = 1e-99; // Temperature. Default is (essentially) 0.
string inputfn;   // Filename for input data

double sum;                              // integral over input function on [Xmin:Xmax] interval
double total, totalabs, pos, neg, fermi; // Results
string out = "total";                    // What to output to STDOUT

void about(ostream &output = cout) {
  output << "Integration tool - GSL" << endl;
}

void usage() {
  cout << "\nUsage: integ [options] <input> [-p|n|a|f]" << endl;
  cout << "-h, --help: show help" << endl;
  cout << "-v: show resolved configuration and verbose diagnostics on standard error" << endl;
  cout << "-vv: also show detailed integration diagnostics" << endl;
  cout << "-V, --version: show project version" << endl;
  cout << "-w: ignore QAG errors (alias for --gsl-error-policy ignore)" << endl;
  cout << "-T: temperature T" << endl;
  cout << "-p: integral over positive X range" << endl;
  cout << "-n: integral over negative X range" << endl;
  cout << "-a: integral over |f|" << endl;
  cout << "-f: integral weighted with Fermi-Dirac function for temperature T" << endl;
  cout << "-i, --interpolation linear|cspline|akima|steffen: interpolation method (default=akima)" << endl;
  cout << "--epsabs value: absolute integration tolerance (default=1e-12)" << endl;
  cout << "--epsrel value: relative integration tolerance (default=1e-8)" << endl;
  cout << "--workspace-limit value: integration workspace size (default=1000)" << endl;
  cout << "--quadrature-rule 15|21|31|41|51|61: Gauss-Kronrod rule (default=15)" << endl;
  cout << "--gsl-error-policy ignore|warn|fail: QAG error policy (default=warn)" << endl;
}

InterpolationMethod interpolation_method = InterpolationMethod::akima;
double epsabs                             = 1e-12;
double epsrel                             = 1e-8;
size_t limit                              = 1000;
QagRule quadrature_rule                   = QagRule::gauss15;
GslErrorPolicy gsl_error_policy           = GslErrorPolicy::warn;

void cmd_line(int argc, char *argv[]) {
  enum LongOption {
    OPT_EPSABS = 256,
    OPT_EPSREL,
    OPT_WORKSPACE_LIMIT,
    OPT_QUADRATURE_RULE,
    OPT_GSL_ERROR_POLICY,
  };
  static const struct option long_options[] = {
    {"help", no_argument, nullptr, 'h'},
    {"interpolation", required_argument, nullptr, 'i'},
    {"epsabs", required_argument, nullptr, OPT_EPSABS},
    {"epsrel", required_argument, nullptr, OPT_EPSREL},
    {"workspace-limit", required_argument, nullptr, OPT_WORKSPACE_LIMIT},
    {"quadrature-rule", required_argument, nullptr, OPT_QUADRATURE_RULE},
    {"gsl-error-policy", required_argument, nullptr, OPT_GSL_ERROR_POLICY},
    {nullptr, 0, nullptr, 0},
  };

  int c;
  while ((c = getopt_long(argc, argv, "hi:vwt:T:pnaf", long_options, nullptr)) != -1) {
    switch (c) {
      case 'h':
        usage();
        exit(EXIT_SUCCESS);
      
      case 'v': ++verbosity; break;

      case 'w': gsl_error_policy = GslErrorPolicy::ignore; break;

      case 't': // case insensitive!
      case 'T':
        T = atof(optarg);
        if (verbosity >= 1) { cerr << "T=" << T << endl; }
        break;

      case 'p': out = "pos"; break;

      case 'n': out = "neg"; break;

      case 'a': out = "abs"; break;

      case 'f': out = "fermi"; break;

      case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;

      case OPT_EPSABS: epsabs = NRG::Tools::parse_finite_double(optarg, "Absolute integration tolerance"); break;

      case OPT_EPSREL: epsrel = NRG::Tools::parse_finite_double(optarg, "Relative integration tolerance"); break;

      case OPT_WORKSPACE_LIMIT: limit = NRG::Tools::parse_positive_size(optarg, "Workspace limit"); break;

      case OPT_QUADRATURE_RULE: quadrature_rule = NRG::Tools::parse_qag_rule(optarg); break;

      case OPT_GSL_ERROR_POLICY: gsl_error_policy = NRG::Tools::parse_gsl_error_policy(optarg); break;

      default: exit(EXIT_FAILURE);
    }
  }

  NRG::Tools::validate_tolerances(epsabs, epsrel);
  NRG::Tools::validate_qag_workspace_limit(limit);

  int remaining = argc - optind;

  if (remaining != 1) {
    about();
    usage();
    exit(1);
  }

  inputfn = string(argv[optind]);
}

// Read data from stream F.
void readtable(istream &F, XYFUNC &v) {
  while (F) {
    if (F.peek() == '#') { // skip comment lines
      string line;
      getline(F, line);
    } else {
      double x, y;
      F >> x >> y;

      if (F.fail()) break;

      assert(std::isfinite(x) && std::isfinite(y));
      v.push_back(make_pair(x, y));
    }
  }

  if (verbosity >= 1) cerr << v.size() << " lines read." << endl;
}

int len;           // number of data points
double Xmin, Xmax; // the interval boundaries

DVEC Xpts, Ypts;

gsl_interp_accel *acc;
gsl_spline *spline;

gsl_integration_workspace *w;

void report_configuration() {
  if (verbosity == 0) return;
  NRG::Tools::ConfigurationReport report("integ");
  report.value("verbosity", verbosity);
  report.value("input.file", inputfn);
  report.value("input.points", len);
  report.resolved("input.lower_bound", Xmin, "sorted input data");
  report.resolved("input.upper_bound", Xmax, "sorted input data");
  report.value("output.mode", out);
  report.value("output.precision", OUTPUT_PRECISION);
  report.value("temperature", T);
  report.value("interpolation", NRG::Tools::interpolation_method_name(interpolation_method));
  report.value("integration.epsabs", epsabs);
  report.value("integration.epsrel", epsrel);
  report.value("integration.workspace_limit", limit);
  report.value("integration.quadrature_rule", static_cast<int>(quadrature_rule));
  report.value("integration.gsl_error_policy", NRG::Tools::gsl_error_policy_name(gsl_error_policy));
  report.write(cerr);
}

void init(XYFUNC &im) {
  if (im.empty()) {
    cerr << "Error: no data points found." << endl;
    exit(1);
  }

  const auto minimum_size = NRG::Tools::interpolation_minimum_size(interpolation_method);
  if (im.size() < minimum_size) {
    cerr << "Error: failed to initialize GSL spline: " << NRG::Tools::interpolation_method_name(interpolation_method)
         << " interpolation requires at least " << minimum_size << " data points." << endl;
    exit(1);
  }

  sort(im.begin(), im.end());

  len = im.size();

  Xmin = im[0].first;
  Xmax = im[len - 1].first;

  if (verbosity >= 1) cerr << "Range: [" << Xmin << " ; " << Xmax << "]" << endl;

  // Xpts are increasing
  Xpts = DVEC(len);
  Ypts = DVEC(len);

  for (int i = 0; i < len; i++) {
    Xpts[i] = im[i].first;
    Ypts[i] = im[i].second;
  }

  report_configuration();

  gsl_set_error_handler_off();
  acc = gsl_interp_accel_alloc();
  if (!acc) {
    cerr << "Error: failed to allocate GSL interpolation accelerator." << endl;
    exit(1);
  }
  const gsl_interp_type *Interp_type = NRG::Tools::gsl_interpolation_type(interpolation_method);
  spline                             = gsl_spline_alloc(Interp_type, len);
  if (!spline) {
    cerr << "Error: failed to allocate GSL spline." << endl;
    exit(1);
  }
  if (const auto status = gsl_spline_init(spline, &Xpts[0], &Ypts[0], len); status != 0) {
    cerr << "Error: failed to initialize GSL spline: " << gsl_strerror(status) << endl;
    exit(1);
  }

  sum = gsl_spline_eval_integ(spline, Xmin, Xmax, acc);

  if (!std::isfinite(sum)) {
    cerr << "Error: Integral is not a finite number." << endl;
    exit(1);
  }

  if (verbosity >= 1) cerr << "Sum=" << sum << endl;

  w = gsl_integration_workspace_alloc(limit);
  if (!w) {
    cerr << "Error: failed to allocate GSL integration workspace." << endl;
    exit(1);
  }
}

inline double f_neg(double X, [[maybe_unused]] void *params) { return (X < 0 ? gsl_spline_eval(spline, X, acc) : 0); }

inline double f_pos(double X, [[maybe_unused]] void *params) { return (X > 0 ? gsl_spline_eval(spline, X, acc) : 0); }

inline double f_total(double X, [[maybe_unused]] void *params) { return gsl_spline_eval(spline, X, acc); }

inline double f_abs(double X, [[maybe_unused]]  void *params) { return fabs(gsl_spline_eval(spline, X, acc)); }

inline double f_fermi(double X, [[maybe_unused]] void *params) {
  double fd = 1.0 / (1.0 + exp(X / T));
  return gsl_spline_eval(spline, X, acc) * fd;
}

void handle_qag(int status, double integral, double integration_error) {
  if (!NRG::Tools::gsl_integration_failed(status, integral, integration_error) || gsl_error_policy == GslErrorPolicy::ignore) return;

  const char *severity = gsl_error_policy == GslErrorPolicy::fail ? "ERROR" : "WARNING";
  if (status != GSL_SUCCESS) cerr << severity << " - qag error: " << status << " -- " << gsl_strerror(status) << endl;
  if (!std::isfinite(integral) || !std::isfinite(integration_error))
    cerr << severity << " - qag returned a non-finite result or error estimate." << endl;

  if (gsl_error_policy == GslErrorPolicy::fail) exit(EXIT_FAILURE);
}

// NOTE about Gauss-Kronrod: The higher-order rules give better accuracy
// for smooth functions, while lower-order rules save time when the
// function contains local difficulties, such as discontinuities. [GSL manual]

// On each iteration the adaptive integration strategy bisects the
// interval with the largest error estimate. The subintervals and their
// results are stored in the memory provided by workspace. The maximum
// number of subintervals is given by limit, which may not exceed the
// allocated size of the workspace. [GSL manual]

double calc(double (*fnc)(double, void *)) {
  gsl_function F;
  F.function = fnc;
  F.params   = nullptr;

  double integral          = std::numeric_limits<double>::quiet_NaN();
  double integration_error = std::numeric_limits<double>::quiet_NaN();
  int status = gsl_integration_qag(&F,   // integrand function
                                   Xmin, // lower integration boundary
                                   Xmax, // upper integration boundary
                                   epsabs, epsrel,
                                   limit,                                      // size of workspace w
                                   NRG::Tools::gsl_qag_rule(quadrature_rule),   // Gauss-Kronrod rule
                                   w,                   // integration workspace
                                   &integral,           // final approximation
                                   &integration_error); // estimate of absolute error

  handle_qag(status, integral, integration_error);

  if (verbosity >= 2) {
    cerr << scientific;
    cerr << "Result=" << integral << endl;
    cerr << "Int. error=" << integration_error << endl;
    cerr.unsetf(ios_base::scientific);
  }

  return integral;
}

void done() {
  gsl_spline_free(spline);
  gsl_interp_accel_free(acc);
  gsl_integration_workspace_free(w);
}

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "integ")) return EXIT_SUCCESS;
  try {
    cmd_line(argc, argv);
    if (verbosity >= 1) about(cerr);

    if (verbosity >= 1) cerr << "T=" << T << endl;

    ifstream Fin;

    Fin.open(inputfn.c_str());
    if (!Fin) {
      cerr << "Can't open " << inputfn << " for reading." << endl;
      exit(2);
    }

    XYFUNC f;
    readtable(Fin, f);
    init(f);
    total    = calc(f_total);
    pos      = calc(f_pos);
    neg      = calc(f_neg);
    totalabs = calc(f_abs);
    fermi    = calc(f_fermi);
    done();

    if (verbosity >= 1) {
      cerr << "Total=" << total << endl;
      cerr << "Positive=" << pos << endl;
      cerr << "Negative=" << neg << endl;
      cerr << "Total|f|=" << totalabs << endl;
      cerr << "Fermi-Dirac weighted=" << fermi << endl;
    }

    cout << setprecision(OUTPUT_PRECISION);
    if (out == "total") cout << total << endl;
    if (out == "pos") cout << pos << endl;
    if (out == "neg") cout << neg << endl;
    if (out == "abs") cout << totalabs << endl;
    if (out == "fermi") cout << fermi << endl;
    return EXIT_SUCCESS;
  } catch (const exception &error) {
    cerr << "Error: " << error.what() << endl;
    return EXIT_FAILURE;
  }
}
