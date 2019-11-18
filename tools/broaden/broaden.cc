// broaden - Finite-temperature raw spectral data broadening tool
// Ljubljana code Rok Zitko, rok.zitko@ijs.si, 2012-2016

// CHANGE LOG
// 30.1.2013 - option for normalization-conserving broadening kernel
// 25.4.2012 - first version based on the 'bw' tool
// 30.4.2012 - more digits of precision
//           - mesh generated starting from the maximum value
// 31.10.2012 - support for complex y values
// 16.11.2012 - increased precision for output
// 10.10.2013 - finite-temperature sum-rules output
// 16.10.2013 - cumulative spectrum output
// 19.5.2014 - Gaussian broadening
// 6.3.2015 - arbitrary accumulation point at lower end (-a switch)
//          - -l, -h filtering switches
// 20.6.2016 - -P, -N positive/negative filtering
//             -A, -B output mesh control
// 10.10.2016 - final Gaussian broadening pass (-f), or derFD pass (-x)
//            - code cleanup

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

bool verbose         = false; // output verbosity level
bool sumrules        = false; // compute the integrals for testing the T!=0 sum rules
bool cumulative      = false; // output of integrated spectral function
double broaden_min   = 1e-7;  // parameters defining the frequency mesh
double broaden_max   = 2.0;
double broaden_ratio = 1.01;
double alpha;          // broadening parameter alpha
double ggamma, dgamma; // broadening parameter gamma (final Gaussian, derFD)
double T;              // temperature parameter
double omega0_ratio;   // omega0=omega0_ratio*T
double omega0;
double accumulation = 0.0;     // accumpulation point for the mesh
bool one            = false;   // For Nz=1, no subdir.
bool normalization  = false;   // What cross-over function to use?
double filterlow    = 0.0;     // filter all input data points with |omega|<filterlow
double filterhigh   = DBL_MAX; // filter all input data points with |omega|>filterhigh
bool keeppositive   = true;    // Keep omega>0 data points when reading input data
bool keepnegative   = true;    // Keep omega<0 data points when reading input data
bool meshpositive   = true;    // Make omega>0 output mesh
bool meshnegative   = true;    // Make omega<0 output mesh

bool gaussian      = false; // Gaussian broadening scheme
bool finalgaussian = false; // Final pass of Gaussian broadening of width ggamma*T
bool finalderfd    = false; // Final pass of broadening with a derivative of FD distribution

// Multi-column support
int nrcol = 1; // Number of columns
int col   = 1; // Which y column are we interested in?

string name; // filename of binary files containing the raw data
int Nz;      // Number of spectra (1..Nz)

const string output_filename     = "spec.dat";
const string cumulative_filename = "cumulative.dat";

void usage(ostream &F = cout) {
  F << "Usage: broaden <name> <Nz> <alpha> <T> [omega0_ratio]" << endl;
  F << endl;
  F << "Optional parameters:" << endl;
  F << " -v -- verbose" << endl;
  F << " -m <min> -- minimal mesh frequency" << endl;
  F << " -M <max> -- maximal mesh frequency" << endl;
  F << " -r <ratio> -- ratio between two consecutive frequency points" << endl;
  F << " -o -- one .dat file" << endl;
  F << " -2 -- use the 2nd column for weight values (complex spectra)" << endl;
  F << " -3 -- use the 3rd column for weight values (complex spectra)" << endl;
  F << " -n -- normalization-conserving broadening kernel" << endl;
  F << " -s -- compute weighted integrals for testing sum-rules" << endl;
  F << " -c -- compute cumulative spectrum" << endl;
  F << " -g -- Gaussian broadening (width alpha)" << endl;
  F << " -f -- final Gaussian broadening pass" << endl;
  F << " -x -- final derFD broadening pass" << endl;
  F << " -a -- accumulation point for the mesh" << endl;
  F << " -l -- filter out low-frequency raw data" << endl;
  F << " -h -- filter out high-frequency raw data" << endl;
  F << " -P -- keep only positive input frequencies" << endl;
  F << " -N -- keep only negative input frequencies" << endl;
  F << " -A -- output only positive frequencies" << endl;
  F << " -B -- output only negative frequencies" << endl;
}

void cmd_line(int argc, char *argv[]) {
  char c;

  while ((c = getopt(argc, argv, "vm:M:r:o23nscgf:x:a:l:h:PNAB")) != -1) {
    switch (c) {
      case 'c': cumulative = true; break;

      case 's': sumrules = true; break;

      case 'v': verbose = true; break;

      case 'm':
        broaden_min = atof(optarg);
        if (verbose) { cout << "broaden_min=" << broaden_min << endl; }
        break;

      case 'M':
        broaden_max = atof(optarg);
        if (verbose) { cout << "broaden_max=" << broaden_max << endl; }
        break;

      case 'r':
        broaden_ratio = atof(optarg);
        if (verbose) { cout << "broaden_ratio=" << broaden_ratio << endl; }
        break;

      case 'o': one = true; break;

      case '2':
        nrcol = 2;
        col   = 1;
        break;

      case '3':
        nrcol = 2;
        col   = 2;
        break;

      case 'n': normalization = true; break;

      case 'g': gaussian = true; break;

      case 'f':
        finalgaussian = true;
        ggamma        = atof(optarg);
        if (verbose) { cout << "final Gaussian with gamma=" << ggamma << endl; }
        break;

      case 'x':
        finalderfd = true;
        dgamma     = atof(optarg);
        if (verbose) { cout << "final derFD with gamma=" << dgamma << endl; }
        break;

      case 'a':
        accumulation = atof(optarg);
        if (verbose) { cout << "accumulation=" << accumulation << endl; }
        break;

      case 'l':
        filterlow = atof(optarg);
        if (verbose) { cout << "filterlow=" << filterlow << endl; }
        break;

      case 'h':
        filterhigh = atof(optarg);
        if (verbose) { cout << "filterhigh=" << filterhigh << endl; }
        break;

      case 'P': keepnegative = false; break;

      case 'N': keeppositive = false; break;

      case 'A': meshnegative = false; break;

      case 'B': meshpositive = false; break;

      default: abort();
    }
  }

  int remaining = argc - optind; // arguments left

  if (remaining != 5 && remaining != 4) {
    usage();
    exit(1);
  }

  name = string(argv[optind]); // Name of spectral density files

  Nz = atoi(argv[optind + 1]); // Number of z-values
  assert(Nz >= 1);

  alpha = atof(argv[optind + 2]); // High-energy broadening parameter
  assert(alpha > 0.0);

  T = atof(argv[optind + 3]); // Temperature
  assert(T > 0.0);

  if (remaining == 5) {
    omega0_ratio = atof(argv[optind + 4]); // omega0/T
    assert(omega0_ratio > 0.0);
    omega0 = omega0_ratio * T;
  }

  if (remaining == 4) {
    omega0_ratio = 1e-9; // Effectively zero
    omega0       = omega0_ratio * T;
  }

  cout << "Processing: " << name << endl;
  if (verbose) { cout << "Nz=" << Nz << " alpha=" << alpha << " T=" << T << " omega0_ratio=" << omega0_ratio << endl; }
}

double **buffers; // binary data buffers
int *sizes;       // sizes of buffers

typedef map<double, double> mapdd;
using vec = vector<double>;

mapdd spec;           // Spectrum
unsigned int nr_spec; // Number of raw spectrum points
vec vfreq, vspec;     // Same info as spectrum, but in vector<double> form

string tostring(int i) {
  ostringstream S;
  S << i;
  return S.str();
}

// Load a file containing binary representation of raw spectral density.
// The grid is not assumed to be uniform.
void load(int i) {
  // i-th z-value defines the name of the directory where the results of
  // the NRG calculation are contained.
  string filename;
  if (one && Nz == 1) {
    filename = name;
  } else {
    filename = tostring(i) + "/" + name;
  }
  ifstream f(filename.c_str(), ios::in | ios::binary);
  if (!f.good() || f.eof() || !f.is_open()) {
    cerr << "Error opening file " << filename << endl;
    exit(1);
  }
  if (verbose) { cout << "Reading " << filename << endl; }

  const int rows = 1 + nrcol; // number of elements in a line

  // Determine the number of records
  f.seekg(0, ios::beg);
  const ios::pos_type begin_pos = f.tellg();
  f.seekg(0, ios::end);
  const ios::pos_type end_pos = f.tellg();
  const long len              = end_pos - begin_pos;
  assert(len % (rows * sizeof(double)) == 0);
  const int nr = len / (rows * sizeof(double)); // number of lines
  if (verbose) { cout << "len=" << len << " nr=" << nr << " data points" << endl; }

  // Allocate the read buffer. The data will be kept in memory for the
  // duration of the calculation!
  auto *buffer = new double[rows * nr];
  f.seekg(0, ios::beg); // Return to the beginning of the file.
  f.read((char *)buffer, len);
  if (f.fail()) {
    cerr << "Error reading " << filename << endl;
    exit(1);
  }
  f.close();

  // Keep record of the the buffer and its size.
  buffers[i] = buffer;
  sizes[i]   = nr;

  if (verbose) {
    // Check normalization.
    double sum = 0.0;
    for (int j = 0; j < nr; j++) sum += buffer[rows * j + col];
    cout << "Weight=" << sum << endl;
  }
}

// Load all the input data.
void read_files() {
  buffers = new double *[Nz + 1];
  sizes   = new int[Nz + 1];

  for (int i = 1; i <= Nz; i++) { load(i); }
}

// Combine data from all NRG runs (z-averaging).
void merge() {
  const int rows = 1 + nrcol; // number of elements in a line

  // Sum weight corresponding to equal frequencies.  Map of
  // (frequency,weight) pairs is used for this purpose.
  for (int i = 1; i <= Nz; i++) {
    for (int l = 0; l < sizes[i]; l++) {
      double &freq      = buffers[i][rows * l];
      double &value     = buffers[i][rows * l + col];
      auto I = spec.find(freq);
      if (I == spec.end()) {
        spec[freq] = value;
      } else {
        I->second += value;
      }
    }
  }

  nr_spec = spec.size();
  if (verbose) { cout << nr_spec << " unique frequencies." << endl; }

  // Normalize weight by 1/Nz, determine total weight, and store the
  // (frequency,weight) data in the form of linear vectors for faster
  // access in the ensuing calculations.
  double sum = 0.0;
  for (auto & I : spec) {
    const double weight = (I.second /= Nz); // Normalize weight on the fly
    const double freq   = I.first;
    vfreq.push_back(freq);
    vspec.push_back(weight);
    sum += weight;
  }
  if (verbose) { cout << "Total weight=" << sum << endl; }
  assert(vfreq.size() == nr_spec && vspec.size() == nr_spec);
}

vec mesh; // Frequency mesh

vec a; // Spectral function
vec c; // Cumulative spectrum = int_{-inf}^omega a(x)dx.

void filter() {
  for (unsigned int j = 0; j < nr_spec; j++) {
    if (abs(vfreq[j]) < filterlow) { vspec[j] = 0.0; }
    if (abs(vfreq[j]) > filterhigh) { vspec[j] = 0.0; }
  }

  for (unsigned int j = 0; j < nr_spec; j++) {
    if (!keeppositive && vfreq[j] >= 0) { vspec[j] = 0.0; }
    if (!keepnegative && vfreq[j] < 0) { vspec[j] = 0.0; }
  }
}

void integrals_for_sumrules() {
  int nr             = vfreq.size();
  double sum         = 0.0;
  double sumpos      = 0.0;
  double sumneg      = 0.0;
  double sumfermi    = 0.0;
  double sumbose     = 0.0;
  double sumfermiinv = 0.0;
  double sumboseinv  = 0.0;
  for (int i = 0; i < nr; i++) {
    const double omega = vfreq[i];
    const double w     = vspec[i];

    sum += w;
    if (omega > 0.0) { sumpos += w; }
    if (omega < 0.0) { sumneg += w; }
    const double f  = 1 / (1 + exp(-omega / T));
    const double b  = 1 / (1 - exp(-omega / T));
    const double fi = 1 / (1 + exp(+omega / T)); // fi=1-f
    const double bi = 1 / (1 - exp(+omega / T)); // bi=1-b

    if (isfinite(f)) { sumfermi += f * w; }
    if (isfinite(f)) { sumbose += b * w; }
    if (isfinite(fi)) { sumfermiinv += fi * w; }
    if (isfinite(bi)) { sumboseinv += bi * w; }
  }
  cout << "Total weight=" << sum << endl;
  cout << "Positive-omega weight=" << sumpos << endl;
  cout << "Negative-omega weight=" << sumneg << endl;
  cout << "Integral with fermionic kernel=" << sumfermi << endl;
  cout << "Integral with bosonic kernel=" << sumbose << endl;
  cout << "Integral with fermionic kernel (omega -> -omega)=" << sumfermiinv << endl;
  cout << "Integral with bosonic kernel (omega -> -omega)=" << sumboseinv << endl;
}

// Create a mesh on which the output spectral function will be computed.
void make_mesh(vec &mesh) {
  assert(broaden_min < broaden_max);
  assert(broaden_ratio > 1.0);

  const double a = accumulation; // accumulation point for the mesh

  for (double z = broaden_max; z > broaden_min; z /= broaden_ratio) {
    double x = z * (broaden_max - a) / broaden_max + a;
    if (meshpositive) mesh.push_back(x);
    if (meshnegative) mesh.push_back(-x);
  }

  sort(mesh.begin(), mesh.end());
}

#ifndef M_SQRTPI
#define M_SQRTPI 1.7724538509055160273
#endif

inline double sqr(const double x) { return x * x; }

// Modified log-Gaussian broadening kernel. For gamma=alpha/4, the
// kernel is symmetric in both arguments.
// e is the energy of the spectral function point being computed.
// ept is the energy of point.
inline double BR_L_orig(const double e, const double ept) {
  if ((e < 0.0 && ept > 0.0) || (e > 0.0 && ept < 0.0)) return 0.0;

  if (ept == 0.0) return 0.0;

  const double gamma = alpha / 4;
  return exp(-sqr(log(e / ept) / alpha - gamma)) / (alpha * abs(e) * M_SQRTPI);
}

// As above, with support for a shifted accumulation point
inline double BR_L_acc(const double e0, const double ept0) {
  double e   = e0;
  double ept = ept0;

  if (e > accumulation && ept > accumulation) {
    const double shift = accumulation;
    e                  = e - shift;
    ept                = ept - shift;
  }

  if (e < -accumulation && ept < -accumulation) {
    const double shift = -accumulation; // note the sign!
    e                  = e - shift;
    ept                = ept - shift;
  }

  if ((e < 0.0 && ept > 0.0) || (e > 0.0 && ept < 0.0)) return 0.0;

  if (ept == 0.0) return 0.0;

  const double gamma = alpha / 4;
  return exp(-sqr(log(e / ept) / alpha - gamma)) / (alpha * abs(e) * M_SQRTPI);
}

#define BR_L BR_L_acc

// Normalized to 1, width omega0. The kernel is symmetric in both
// arguments.
inline double BR_G(const double e, const double ept) { return exp(-sqr((e - ept) / omega0)) / (omega0 * M_SQRTPI); }

inline double BR_G_alpha(const double e, const double ept) {
  const double width = alpha;
  return exp(-sqr((e - ept) / width)) / (width * M_SQRTPI);
}

inline double BR_h0(const double x) {
  const double absx = abs(x);
  if (absx > omega0) {
    return 1.0;
  } else {
    return exp(-sqr(log(absx / omega0) / alpha));
  }
}

inline double BR_h(const double e, const double ept) {
  double h;
  if (normalization) {
    // Cross-over funciton as proposed by A. Weichselbaum et al.
    h = BR_h0(ept);
  } else {
    // Note: this is DIFFERENT from the broadening kernel proposed by
    // by Weichselbaum et al. This BR_h is a function of 'e', not
    // of 'ept'! This breaks the normalization at finite temperatures.
    // On the other hand, it gives nicer spectra, especially in
    // combination with the self-energy trick.
    h = BR_h0(e);
  }
  return h;
}

// e - output energy
// ept - energy of the delta peak (data point)
double bfnc(const double e, const double ept) {
  if (gaussian) {
    return BR_G_alpha(e, ept);
  } else {
    const double part_l = BR_L(e, ept);
    const double part_g = BR_G(e, ept);
    const double h      = BR_h(e, ept);
    assert(h >= 0.0 && h <= 1.0);
    return part_l * h + part_g * (1.0 - h);
  }
}

// Do the broadening
void broaden(const vec &mesh, vec &a) {
  const int nr_mesh = mesh.size();

  if (verbose) { cout << "Broadening. nr_mesh=" << nr_mesh << endl; }

  a.resize(nr_mesh);

  for (int i = 0; i < nr_mesh; i++) {
    const double outputfreq = mesh[i];
    a[i]                    = 0.0; // clear!
    for (unsigned int j = 0; j < nr_spec; j++) { a[i] += vspec[j] * bfnc(outputfreq, vfreq[j]); }
  }
}

const double R_1_SQRT2PI = 1.0 / sqrt(2.0 * M_PI);

double gaussian_kernel(const double x, const double y, const double sigma) {
  const double d = (x - y) / sigma;
  return R_1_SQRT2PI * exp(-d * d / 2.0) / sigma;
}

// Derivative of the Fermi-Dirac function: -d/dw f_FD(x-y) with
// effective temperature sigma.
double derfd_kernel(const double x, const double y, const double sigma) {
  const double d = (x - y) / sigma;
  return 1.0 / ((1.0 + cosh(d)) * 2.0 * sigma);
}

const double cutoff_ratio = 100;

void gaussian_convolve(const vec &mesh, vec &a, const double sigma) {
  const int nr_mesh = mesh.size();

  if (verbose) { cout << "Convolving with a Gaussian. sigma=" << sigma << endl; }

  vec b(a); // source

  a[0]           = b[0];
  a[nr_mesh - 1] = b[nr_mesh - 1];

  for (int i = 1; i < nr_mesh - 1; i++) {
    const double x = mesh[i];

    if (abs(x) < cutoff_ratio * sigma) {
      double sum = 0.0, sum_kernel = 0.0;
      for (unsigned int j = 1; j < nr_mesh - 1; j++) {
        const double y     = mesh[j];
        const double width = (mesh[j + 1] + mesh[j]) / 2.0 - (mesh[j - 1] + mesh[j]) / 2.0;
        const double kw    = gaussian_kernel(x, y, sigma) * width;
        assert(isfinite(kw));
        sum += b[j] * kw;
        sum_kernel += kw;
      }
      a[i] = sum;
    } else {
      // High temperatures: convolution not necessary
      a[i] = b[i];
    }
  }
}

void derfd_convolve(const vec &mesh, vec &a, const double sigma) {
  const int nr_mesh = mesh.size();

  if (verbose) { cout << "Convolving with a derivative of the Fermi-Dirac distribution. sigma=" << sigma << endl; }

  vec b(a); // source

  a[0]           = b[0];
  a[nr_mesh - 1] = b[nr_mesh - 1];

  for (int i = 1; i < nr_mesh - 1; i++) {
    const double x = mesh[i];

    if (abs(x) < cutoff_ratio * sigma) {
      double sum = 0.0, sum_kernel = 0.0;
      for (unsigned int j = 1; j < nr_mesh - 1; j++) {
        const double y     = mesh[j];
        const double width = (mesh[j + 1] + mesh[j]) / 2.0 - (mesh[j - 1] + mesh[j]) / 2.0;
        const double kw    = derfd_kernel(x, y, sigma) * width;
        assert(isfinite(kw));
        sum += b[j] * kw;
        sum_kernel += kw;
      }
      a[i] = sum;
    } else {
      // High temperatures: convolution not necessary
      a[i] = b[i];
    }
  }
}

// Cumulative spectrum
void calc_cumulative(const vec &mesh, vec &c) {
  const int nr_mesh = mesh.size();

  if (verbose) cout << "Calculating cumulative spectrum." << endl;

  c.resize(nr_mesh);

  double sum = 0.0;
  int j      = 0;
  for (int i = 0; i < nr_mesh; i++) {
    const double max_freq = mesh[i];

    while (vfreq[j] < max_freq) {
      sum += vspec[j];
      j++;
    }

    c[i] = sum;
  }

  cout << "End sum=" << sum << endl;
}

// Increased 12->18, Nov 16 2012
const int SAVE_PREC = 18; // Precision for output to the file
const int COUT_PREC = 18; // Precision for verbose reporting on console

// Save a map of (double,double) pairs to a file.
void save(const string filename, const vec &x, const vec &y) {
  if (verbose) { cout << "Saving " << filename << endl; }

  ofstream F(filename.c_str());
  if (!F) {
    cerr << "Failed to open " << filename << " for writing." << endl;
    exit(1);
  }

  F << setprecision(SAVE_PREC);

  assert(x.size() == y.size());
  unsigned int nr = x.size();
  for (unsigned int i = 0; i < nr; i++) { F << x[i] << " " << y[i] << endl; }
}

// Estimate the weight using the trapezoidal rule.
// The x array must be sorted.
double trapez(const vec &x, const vec &y) {
  double weight = 0.0;

  assert(x.size() == y.size());
  unsigned int nr = x.size();
  for (unsigned int i = 1; i < nr; i++) {
    assert(x[i] >= x[i - 1]);
    weight += (y[i - 1] + y[i]) / 2 * (x[i] - x[i - 1]);
  }

  return weight;
}

void check_normalizations(const vec &x) {
  vec y(x.size());

  // For a range of delta peak positions, compute the total weight
  // of the resulting broadened spectra.
  for (double z = 1e-10; z < 1.0; z *= 10) {
    unsigned int nr = x.size();
    for (unsigned int i = 0; i < nr; i++) { y[i] = bfnc(x[i], z); }

    double weight = trapez(x, y);

    cout << "z=" << z << " weight=" << weight << endl;
  }
}

int main(int argc, char *argv[]) {
  cout << "broaden - finite-temperature broadening tool" << endl;
  cout << setprecision(COUT_PREC);

  cmd_line(argc, argv);
  read_files();
  merge();
  filter();

  if (sumrules) integrals_for_sumrules();

  make_mesh(mesh);

  if (verbose) check_normalizations(mesh);

  broaden(mesh, a);

  if (finalgaussian) gaussian_convolve(mesh, a, ggamma * T);

  if (finalderfd) derfd_convolve(mesh, a, dgamma * T);

  double weight = trapez(mesh, a);
  cout << "Estimated weight (trapezoidal rule)=" << weight << endl;

  save(output_filename, mesh, a);

  if (cumulative) {
    calc_cumulative(mesh, c);
    save(cumulative_filename, mesh, c);
  }
}
