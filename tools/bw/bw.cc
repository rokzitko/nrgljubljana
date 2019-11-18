// bw - Adaptive broadening of spectral functions obtained using NRG.
// Ljubljana code Rok Zitko, rok.zitko@ijs.si, Mar, Aug 2009, Jun 2010

// Approach from A. Freyn, S. Florens, PRB 79, 121102(R) (2009).

// CHANGE LOG
// 24.9.2009 - bmin - minimal b(w)
//           - dynamic mesh
// 28.9.2009 - option 'one'
// 12.1.2010 - missing header added
// 14.6.2010 - rescale 'avg' by q in calc_b()
// 22.10.2010 - removed some overzealous assertion checks

#define VERSION "0.2.5"

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

bool verbose            = false; // output verbosity level
bool veryverbose        = false; // horribly detailed output
double broaden_min      = 1e-7;  // parameters defining the frequency mesh
double broaden_max      = 2.0;
double broaden_ratio    = 1.01;
int nr_iter             = 10;    // number of iterations to perform
double trim             = -999;  // the lowest (abs) frequency for output
bool enforce_positivity = false; // Enforce A(w) to be positive
double b0;                       // initial broadening parameter, b_0
double bmin     = -999;          // the smallest b(omega) allowed
double q        = 1.0;           // regularization parameter
bool savemore   = false;         // save intermediate results for A(w) and b(w)
bool saveall    = false;         // save all quantities (integrated DOS, etc.)
double dyn_mesh = -999;          // if > 0, the mesh is dynamically refined
bool one        = false;         // For Nz=1, no subdir.

string name;      // filename of binary files containing the raw data
int Nz;           // Number of spectra (1..Nz)
double **buffers; // binary data buffers
int *sizes;       // sizes of buffers

typedef map<double, double> mapdd;
using vec = vector<double>;

mapdd spec;           // Spectrum
unsigned int nr_spec; // Number of raw spectrum points
vec vfreq, vspec;     // Same info as spectrum, but in vector<double> form
mapdd intspec;        // Integrated spectrum

vec b; // Broadening coefficients (associated with data points!)

vec mesh; // Frequency mesh

vec a;                // Spectral function [current approximation]
vec inta, intb, intc; // Integrated spectral functions

void usage(ostream &F = cout) { F << "Usage: bw <name> <b0> <Nz>\n"; }

void cmd_line(int argc, char *argv[]) {
  char c;

  while ((c = getopt(argc, argv, "vVm:M:r:n:t:pq:sSx:d:o")) != -1) {
    switch (c) {
      case 'v': verbose = true; break;

      case 'V': veryverbose = true; break;

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

      case 'n':
        nr_iter = atoi(optarg);
        if (verbose) { cout << "nr_iter=" << nr_iter << endl; }
        break;

      case 't':
        trim = atof(optarg);
        if (verbose) { cout << "trim=" << trim << " trim/broaden_min=" << trim / broaden_min << endl; }
        break;

      case 'p':
        enforce_positivity = true;
        if (verbose) { cout << "positivity of the spectral function will be enforced" << endl; }
        break;

      case 'q':
        q = atof(optarg);
        if (verbose) { cout << "regularization q=" << q << endl; }
        break;

      case 's': savemore = true; break;

      case 'S': saveall = true; break;

      case 'x':
        bmin = atof(optarg);
        if (verbose) { cout << "minimal b(w)=" << bmin << endl; }
        break;

      case 'd':
        dyn_mesh = atof(optarg);
        if (verbose) { cout << "dyn_mesh=" << dyn_mesh << endl; }
        break;

      case 'o': one = true; break;

      default: abort();
    }
  }

  int remaining = argc - optind; // arguments left

  if (remaining != 3) {
    usage();
    exit(1);
  }
  name = string(argv[optind]);   // Name of spectral density files
  b0   = atof(argv[optind + 1]); // Parameter b0 (asymptotic broadening)
  assert(b0 > 0.0);
  Nz = atoi(argv[optind + 2]); // Number of z-values
  assert(Nz >= 1);

  cout << "Processing: " << name << endl;
  cout << "b0=" << b0 << " Nz=" << Nz << endl;
}

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

  // Determine the number of records
  f.seekg(0, ios::beg);
  const ios::pos_type begin_pos = f.tellg();
  f.seekg(0, ios::end);
  const ios::pos_type end_pos = f.tellg();
  const long len              = end_pos - begin_pos;
  assert(len % (2 * sizeof(double)) == 0);
  const int nr = len / (2 * sizeof(double)); // number of pairs of double
  if (verbose) { cout << "len=" << len << " nr=" << nr << " data points" << endl; }

  // Allocate the read buffer. The data will be kept in memory for the
  // duration of the calculation!
  auto *buffer = new double[2 * nr];
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
    // Check normalization to 1.
    double sum = 0.0;
    for (int j = 0; j < nr; j++) sum += buffer[2 * j + 1];
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
  // Sum weight corresponding to equal frequencies.  Map of
  // (frequency,weight) pairs is used for this purpose.
  for (int i = 1; i <= Nz; i++) {
    for (int l = 0; l < sizes[i]; l++) {
      double &freq      = buffers[i][2 * l];
      double &value     = buffers[i][2 * l + 1];
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

// Calculate integrated spectral function of (a) defined on the mesh (mesh)
// and store the results in a vector (inta). Trapezoid rule is used to
// perform the integration.
void integrate_a(const vec &a, const vec &mesh, vec &inta) {
  assert(a.size() == mesh.size());

  double sum   = 0.0;
  const int nr = a.size();
  inta.resize(nr);
  inta[0] = 0.0;

  for (int i = 1; i < nr; i++) {
    assert(mesh[i] > mesh[i - 1]);
    sum += (a[i] + a[i - 1]) * (mesh[i] - mesh[i - 1]) / 2.0; // Trapezoid rule
    inta[i] = sum;                                            // Integrated spectrum up to current freq
  }

  if (verbose) { cout << "Total weight=" << sum << endl; }
}

// inta = \int [-infty, omega]
// intb = \int [0, omega]
// intc = \int [+infty, omega]
void combinations(const vec &mesh, const vec &inta, vec &intb, vec &intc) {
  const int nr = mesh.size();
  intb.resize(nr);
  intc.resize(nr);

  // Determine the value at omega=0
  double omega0 = -1;
  for (int i = 1; i < nr; i++) {
    if (mesh[i - 1] < 0.0 && mesh[i] > 0.0) {
      omega0 = (inta[i - 1] + inta[i]) / 2.0;
      break;
    }
  }

  // The value at +infty
  const double omegainf = inta[nr - 1];

  for (int i = 0; i < nr; i++) {
    intb[i] = inta[i] - omega0;
    intc[i] = omegainf - inta[i];
  }
}

// Create a mesh on which the output spectral function will be computed.
void make_mesh(vec &mesh) {
  assert(broaden_min < broaden_max);
  assert(broaden_ratio > 1.0);

  double z;
  for (z = broaden_min; z < broaden_max; z *= broaden_ratio) {
    mesh.push_back(z);
    mesh.push_back(-z);
  }

  // One more point to ensure that the point 'broaden_max' is part of the
  // full frequency interval.
  z *= broaden_ratio;
  mesh.push_back(z);
  mesh.push_back(-z);

  if (verbose) { cout << "Maximal frequency=" << z << endl; }

  sort(mesh.begin(), mesh.end());
}

// Create an initial approximation for the frequency dependent broadening
// function b(omega).
void initial_b(vec &b) {
  b.resize(nr_spec);
  for (unsigned int i = 0; i < nr_spec; i++) { b[i] = b0; }
}

// The modified log-Gaussian broadening function.
// E is the output energy, omega the energy of a delta-peak contribution.
// Normalization issue: since alpha is associated with a particular
// delta peak, the broadened function is still normalized to one since
// alpha is a constant [for that particular contribution].
double bfnc(double E, double alpha, double omega) {
  const double div = alpha * abs(E) * sqrt(M_PI);

  const double gamma = alpha / 4.0;
  const double ln    = log(E / omega);

  const double eksponent = -pow(ln / alpha - gamma, 2);

  const double res = exp(eksponent) / div;

  return res;
}

void refine_mesh(vec &mesh, const vec &a) {
  assert(mesh.size() == a.size());
  assert(mesh.size() >= 6); // minimal meaningful mesh size

  const int nr = mesh.size();

  vec newmesh;

  // Negative part of the spectrum
  newmesh.push_back(mesh[0]);
  double newpt1 = -sqrt(mesh[0] * mesh[1]);
  if (mesh[0] < newpt1 && newpt1 < mesh[1]) newmesh.push_back(newpt1); // always refine here, if still possible!
  newmesh.push_back(mesh[1]);
  for (int i = 2; i < nr && mesh[i] < 0.0; i++) {
    const double deriv       = (a[i - 2] - a[i - 1]) / (mesh[i - 2] - mesh[i - 1]);
    const double y_predicted = a[i - 1] + deriv * (mesh[i] - mesh[i - 1]);
    const double y           = a[i];
    const double error_ratio = abs((y - y_predicted) / y);
    if (error_ratio > dyn_mesh) {
      // add one additional mesh point!
      // geometric average!
      double newpt = -sqrt(mesh[i] * mesh[i - 1]);
      // Only add a point if it is actually distinct from its
      // neighbours. This check is necessary in the presence of
      // singularities.
      if (mesh[i - 1] < newpt && newpt < mesh[i]) newmesh.push_back(newpt);
    }
    newmesh.push_back(mesh[i]);
  }

  // Positive part of the spectrum
  newmesh.push_back(mesh[nr - 1]);
  double newpt2 = sqrt(mesh[nr - 1] * mesh[nr - 2]);
  if (mesh[nr - 2] < newpt2 && newpt2 < mesh[nr - 1]) newmesh.push_back(newpt2); // always refine here!
  newmesh.push_back(mesh[nr - 2]);
  for (int i = (nr - 1) - 2; i > 0 && mesh[i] > 0.0; i--) {
    const double deriv       = (a[i + 2] - a[i + 1]) / (mesh[i + 2] - mesh[i + 1]);
    const double y_predicted = a[i + 1] + deriv * (mesh[i] - mesh[i + 1]);
    const double y           = a[i];
    const double error_ratio = abs((y - y_predicted) / y);
    if (error_ratio > dyn_mesh) {
      double newpt = +sqrt(mesh[i] * mesh[i + 1]);
      if (mesh[i] < newpt && newpt < mesh[i + 1]) newmesh.push_back(newpt);
    }
    newmesh.push_back(mesh[i]);
  }

  sort(newmesh.begin(), newmesh.end());

  // Trick: avoid copying!
  mesh.swap(newmesh);
}

// Note: vector 'a' is resized to match the length of vector 'mesh'. The
// mesh can thus dynamically change from iteration to iteration.
void broaden(const vec &mesh, vec &a, const vec &b) {
  const int nr_mesh = mesh.size();

  if (verbose) { cout << "Broadening. nr_mesh=" << nr_mesh << endl; }

  a.resize(nr_mesh);

  for (int i = 0; i < nr_mesh; i++) {
    const double &outputfreq = mesh[i];
    a[i]                     = 0.0;
    for (unsigned int j = 0; j < nr_spec; j++) {
      const double &omega = vfreq[j];
      bool sameSign       = (outputfreq < 0.0) == (omega < 0.0);
      if (sameSign) { a[i] += vspec[j] * bfnc(outputfreq, b[j], omega); }
    }

    // Enforce positivity
    if (a[i] < 0.0 && enforce_positivity) {
      if (veryverbose) { cout << "Warning: a(" << outputfreq << ")=" << a[i] << endl; }
      a[i] = 1e-16;
    }
  }
}

// Save a map of (double,double) pairs to a file.
void save(const string filename0, const mapdd &m, int iter = 0, double trim = 0.0) {
  const string filename = filename0 + (iter > 0 ? tostring(iter) : "") + ".dat";

  if (verbose) { cout << "Saving " << filename << endl; }

  ofstream F(filename.c_str());
  if (!F) {
    cerr << "Failed to open " << filename << " for writing." << endl;
    exit(1);
  }

  for (auto I : m) {
    if (trim != 0.0 && abs(I.first) < trim) { continue; }
    F << I.first << " " << I.second << endl;
  }
}

// Save pairs taken respectively from (mesh) and (data).
void save(const string filename0, const vec &mesh, const vec &data, int iter = 0, double trim = 0.0) {
  const string filename = filename0 + (iter > 0 ? tostring(iter) : "") + ".dat";
  if (verbose) { cout << "Saving " << filename << endl; }

  ofstream F(filename.c_str());
  if (!F) {
    cerr << "Failed to open " << filename << " for writing." << endl;
    exit(1);
  }

  assert(mesh.size() == data.size());
  int nr = mesh.size();
  for (int i = 0; i < nr; i++) {
    const double x = mesh[i];
    const double y = data[i];
    if (trim != 0.0 && abs(x) < trim) { continue; }
    if (isfinite(data[i])) {
      F << x << " " << y << endl;
    } else {
      if (veryverbose) { cout << "Warning: " << i << " not finite." << endl; }
    }
  }
}

// Eq. (7) in the paper of Freyn & Florens
void calc_b(const vec &logd1, const vec &logd2, vec &bb) {
  const unsigned int n = logd1.size();
  bb.resize(n);

  for (unsigned int i = 0; i < n; i++) {
    const double part1 = 1.0 / (q + abs(logd1[i]));
    const double part2 = 1.0 / (q + abs(logd2[i]));
    const double avg   = q * (part1 + part2) / 2.0;
    bb[i]              = b0 * avg;
  }
}

// Calculate a logarithmic derivative of (inta).
// We take dlog f/dlog omega = Delta (log f)/Delta (log omega)
// = (log f_1-log f_0)/(log omega_1-log omega_0)
// = log(f_1/f_0) / log(omega_1/omega_0).
void calc_deriv(const vec &inta, vec &deriv) {
  if (verbose) { cout << "Calculating a logarithmic derivative." << endl; }

  int nr = inta.size();
  deriv.resize(nr);
  deriv[0] = 0.0;

  for (int i = 1; i < nr; i++) {
    const double &w0 = mesh[i - 1];
    const double &w1 = mesh[i];
    const double &f0 = inta[i - 1];
    const double &f1 = inta[i];
    assert(w1 > w0);
    const double log1 = log(f1 / f0);
    if (!isfinite(log1)) {
      if (veryverbose) { cout << "Warning: log1 not finite at w0=" << w0 << endl; }
    }

    const double log2 = log(abs(w1 / w0));
    assert(isfinite(log2));

    const double d = log1 / log2;

    if (isfinite(d)) {
      deriv[i] = d;
    } else {
      deriv[i] = 0.0;
    }
  }
}

typedef pair<double, double> Pair;
using Vec = vector<Pair>;

// Linear interpolation class
class LinInt {
  protected:
  Vec vec;               // tabulated data
  int len{};               // length of vec
  int index{};             // index of the interval where last x was found
  double x0{}, x1{};         // last x was in [x0:x1]
  double f0{}, f1{};         // f(x0), f(x1)
  double deriv{};          // (f1-f0)/(x1-x0)
  bool newintegral_flag{}; // set to true when we switch to a new interval
  double xmin{}, xmax{};     // lowest and highest x contained in vec
  double fxmin{}, fxmax{};   // f(xmin), f(xmax)

  public:
  LinInt()= default;;
  LinInt(Vec &in_vec) : vec(in_vec) {
    len              = vec.size();
    index            = -1;
    newintegral_flag = false;
    xmin             = vec.front().first;
    fxmin            = vec.front().second;
    xmax             = vec.back().first;
    fxmax            = vec.back().second;
  };
  void findindex(double x);
  double operator()(double x);
};

// Serach for the interval so that x is contained in [x0:x1].
void LinInt::findindex(double x) {
  if (index == -1) {
    // When interpolation class object is constructed, index is
    // initialized to -1. We have no prior knowledge, thus we search in
    // the full interval.
    for (int i = 0; i < len - 1; i++) {
      x0 = vec[i].first;
      x1 = vec[i + 1].first;
      if (x0 <= x && x <= x1) {
        index = i;
        break;
      }
    }
  } else {
    if (x >= x1) {
      for (int i = index + 1; i < len - 1; i++) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else {
      for (int i = index - 1; i >= 0; i--) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    }
  }

  if (!(0 <= index && index < len && x0 <= x && x <= x1)) {
    cerr << "findindex() error."
         << " x=" << x << " x0=" << x0 << " x1=" << x1 << " index=" << index << " len=" << len << endl;
    exit(1);
  }

  f0             = vec[index].second;     // f(x0)
  f1             = vec[index + 1].second; // f(x1)
  double Delta_y = f1 - f0;
  double Delta_x = x1 - x0;
  deriv          = Delta_y / Delta_x;
}

// Return y(x) using linear interpolation between the tabulated values.
double LinInt::operator()(double x) {
  // Extrapolate if necessary
  if (x <= xmin) { return fxmin; }
  if (x >= xmax) { return fxmax; }

  if (index == -1 || !(x0 <= x && x < x1)) {
    newintegral_flag = true;
    findindex(x);
  }

  double dx = x - x0;
  return f0 + deriv * dx;
}

void recalc_b(const vec &mesh, const vec &bpos, const vec &bneg, vec &b) {
  if (verbose) { cout << "Recalculating b(omega)" << endl; }

  const unsigned int sizepos = bpos.size();
  Vec Vecbpos(sizepos);
  for (unsigned int i = 0; i < sizepos; i++) {
    Vecbpos[i].first  = mesh[i];
    Vecbpos[i].second = bpos[i];
  }
  LinInt fpos(Vecbpos);

  const unsigned int sizeneg = bneg.size();
  Vec Vecbneg(sizeneg);
  for (unsigned int i = 0; i < sizeneg; i++) {
    Vecbneg[i].first  = mesh[i];
    Vecbneg[i].second = bneg[i];
  }
  LinInt fneg(Vecbneg);

  for (unsigned int i = 0; i < nr_spec; i++) {
    const double new_b = (vfreq[i] < 0.0 ? fneg(vfreq[i]) : fpos(vfreq[i]));
    assert(new_b > 0.0);
    b[i] = (new_b > bmin ? new_b : bmin);
  }
}

// Set default values which depend on the chosen values of other
// parameters.
void defaults() {
  if (dyn_mesh > 0.0) {
    // Dynamic mesh
    if (bmin < 0.0) { bmin = 0.0; }
  } else {
    // Fixed mesh
    if (bmin < 0.0) { bmin = broaden_ratio - 1.0; }
    assert(bmin <= b0);
  }

  if (trim < 0.0) { trim = broaden_min * 10; }
}

int main(int argc, char *argv[]) {
  cout << "bw - Adaptive broadening tool - " << VERSION << endl;
  cout << "Rok Zitko, rok.zitko@ijs.si, 2009-2010" << endl;
  cout << setprecision(16);

  cmd_line(argc, argv);
  defaults();
  read_files();
  merge();

  make_mesh(mesh);
  initial_b(b);

  if (savemore) { save("b", vfreq, b, 0, trim); }

  for (int iter = 1; iter <= nr_iter; iter++) {
    broaden(mesh, a, b);

    if (savemore || iter == nr_iter) { save("a", mesh, a, iter, trim); }

    integrate_a(a, mesh, inta);
    if (saveall) { save("inta", mesh, inta, iter); }

    // inta = \int [-infty, omega]
    // intb = \int [0, omega]
    // intc = \int [+infty, omega]
    combinations(mesh, inta, intb, intc);
    if (saveall) {
      save("intb", mesh, intb, iter);
      save("intc", mesh, intc, iter);
    }

    vec deriva, derivb, derivc;

    calc_deriv(inta, deriva);
    calc_deriv(intb, derivb);
    calc_deriv(intc, derivc);

    if (saveall) {
      save("deriva", mesh, deriva, iter);
      save("derivb", mesh, derivb, iter);
      save("derivc", mesh, derivc, iter);
    }

    vec bpos;
    calc_b(derivb, derivc, bpos);

    vec bneg;
    calc_b(deriva, derivb, bneg);

    recalc_b(mesh, bpos, bneg, b);

    if (saveall) {
      save("bpos", mesh, bpos, iter);
      save("bneg", mesh, bneg, iter);
    }
    if (savemore) { save("b", vfreq, b, iter, trim); }

    if (dyn_mesh > 0.0) { refine_mesh(mesh, a); }
  }
}
