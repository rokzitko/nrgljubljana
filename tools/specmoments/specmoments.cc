// specmoments - Compute spectral moments from raw spectal data
// Ljubljana code, Rok Zitko, rok.zitko@ijs.si, Aug 2015

// CHANGE LOG
// 25.8.2015 - first version

#define VERSION "1.0"

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

bool verbose = false; // output verbosity
bool one = false; // For Nz=1, no subdirs.
double filterlow = 0.0; // filter all input data points with |omega|<filterlow
double filterhigh = DBL_MAX; // filter all input data points with |omega|>filterhigh
bool posonly = false; // keep only positive omega
bool negonly = false; // keep only negative omega
bool central = false; // central or raw moments?

// Multi-column support
int nrcol = 1; // Number of columns
int col = 1; // Which y column are we interested in?

string name;
int Nz;

const int COUT_PREC = 16;

void about(ostream &F = cout)
{
  F << "specmoments - Spectral moment calculation tool - " << VERSION << endl;
  F << "Rok Zitko, rok.zitko@ijs.si, 2015" << endl;
}

void usage(ostream &F = cout)
{
  F << "Usage: specmoments <name> <Nz>\n";
  F << endl;
  F << "Optional parameters:" << endl;
  F << " -v -- verbose" << endl;
  F << " -o -- one .dat file" << endl;
  F << " -l -- filter out low-frequency raw data" << endl;
  F << " -h -- filter out high-frequency raw data" << endl;
  F << " -p -- positive frequencies only" << endl;
  F << " -n -- negative frequencies only" << endl;
  F << " -c -- central moments" << endl;
  F << endl;
  F << "Output: <mom0> <mom1> <mom2> <mom3> <mom4>" << endl;
}

void cmd_line(int argc, char *argv[])
{
   char c;
   while ((c = getopt(argc, argv, "vol:h:pnc")) != -1) {
     switch (c) {
     case 'v':
       verbose = true;
       about(); // only if verbose is on
       break;
     case 'o':
       one = true;
       break;
     case 'l':
       filterlow = atof(optarg);
       if (verbose) {
	 cout << "filterlow=" << filterlow << endl;
       }
       break;
     case 'h':
       filterhigh = atof(optarg);
       if (verbose) {
	 cout << "filterhigh=" << filterhigh << endl;
       }
       break;
     case 'p':
       posonly = true;
       if (verbose) {
	 cout << "positive part" << endl;
       }
       break;
     case 'n':
       negonly = true;
       if (verbose) {
	 cout << "negative part" << endl;
       }
       break;
     case 'c':
       central = true;
       break;
     default:
       abort();
     }
   }
   int remaining = argc-optind;
   if (remaining != 2) {
     if (!verbose)
       about();
     usage();
     exit(1);
   }
   name = string(argv[optind]);
   Nz = atoi(argv[optind+1]);
   assert(Nz >= 1);
   if (verbose)
     cout << "Processing: " << name << " Nz=" << Nz << endl;
}

double **buffers; // binary data buffers
int *sizes; // sizes of buffers

typedef map<double, double> mapdd;
typedef vector<double> vec;

mapdd spec; // Spectrum
unsigned int nr_spec; // Number of raw spectrum points
vec vfreq, vspec; // Same info as spectrum, but in vector<double> form

string tostring(int i)
{
  ostringstream S;
  S << i;
  return S.str();
}

// Load a file containing binary representation of raw spectral density.
// The grid is not assumed to be uniform.
void load(int i)
{
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
  if (verbose)
    cout << "Reading " << filename << endl;
  const int rows = 1+nrcol; // number of elements in a line
  // Determine the number of records
  f.seekg(0, ios::beg);
  const ios::pos_type begin_pos = f.tellg();
  f.seekg(0, ios::end);
  const ios::pos_type end_pos = f.tellg();
  const long len = end_pos-begin_pos;
  assert(len % (rows*sizeof(double)) == 0);
  const int nr = len / (rows*sizeof(double)); // number of lines
  if (verbose) 
    cout << "len=" << len
	 << " nr=" << nr << " data points" << endl;
  // Allocate the read buffer. The data will be kept in memory for the
  // duration of the calculation!
  double *buffer = new double[rows*nr];
  f.seekg(0, ios::beg); // Return to the beginning of the file.
  f.read((char*)buffer, len);
  if (f.fail()) {
    cerr << "Error reading " << filename << endl;
    exit(1);
  }
  f.close();
  // Keep record of the the buffer and its size.
  buffers[i] = buffer;
  sizes[i] = nr;
  if (verbose) {
    // Check normalization.
    double sum = 0.0;
    for (int j = 0; j < nr; j++)
      sum += buffer[rows*j+col];
    cout << "Weight=" << sum << endl;
  }
}

// Load all the input data.
void read_files()
{
  buffers = new double*[Nz+1];
  sizes = new int[Nz+1];
  for (int i = 1; i <= Nz; i++)
    load(i);
}

// Combine data from all NRG runs (z-averaging).
void merge()
{
  const int rows = 1+nrcol; // number of elements in a line
  // Sum weight corresponding to equal frequencies.  Map of
  // (frequency,weight) pairs is used for this purpose.
  for (int i = 1; i <= Nz; i++) {
    for (int l = 0; l < sizes[i]; l++) {
      double & freq = buffers[i][rows*l];
      double & value = buffers[i][rows*l+col];
      mapdd::iterator I = spec.find(freq);
      if (I == spec.end()) {
	spec[freq] = value;
      } else {
	I->second += value;
      }
    }
  }
  nr_spec = spec.size();
  if (verbose) {
    cout << nr_spec << " unique frequencies." << endl;
  }
  // Normalize weight by 1/Nz, determine total weight, and store the
  // (frequency,weight) data in the form of linear vectors for faster
  // access in the ensuing calculations.
  double sum = 0.0;
  for (mapdd::iterator I = spec.begin(); I != spec.end(); I++) {
    const double weight = (I->second /= Nz); // Normalize weight on the fly
    const double freq = I->first;
    vfreq.push_back(freq);
    vspec.push_back(weight);
    sum += weight;
  }
  if (verbose) {
    cout << "Total weight=" << sum << endl;
  }
  assert(vfreq.size() == nr_spec && vspec.size() == nr_spec);
}

void filter()
{
  for (unsigned int j = 0; j < nr_spec; j++) {
    if (abs(vfreq[j]) < filterlow)
      vspec[j] = 0.0;
    if (abs(vfreq[j]) > filterhigh)
      vspec[j] = 0.0;
    if (vfreq[j] < 0.0 && posonly)
      vspec[j] = 0.0;
    if (vfreq[j] > 0.0 && negonly)
      vspec[j] = 0.0;
  }
}

// Moments about 0
void moments()
{
   double mom0 = 0.0, mom1 = 0.0, mom2 = 0.0, 
     mom3 = 0.0, mom4 = 0.0;
   for (unsigned int j = 0; j < nr_spec; j++) {
      mom0 += vspec[j];
      mom1 += vspec[j] * vfreq[j];
      mom2 += vspec[j] * pow(vfreq[j], 2);
      mom3 += vspec[j] * pow(vfreq[j], 3);
      mom4 += vspec[j] * pow(vfreq[j], 4);
   }
   // Important: no normalization!
   if (verbose) {
      cout << endl;
      cout << "0. raw moment = " << mom0 << endl;
      cout << "1. raw moment = " << mom1 << endl;
      cout << "2. raw moment = " << mom2 << endl;
      cout << "3. raw moment = " << mom3 << endl;
      cout << "4. raw moment = " << mom4 << endl;
   } else {
      cout << mom0 << " " << mom1 << " " << mom2 << " " << mom3 << " " << mom4;
      cout << endl;
   }   
}

// Moments about the mean
void central_moments()
{
   double mom0 = 0.0, mom1 = 0.0, mom2 = 0.0, 
     mom3 = 0.0, mom4 = 0.0;   
   for (unsigned int j = 0; j < nr_spec; j++) {
      mom0 += vspec[j];
      mom1 += vspec[j] * vfreq[j];
      mom2 += vspec[j] * pow(vfreq[j], 2);
      mom3 += vspec[j] * pow(vfreq[j], 3);
      mom4 += vspec[j] * pow(vfreq[j], 4);
   }
   if (verbose) {
     cout << endl;
     cout << "0. raw moment = " << mom0 << endl;
     cout << "1. raw moment = " << mom1 << endl;
     cout << "2. raw moment = " << mom2 << endl;
     cout << "3. raw moment = " << mom3 << endl;
     cout << "4. raw moment = " << mom4 << endl;
   }
   const double weight = mom0;
   const double mean = mom1/weight;
   double cmom2 = 0.0, cmom3 = 0.0, cmom4 = 0.0;
   for (unsigned int j = 0; j < nr_spec; j++) {
      cmom2 += vspec[j] * pow(vfreq[j]-mean, 2);
      cmom3 += vspec[j] * pow(vfreq[j]-mean, 3);
      cmom4 += vspec[j] * pow(vfreq[j]-mean, 4);
   }
   // Normalize by total weight
   cmom2 /= weight;
   cmom3 /= weight;
   cmom4 /= weight;
   const double sigma = sqrt(cmom2); // std. deviation
   const double skewness = cmom3/pow(sigma, 3);
   const double kurtosis = cmom4/pow(sigma, 4) - 3.0;
   if (verbose) {
      cout << endl;
      cout << "0. moment (weight) = " << mom0 << endl;
      cout << "1. moment          = " << mom1 << endl;
      cout << endl;
      cout << "2. central moment  = " << cmom2 << endl;
      cout << "3. central moment  = " << cmom3 << endl;
      cout << "4. central moment  = " << cmom4 << endl;
      cout << endl;
      cout << "mean               = " << mean << endl;
      cout << "sigma (std dev)    = " << sigma << endl;
      cout << "skewness           = " << skewness << endl;
      cout << "kurtosis           = " << kurtosis << endl;
   } else {
      cout << mom0 << " " << mom1 << " " << cmom2 << " " << cmom3 << " " << cmom4;
      cout << endl;
   }  
}
	 
int main(int argc, char *argv[])
{
  cmd_line(argc, argv);
  cout << setprecision(COUT_PREC);
  read_files();
  merge();
  filter();
  if (central)
    central_moments();
  else
    moments();
}
