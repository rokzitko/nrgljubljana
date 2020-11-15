#ifndef _binavg_binavg_hpp_
#define _binavg_binavg_hpp_

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
#include <unistd.h> // getopt

namespace NRG::BinAvg {

inline std::string tostring(const int i) {
  std::ostringstream S;
  S << i;
  return S.str();
}

class BinAvg {
 private:   
   bool verbose = false; // output verbosity level
   bool one     = false; // For Nz=1, no subdir.
   int nrcol = 1; // Number of columns
   int col   = 1; // Which y column are we interested in?
   std::string name; // filename of binary files containing the raw data
   int Nz;           // Number of spectra (1..Nz)
   double **buffers; // binary data buffers
   int *sizes;       // sizes of buffers

   typedef std::map<double, double> mapdd;
   using vec = std::vector<double>;
   mapdd spec;           // Spectrum
   unsigned int nr_spec; // Number of raw spectrum points
   vec vfreq, vspec;     // Same info as spectrum, but in vector<double> form

   void usage(std::ostream &F = std::cout) {
     F << "Usage: binavg <name> <Nz>" << std::endl;
     F << std::endl;
     F << "Optional parameters:" << std::endl;
     F << " -h -- show help" << std::endl;
     F << " -v -- verbose" << std::endl;
     F << " -o -- one .dat file" << std::endl;
     F << " -2 -- use the 2nd column for weight values (complex spectra)" << std::endl;
     F << " -3 -- use the 3rd column for weight values (complex spectra)" << std::endl;
   }

   void cmd_line(int argc, char *argv[]) {
     char c;
     while (c = getopt(argc, argv, "hvo23"), c != -1) {
       switch (c) {
       case 'h':
         usage();
         exit(EXIT_SUCCESS);
       case 'v': verbose = true; break;
       case 'o': one = true; break;
       case '2':
         nrcol = 2;
         col   = 1;
         break;
       case '3':
         nrcol = 2;
         col   = 2;
         break;
       default: abort();
       }
     }
     int remaining = argc - optind; // arguments left
     if (remaining != 2) {
       usage();
       exit(1);
     }
     name = std::string(argv[optind]); // Name of spectral density files
     Nz = atoi(argv[optind + 1]); // Number of z-values
     assert(Nz >= 1);
     std::cout << "Processing: " << name << std::endl;
     std::cout << "Nz=" << Nz << std::endl;
   }
   
   // Load a file containing binary representation of raw spectral density. The grid is not assumed to be uniform.
   void load(const int i) {
     // i-th z-value defines the name of the directory where the results of the NRG calculation are contained.
     std::string filename = one && Nz == 1 ? name : tostring(i) + "/" + name;
     std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
     if (!f.good() || f.eof() || !f.is_open()) {
       std::cerr << "Error opening file " << filename << std::endl;
       exit(1);
     }
     if (verbose) { std::cout << "Reading " << filename << std::endl; }
     const int rows = 1 + nrcol; // number of elements in a line
     // Determine the number of records
     f.seekg(0, std::ios::beg);
     const auto begin_pos = f.tellg();
     f.seekg(0, std::ios::end);
     const auto end_pos   = f.tellg();
     const long len       = end_pos - begin_pos;
     assert(len % (rows * sizeof(double)) == 0);
     const int nr = len / (rows * sizeof(double)); // number of lines
     if (verbose) { std::cout << "len=" << len << " nr=" << nr << " data points" << std::endl; }
     // Allocate the read buffer. The data will be kept in memory for the duration of the calculation!
     auto *buffer = new double[rows * nr];
     f.seekg(0, std::ios::beg); // Return to the beginning of the file.
     f.read((char *)buffer, len);
     if (f.fail()) {
       std::cerr << "Error reading " << filename << std::endl;
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
       std::cout << "Weight=" << sum << std::endl;
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
     // Sum weight corresponding to equal frequencies.  Map of (frequency,weight) pairs is used for this purpose.
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
     if (verbose) { std::cout << nr_spec << " unique frequencies." << std::endl; }
     // Normalize weight by 1/Nz, determine total weight, and store the (frequency,weight) data in the form of linear
     // vectors for faster access in the ensuing calculations.
     double sum = 0.0;
     for (auto & I : spec) {
       const double weight = (I.second /= Nz); // Normalize weight on the fly
       const double freq   = I.first;
       vfreq.push_back(freq);
       vspec.push_back(weight);
       sum += weight;
     }
     if (verbose) { std::cout << "Total weight=" << sum << std::endl; }
     assert(vfreq.size() == nr_spec && vspec.size() == nr_spec);
   }

   // Save a map of (double,double) pairs to a BINARY file.
   void save_binary(const std::string filename, const vec &x, const vec &y, const int SAVE_PREC = 18) {
     if (verbose) { std::cout << "Saving " << filename << std::endl; }
     std::ofstream F(filename.c_str(), std::ios::out | std::ios::binary);
     if (!F) {
       std::cerr << "Failed to open " << filename << " for writing." << std::endl;
       exit(1);
     }
     F << std::setprecision(SAVE_PREC);
     assert(x.size() == y.size());
     unsigned int nr = x.size();
     for (unsigned int i = 0; i < nr; i++) {
       const double vx = x[i];
       F.write((char *)&vx, sizeof(double));
       const double vy = y[i];
       F.write((char *)&vy, sizeof(double));
     }
   }
   
 public:
   BinAvg(int argc, char *argv[]) {
     cmd_line(argc, argv);
   }
   void calc() {
     read_files();
     merge();
     save_binary(one ? "spec.bin" : name, vfreq, vspec);
   }
};

} // namespace

#endif
