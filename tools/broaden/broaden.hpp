// broaden - Finite-temperature raw spectral data broadening tool

#ifndef _broaden_broaden_hpp_
#define _broaden_broaden_hpp_

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <atomic>
#include <chrono>
#include <cstddef>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <cassert>
#include <cfloat>
#include <exception>
#include <utility>
#include <vector>
#include <map>
#include <string>
#include <algorithm>
#include <functional>
#include <ios>
#include <istream>
#include <limits>
#include <mutex>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <numeric>
#include <thread>

#include <unistd.h>
#include <getopt.h>

#include <range/v3/all.hpp>

#include "broadening.hpp"
#include "misc.hpp"
#include "basicio.hpp"

#include "../common/diagnostics.hpp"
#include "../common/gsl_config.hpp"
#include "../common/output_file.hpp"
#include "../common/parallel.hpp"

namespace NRG::Broaden {

const double R_1_SQRT2PI = 1.0 / std::sqrt(2.0 * M_PI); // sqrt not constexpr on all platforms (yet)

using std::abs; // important!!
using std::size_t;

inline std::string tostring(const int i) {
  std::ostringstream S;
  S << i;
  return S.str();
}

inline double gaussian_kernel(const double x, const double y, const double sigma) {
  const auto d = (x - y) / sigma;
  return R_1_SQRT2PI * std::exp(-d * d / 2.0) / sigma;
}

// Derivative of the Fermi-Dirac function: -d/dw f_FD(x-y) with effective temperature sigma.
inline double derfd_kernel(const double x, const double y, const double sigma) {
  const auto d = (x - y) / sigma;
  return 1.0 / ((1.0 + std::cosh(d)) * 2.0 * sigma);
}

template<typename FNC>
void for_each_index(const size_t count, const size_t requested_workers, FNC &&calculate) {
  if (count == 0) return;
  if (requested_workers == 0) throw std::invalid_argument("Worker count must be positive.");
  const auto worker_count = std::min(count, requested_workers);
  if (worker_count == 1) {
    for (size_t index = 0; index < count; ++index) calculate(index);
    return;
  }

  std::atomic<size_t> next_index{0};
  std::atomic<bool> stopped{false};
  std::exception_ptr failure;
  auto failure_index = std::numeric_limits<size_t>::max();
  std::mutex failure_mutex;
  auto run_worker = [&] {
    while (!stopped.load(std::memory_order_relaxed)) {
      const auto index = next_index.fetch_add(1, std::memory_order_relaxed);
      if (index >= count) break;
      try {
        calculate(index);
      } catch (...) {
        {
          const std::lock_guard lock(failure_mutex);
          if (!failure || index < failure_index) {
            failure = std::current_exception();
            failure_index = index;
          }
        }
        stopped.store(true, std::memory_order_relaxed);
      }
    }
  };

  std::vector<std::thread> workers;
  workers.reserve(worker_count);
  try {
    for (size_t worker = 0; worker < worker_count; ++worker) workers.emplace_back(run_worker);
  } catch (...) {
    stopped.store(true, std::memory_order_relaxed);
    for (auto &worker : workers) worker.join();
    throw;
  }
  for (auto &worker : workers) worker.join();
  if (failure) std::rethrow_exception(failure);
}

template<typename S, typename T, typename FNC>
void convolve(const std::vector<S> &mesh, std::vector<T> &a, const double sigma, 
              const size_t jobs, FNC kernel, const double cutoff_ratio = 100) {
  const auto nr_mesh = mesh.size();
  if (nr_mesh == 0) throw std::runtime_error("Output mesh is empty.");
  if (sigma <= 0.0) throw std::invalid_argument("Convolution width must be greater than 0.");
  auto b(a); // source
  for_each_index(nr_mesh, jobs, [&](const size_t i) {
    if (i == 0 || i + 1 == nr_mesh) {
      a[i] = b[i];
      return;
    }
    const auto x = mesh[i];
    if (std::abs(x) < cutoff_ratio * sigma) {
      auto sum = 0.0;
      auto sum_kernel = 0.0;
      for (auto j = 1; j < nr_mesh - 1; j++) {
        const auto y     = mesh[j];
        const auto width = (mesh[j + 1] + mesh[j]) / 2.0 - (mesh[j - 1] + mesh[j]) / 2.0;
        const auto kw    = kernel(x, y, sigma) * width;
        assert(std::isfinite(kw));
        sum += b[j] * kw;
        sum_kernel += kw;
      }
      a[i] = sum;
    } else {
      // High temperatures: convolution not necessary
      a[i] = b[i];
    }
  });
}

template<typename S, typename T>
void save(const std::string &filename, const std::vector<S> &x, const std::vector<T> &y, 
          bool verbose, const int SAVE_PREC = 18) {
  if (verbose) { std::cerr << "Saving " << filename << std::endl; }
  std::ofstream F(filename.c_str());
  if (!F) throw std::runtime_error("Failed to open " + filename + " for writing.");
  F << std::setprecision(SAVE_PREC);
  assert(x.size() == y.size());
  const auto nr = x.size();
  for (auto i = 0; i < nr; i++) F << x[i] << " " << y[i] << '\n';
  NRG::Tools::finish_output(F, filename);
  try {
    F.close();
  } catch (const std::ios_base::failure &error) {
    throw std::runtime_error(filename + ": output close failed: " + error.what());
  }
  if (!F) throw std::runtime_error(filename + ": output close failed.");
}

// Estimate the weight using the trapezoidal rule. The x array must be sorted.
template<typename T> T trapez(const std::vector<T> &x, const std::vector<T> &y) {
  T weight = 0.0;
  assert(x.size() == y.size());
  const auto nr = x.size();
  for (auto i = 1; i < nr; i++) {
    assert(x[i] >= x[i-1]);
    weight += (y[i-1] + y[i]) / 2 * (x[i] - x[i-1]);
  }
  return weight;
}

// Create a mesh on which the output spectral function will be computed. a is the accumulation point of the mesh.
auto make_mesh(const double min, const double max, const double ratio, 
               const double a, const bool add_positive = true, const bool add_negative = false) {
  if (!std::isfinite(ratio) || ratio <= 1.0)
    throw std::invalid_argument("broaden_ratio must be finite and greater than 1.");
  if (!std::isfinite(max) || max <= 0.0)
    throw std::invalid_argument("broaden_max must be finite and greater than 0.");
  if (!std::isfinite(min) || min <= 0.0)
    throw std::invalid_argument("broaden_min must be finite and greater than 0.");
  if (min >= max) throw std::invalid_argument("broaden_min must be smaller than broaden_max.");
  const auto rescale_factor = (max-a)/max;
  std::vector<double> mesh;
  for (double z = max; z > min;) {
    const auto x = a + z * rescale_factor;
    if (add_positive) mesh.push_back(x);
    if (add_negative) mesh.push_back(-x);
    const double next = z / ratio;
    if (!(next < z)) throw std::runtime_error("Frequency mesh generation failed to make progress.");
    z = next;
  }
  std::sort(mesh.begin(), mesh.end());
  return mesh;
}

auto load_mesh(const std::string &filename, bool verbose = false) {
  auto l = readtable<double,double>(filename, false);
  std::vector<double> mesh;
  for (const auto &first : l | ranges::views::transform([](const auto& pair) { return pair.first; }))
    mesh.push_back(first);
  if (verbose) std::cerr << "Reading " << filename << "; " << mesh.size() << " lines read." << std::endl;
  return mesh;
}

auto file_size(std::ifstream &f) {
  f.seekg(0, std::ios::beg);
  const auto begin_pos = f.tellg();
  f.seekg(0, std::ios::end);
  const auto end_pos   = f.tellg();
  return end_pos - begin_pos;
}
   
auto nr_doubles(std::ifstream &f) {
  const auto len = file_size(f);
  if (len % sizeof(double) != 0) throw std::runtime_error("Binary input size is not a multiple of sizeof(double).");
  return len/sizeof(double);
}
  
auto nr_rows(std::ifstream &f, const int cols) {
  const auto d = nr_doubles(f);
  if (d % cols != 0) throw std::runtime_error("Binary input does not contain a whole number of rows.");
  return d / cols;
}
   
// Load a file containing binary representation of raw spectral density. The grid is not assumed to be uniform.
auto load(const std::string &filename, const int nrcol, const bool verbose) {
  std::ifstream f(filename.c_str(), std::ios::in | std::ios::binary);
  if (!f.good() || f.eof() || !f.is_open()) 
    throw std::runtime_error("Error opening file " + filename + " for reading.");
  if (verbose) { std::cerr << "Reading " << filename << std::endl; }
  const auto len = file_size(f); // in bytes
  const auto cols = 1 + nrcol;
  const auto rows = nr_rows(f, cols);
  if (verbose) { std::cerr << "len=" << len << ", " << rows << " data points" << std::endl; }
  auto buffer = std::vector<double>(cols*rows);
  f.seekg(0, std::ios::beg); // Return to the beginning of the file.
  f.read((char *)buffer.data(), len);
  if (f.fail()) throw std::runtime_error("Error reading " + filename);
  f.close();
  return buffer;
}

class Broaden {
 private:
    bool verbose         = false; // output verbosity level
    bool veryverbose     = false;
   bool sumrules        = false; // compute the integrals for testing the T!=0 sum rules
   bool cumulative      = false; // output of integrated spectral function
   double broaden_min   = 1e-7;  // parameters defining the frequency mesh
   double broaden_max   = 2.0;
   double broaden_ratio = 1.01;
   double alpha;          // broadening parameter alpha
    double ggamma = 0.0, dgamma = 0.0; // broadening parameter gamma (final Gaussian, derFD)
   double T;              // temperature parameter
    double omega0_ratio;   // omega0=omega0_ratio*T
    double omega0;
    bool omega0_ratio_defaulted = false;
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
   int nrcol = 1; // Number of columns
   int col   = 1; // Which y column are we interested in?
   std::string name; // filename of binary files containing the raw data
   int Nz;      // Number of spectra (1..Nz)
   const std::string output_filename     = "spec.dat";
   const std::string cumulative_filename = "cumulative.dat";
   std::vector<std::vector<double>> buffers; // binary data buffers
   using mapdd = std::map<double, double>;
   using vec = std::vector<double>;
   mapdd spec;           // Spectrum
   unsigned int nr_spec; // Number of raw spectrum points
   vec vfreq, vspec;     // Same info as spectrum, but in vector<double> form
   vec mesh; // Frequency mesh
   vec a;    // Spectral function
   vec c;    // Cumulative spectrum = int_{-inf}^omega a(x)dx.
   std::string mesh_filename = "";
   size_t jobs = 1;
   std::string jobs_source = "default";
   size_t actual_workers = 0;
   std::chrono::steady_clock::time_point wall_start;
   std::clock_t cpu_start;
   
   void usage(std::ostream &F = std::cout) {
     F << "Usage: broaden <name> <Nz> <alpha> <T> [omega0_ratio]" << std::endl;
     F << std::endl;
     F << "Optional parameters:" << std::endl;
     F << " -h -- show help (when used as sole cmd line switch)" << std::endl;
      F << " -v -- print resolved configuration and enable verbose diagnostics" << std::endl;
      F << " -vv -- enable very verbose diagnostics" << std::endl;
      F << " -V, --version -- show version and exit" << std::endl;
      F << " -j <n>, --jobs <n> -- worker count (default: first OMP_NUM_THREADS value, or 1)" << std::endl;
     F << " -m <min> -- minimal mesh frequency" << std::endl;
     F << " -M <max> -- maximal mesh frequency" << std::endl;
     F << " -r <ratio> -- ratio between two consecutive frequency points" << std::endl;
     F << " -o -- one .dat file" << std::endl;
     F << " -2 -- use the 2nd column for weight values (complex spectra)" << std::endl;
     F << " -3 -- use the 3rd column for weight values (complex spectra)" << std::endl;
     F << " -n -- normalization-conserving broadening kernel" << std::endl;
     F << " -s -- compute weighted integrals for testing sum-rules" << std::endl;
     F << " -c -- compute cumulative spectrum" << std::endl;
     F << " -g -- Gaussian broadening (width alpha)" << std::endl;
     F << " -f -- final Gaussian broadening pass" << std::endl;
     F << " -x -- final derFD broadening pass" << std::endl;
     F << " -a -- accumulation point for the mesh" << std::endl;
     F << " -l -- filter out low-frequency raw data" << std::endl;
     F << " -h -- filter out high-frequency raw data" << std::endl;
     F << " -P -- keep only positive input frequencies" << std::endl;
     F << " -N -- keep only negative input frequencies" << std::endl;
     F << " -A -- output only positive frequencies" << std::endl;
     F << " -B -- output only negative frequencies" << std::endl;
     F << " -L <filename> -- load the frequency mesh from a file" << std::endl;
   }

   void resolve_jobs(const std::optional<size_t> requested) {
     if (requested) {
       jobs = *requested;
       jobs_source = "--jobs";
     } else if (const auto *environment = std::getenv("OMP_NUM_THREADS")) {
       jobs = NRG::Tools::parse_worker_count(environment, "OMP_NUM_THREADS");
       jobs_source = "OMP_NUM_THREADS";
     } else {
       jobs = 1;
       jobs_source = "default";
     }
   }

   void cmd_line(int argc, char *argv[]) {
     if (argc == 2 && std::string(argv[1]) == "-h") {
       usage();
       std::exit(EXIT_SUCCESS);
     }
     static const option long_options[] = {
       {"jobs", required_argument, nullptr, 'j'},
       {nullptr, 0, nullptr, 0}
     };
     std::optional<size_t> requested_jobs;
     int c_;
     while ((c_ = getopt_long(argc, argv, "vj:m:M:r:o23nscgf:x:a:l:h:PNABL:", long_options, nullptr)) != -1) {
       switch (c_) {
       case 'c': cumulative = true; break;
       case 's': sumrules = true; break;
       case 'v':
         if (verbose) veryverbose = true;
         verbose = true;
         break;
       case 'j': requested_jobs = NRG::Tools::parse_positive_size(optarg, "Worker count"); break;
        case 'm':
          broaden_min = atof(optarg);
          break;
        case 'M':
          broaden_max = atof(optarg);
          break;
        case 'r':
          broaden_ratio = atof(optarg);
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
          break;
        case 'x':
          finalderfd = true;
          dgamma     = atof(optarg);
          break;
        case 'a':
          accumulation = atof(optarg);
          break;
        case 'l':
          filterlow = atof(optarg);
          break;
        case 'h':
          filterhigh = atof(optarg);
          break;
       case 'P': keepnegative = false; break;
       case 'N': keeppositive = false; break;
       case 'A': meshnegative = false; break;
       case 'B': meshpositive = false; break;
       case 'L':
         mesh_filename = std::string(optarg);
         break;
       default: throw std::invalid_argument("Unknown command-line option.");
       }
     }
     resolve_jobs(requested_jobs);
     auto remaining = argc - optind; // arguments left
     if (remaining != 5 && remaining != 4) {
       usage();
       std::exit(1);
     }
     name = std::string(argv[optind]); // Name of spectral density files
     Nz = atoi(argv[optind + 1]); // Number of z-values
     if (!(Nz >= 1)) throw std::invalid_argument("Nz must be greater than or equal to 1.");
     alpha = atof(argv[optind + 2]); // High-energy broadening parameter
     if (!(alpha > 0.0)) throw std::invalid_argument("alpha must be greater than 0.");
     T = atof(argv[optind + 3]); // Temperature
     if (!(T > 0.0)) throw std::invalid_argument("T must be greater than 0.");
     if (remaining == 5) {
       omega0_ratio = atof(argv[optind + 4]); // omega0/T
       if (!(omega0_ratio > 0.0)) throw std::invalid_argument("omega0_ratio must be greater than 0.");
       omega0 = omega0_ratio * T;
     }
      if (remaining == 4) {
        omega0_ratio = 1e-9; // Effectively zero
        omega0       = omega0_ratio * T;
        omega0_ratio_defaulted = true;
     }
     if (finalgaussian && ggamma <= 0.0) throw std::invalid_argument("Final Gaussian width must be greater than 0.");
      if (finalderfd && dgamma <= 0.0) throw std::invalid_argument("Final derFD width must be greater than 0.");
      if (!meshpositive && !meshnegative) throw std::invalid_argument("Output mesh cannot be empty.");
      std::cout << "Processing: " << name << std::endl;
     }

    auto kernel_mode_name() const {
      return gaussian ? "gaussian" : normalization ? "hybrid_peak_frequency" : "hybrid_output_frequency";
    }

    void resolve_mesh() {
      if (mesh_filename != "") {
        mesh = load_mesh(mesh_filename, verbose);
      } else {
        mesh = make_mesh(broaden_min, broaden_max, broaden_ratio, accumulation, meshpositive, meshnegative);
      }
      if (mesh.empty()) throw std::runtime_error("Output mesh is empty.");
    }

    void report_configuration() const {
      NRG::Tools::ConfigurationReport report("broaden");
      report.value("verbosity", veryverbose ? 2 : 1);
      report.value("input_name", name);
      report.value("Nz", Nz);
      report.value("one_file_mode", one);
      report.value("input_layout", one && Nz == 1 ? "direct" : "numbered_subdirectories");
      report.value("input_columns", 1 + nrcol);
      report.value("weight_column", col + 1);
      report.value("alpha", alpha);
      report.value("temperature", T);
      if (omega0_ratio_defaulted) {
        report.resolved("omega0_ratio", omega0_ratio, "omitted positional argument");
      } else {
        report.value("omega0_ratio", omega0_ratio);
      }
      report.resolved("omega0", omega0, "omega0_ratio*temperature");
      report.value("kernel_mode", kernel_mode_name());
      if (jobs_source == "--jobs")
        report.value("jobs", jobs);
      else
        report.resolved("jobs", jobs, jobs_source);
      report.value("jobs.source", jobs_source);
      report.value("normalization_mode", normalization);
      report.value("sum_rules", sumrules);
      report.value("cumulative", cumulative);
      report.value("filter_low", filterlow);
      report.value("filter_high", filterhigh);
      report.value("keep_positive_input", keeppositive);
      report.value("keep_negative_input", keepnegative);
      report.value("accumulation", accumulation);
      report.value("final_gaussian", finalgaussian);
      if (finalgaussian) {
        report.value("final_gaussian_gamma", ggamma);
        report.resolved("final_gaussian_sigma", ggamma * T, "gamma*temperature");
      }
      report.value("final_derfd", finalderfd);
      if (finalderfd) {
        report.value("final_derfd_gamma", dgamma);
        report.resolved("final_derfd_sigma", dgamma * T, "gamma*temperature");
      }
      report.value("mesh_source", mesh_filename.empty() ? "generated" : "loaded");
      if (mesh_filename.empty()) {
        report.value("broaden_min", broaden_min);
        report.value("broaden_max", broaden_max);
        report.value("broaden_ratio", broaden_ratio);
        report.value("mesh_positive", meshpositive);
        report.value("mesh_negative", meshnegative);
      } else {
        report.value("mesh_file", mesh_filename);
      }
      const auto [mesh_min, mesh_max] = std::minmax_element(mesh.begin(), mesh.end());
      report.value("mesh_points", mesh.size());
      report.value("mesh_min", *mesh_min);
      report.value("mesh_max", *mesh_max);
      report.value("mesh_sorted", std::is_sorted(mesh.begin(), mesh.end()));
      report.value("output_file", output_filename);
      report.value("output_precision", 18);
      report.value("cumulative_file", cumulative ? cumulative_filename : "disabled");
      report.write(std::cerr);
    }

     void report_timing() const {
       const auto wall_seconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - wall_start).count();
       const auto cpu_end = std::clock();
       std::optional<double> cpu_seconds;
       if (cpu_start != static_cast<std::clock_t>(-1) && cpu_end != static_cast<std::clock_t>(-1)
           && cpu_end >= cpu_start)
         cpu_seconds = static_cast<double>(cpu_end - cpu_start) / static_cast<double>(CLOCKS_PER_SEC);

       std::ostringstream report;
       report << std::setprecision(6) << "Time elapsed: " << wall_seconds << " s\n";
       if (verbose) {
         const auto throughput = wall_seconds > 0.0 ? static_cast<double>(mesh.size()) / wall_seconds : 0.0;
         const auto passes = 1 + static_cast<int>(finalgaussian) + static_cast<int>(finalderfd);
         report << "Performance: wall=" << wall_seconds << "s";
         if (cpu_seconds && wall_seconds > 0.0)
           report << " cpu=" << *cpu_seconds << "s effective_parallelism=" << *cpu_seconds / wall_seconds << "x";
         else
           report << " cpu=n/a effective_parallelism=n/a";
         report << " throughput=" << throughput << " points/s workers=" << actual_workers
                << " passes=" << passes << " kernel=" << kernel_mode_name() << '\n';
       }
       std::cout << report.str();
       NRG::Tools::finish_output(std::cout, "<stdout>");
     }

    void check_buffer_normalisation(const std::vector<double> &buffer, const int col_) {
     const auto cols = 1 + nrcol;
     const auto rows = buffer.size()/cols;
     auto sum = 0.0;
     for (auto j = 0; j < rows; j++) sum += buffer[cols * j + col_];
      std::cerr << "Weight=" << sum << std::endl;
   }
   // Load all the input data.
   void read_files() {
     buffers.resize(Nz+1);
     for (auto i = 1; i <= Nz; i++) {
       buffers[i] = load(one && Nz == 1 ? name : tostring(i) + "/" + name, nrcol, verbose);
       if (verbose) check_buffer_normalisation(buffers[i], col);
     }
   }
   // Combine data from all NRG runs (z-averaging).
   void merge() {
     const auto cols = 1 + nrcol; // number of elements in a line
     // Sum weight corresponding to equal frequencies.  Map of (frequency,weight) pairs is used for this purpose.
     for (auto i = 1; i <= Nz; i++) {
       for (auto l = 0; l < buffers[i].size()/cols; l++) {
         auto freq  = buffers[i][cols * l];
         auto value = buffers[i][cols * l + col];
         spec[freq] += value;
       }
     }
     nr_spec = spec.size();
      if (verbose) { std::cerr << nr_spec << " unique frequencies." << std::endl; }
     // Normalize weight by 1/Nz, determine total weight, and store the (frequency,weight) data in the form of linear
     // vectors for faster access in the ensuing calculations.
     auto sum = 0.0;
     for (auto & I : spec) {
       const auto weight = (I.second /= Nz); // Normalize weight on the fly
       const auto freq   = I.first;
       vfreq.push_back(freq);
       vspec.push_back(weight);
       sum += weight;
     }
      if (verbose) { std::cerr << "Total weight=" << sum << std::endl; }
     assert(vfreq.size() == nr_spec && vspec.size() == nr_spec);
   }

   void filter() {
     for (auto j = 0; j < nr_spec; j++) {
       if (std::abs(vfreq[j]) < filterlow) { vspec[j] = 0.0; }
       if (std::abs(vfreq[j]) > filterhigh) { vspec[j] = 0.0; }
     }
     for (auto j = 0; j < nr_spec; j++) {
       if (!keeppositive && vfreq[j] >= 0) { vspec[j] = 0.0; }
       if (!keepnegative && vfreq[j] < 0) { vspec[j] = 0.0; }
     }
   }

   void integrals_for_sumrules() {
     const auto nr    = vfreq.size();
     auto sum         = 0.0;
     auto sumpos      = 0.0;
     auto sumneg      = 0.0;
     auto sumfermi    = 0.0;
     auto sumbose     = 0.0;
     auto sumfermiinv = 0.0;
     auto sumboseinv  = 0.0;
     for (int i = 0; i < nr; i++) {
       const auto omega = vfreq[i];
       const auto w     = vspec[i];
       sum += w;
       if (omega > 0.0) { sumpos += w; }
       if (omega < 0.0) { sumneg += w; }
       const auto f  = 1 / (1 + std::exp(-omega / T));
       const auto b  = 1 / (1 - std::exp(-omega / T));
       const auto fi = 1 / (1 + std::exp(+omega / T)); // fi=1-f
       const auto bi = 1 / (1 - std::exp(+omega / T)); // bi=1-b
       if (std::isfinite(f))  { sumfermi    += f * w;  }
        if (std::isfinite(b))  { sumbose     += b * w;  }
       if (std::isfinite(fi)) { sumfermiinv += fi * w; }
       if (std::isfinite(bi)) { sumboseinv  += bi * w; }
     }
     std::cout << "Total weight=" << sum << std::endl;
     std::cout << "Positive-omega weight=" << sumpos << std::endl;
     std::cout << "Negative-omega weight=" << sumneg << std::endl;
     std::cout << "Integral with fermionic kernel=" << sumfermi << std::endl;
     std::cout << "Integral with bosonic kernel=" << sumbose << std::endl;
     std::cout << "Integral with fermionic kernel (omega -> -omega)=" << sumfermiinv << std::endl;
     std::cout << "Integral with bosonic kernel (omega -> -omega)=" << sumboseinv << std::endl;
   }

    // e - output energy
    // ept - energy of the delta peak (data point)
    auto bfnc(const double e, const double ept) const {
      if (gaussian) {
        return NRG::Broadening::gaussian(e, ept, alpha);
      } else {
        const auto mode = normalization ? NRG::Broadening::CrossoverMode::peak_frequency
                                        : NRG::Broadening::CrossoverMode::output_frequency;
        return NRG::Broadening::hybrid_kernel(e, ept, alpha, omega0, mode, accumulation);
      }
    }

   vec broaden(const vec &mesh_) {
      if (verbose) { std::cerr << "Broadening. Number of mesh points = " << mesh_.size() << std::endl; }
     vec result(mesh_.size());
     for_each_index(mesh_.size(), jobs, [&](const size_t index) {
       const auto m = mesh_[index];
       result[index] = std::transform_reduce(vspec.begin(), vspec.end(), vfreq.begin(), 0.0, std::plus<>(),
                                             [this,m](const auto weight, const auto freq) {
                                               return weight * bfnc(m, freq);
                                             });
     });
     return result;
   }
   
   // Cumulative spectrum
    void calc_cumulative(const vec &mesh_, vec &c_) {
      const auto nr_mesh = mesh_.size();
       if (verbose) std::cerr << "Calculating cumulative spectrum." << std::endl;
      c_.resize(nr_mesh);
      auto sum = 0.0;
      size_t j = 0;
      for (auto i = 0; i < nr_mesh; i++) {
        const auto max_freq = mesh_[i];
        while (j < vfreq.size() && vfreq[j] < max_freq) {
          sum += vspec[j];
          j++;
        }
       c_[i] = sum;
     }
     std::cout << "End sum=" << sum << std::endl;
   }

   void check_normalizations(const vec &x) {
     vec y(x.size());
     // For a range of delta peak positions, compute the total weight of the resulting broadened spectra.
     for (auto z = 1e-10; z < 1.0; z *= 10) {
       const auto nr = x.size();
       for (auto i = 0; i < nr; i++) { y[i] = bfnc(x[i], z); }
        std::cerr << "z=" << z << " weight=" << trapez(x, y) << std::endl;
     }
   }
 public:
    Broaden(int argc, char *argv[])
        : wall_start{std::chrono::steady_clock::now()}, cpu_start{std::clock()} {
      std::cout << "broaden - finite-temperature broadening tool" << std::endl;
      cmd_line(argc, argv);
    }
    void calc() {
      read_files();
      merge();
      filter();
      if (sumrules) integrals_for_sumrules();
      resolve_mesh();
      actual_workers = std::min(jobs, mesh.size());
      if (verbose) report_configuration();
      if (verbose) check_normalizations(mesh);
      a = broaden(mesh);
     if (finalgaussian) convolve(mesh, a, ggamma * T, jobs, gaussian_kernel);
     if (finalderfd) convolve(mesh, a, dgamma * T, jobs, derfd_kernel);
     std::cout << "Estimated weight (trapezoidal rule)=" << trapez(mesh, a) << std::endl;
     save(output_filename, mesh, a, verbose);
     if (cumulative) {
       calc_cumulative(mesh, c);
       save(cumulative_filename, mesh, c, verbose);
     }
     report_timing();
    }
};

} // namespace

#endif
