// Numerical renormalization group
// (c) Rok Zitko, rok.zitko@ijs.si, 2005-2019

#ifndef _nrg_h_
#define _nrg_h_

#include <fstream>
#include <sstream>
#include <iostream>
#include <iomanip>

#include <complex>
#include <unordered_map>
#include <map>
#include <list>
#include <deque>
#include <set>

#include <numeric>
#include <math.h> // isfinite
#include <utility>
#include <algorithm>
#include <functional>
#include <limits>

#include <cassert>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cstring>

#include <unistd.h>
#include <cstdlib> // mkdtemp

#include "portabil.h"

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <boost/lexical_cast.hpp>
#include <boost/io/ios_state.hpp>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/serialization/vector.hpp>

//#include <boost/serialization/map.hpp>
//#include <boost/serialization/list.hpp>
//#include <boost/serialization/complex.hpp>
//#include <boost/serialization/base_object.hpp>
//#include <boost/serialization/utility.hpp>
//#include <boost/serialization/list.hpp>

using namespace std;
using namespace boost::numeric;
using namespace boost::numeric::ublas;

// Numeric bindings to BLAS/LAPACK
#include <boost/numeric/bindings/traits/ublas_vector.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/atlas/cblas.hpp>

namespace atlas = boost::numeric::bindings::atlas;

// This must be set manually through CXXFLAG=-DNRG_MPI
#ifdef NRG_MPI
#include <boost/optional.hpp>
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/collectives.hpp>
namespace mpi = boost::mpi;
mpi::environment *mpienv;
mpi::communicator *mpiw;
int myrank;
const int TAG_HELLO = 1;
const int TAG_EXIT = 2;
const int TAG_DIAG = 3;
const int TAG_SYNC = 4;
const int TAG_MATRIX = 5;
const int TAG_INVAR = 6;
const int TAG_EIGEN = 7;
const int TAG_MATRIX_SIZE = 8;
const int TAG_MATRIX_LINE = 9;
const int TAG_EIGEN_INT = 10;
const int TAG_EIGEN_VEC = 11;
const int TAG_EIGEN_RMAXVALS = 12;
#else
int myrank = 0; // in order to simplify the code
#endif // NRG_MPI

// Support for compiler dependant optimizations

#ifdef __GNUC__
#define PUREFNC __attribute__ ((pure))
#define CONSTFNC __attribute__ ((const))
#else
#define PUREFNC
#define CONSTFNC
#endif

// Hard coded & default filenames

// Spectral densities and temperature dependant quantities
#define FN_ENERGIES_NRG "energies.nrg"
#define FN_ENERGIES_DMNRG "energies.dmnrg"
#define FN_ANNOTATED "annotated.dat"

#define FN_TD "td"
#define FN_CUSTOM "custom"
#define FN_CUSTOMFDM "customfdm"
#define FN_CUSTOMSQ "customsq"

// Density-matrix approach
#define FN_UNITARY "unitary"
#define FN_RHO "rho"
#define FN_RHOFDM "rhofdm"

#define FN_SUBSPACES "subspaces.dat"

#define T_WIDTH 12

// Warning: not thread safe!
class Timing {
private:
  double all, timer;
  bool running;
  map<string, double> t;
  bool reportonexit;

public:  
  Timing(bool _reportonexit = true) : reportonexit(_reportonexit) {
     all = gettime();
     timer = gettime();
     running = false;
  }
  void start(void) {
    my_assert(!running);
    running = true;
    timer = gettime();
  }
  double stop(void) {
    my_assert(running);
    running = false;
    double end = gettime();
    return end-timer;
  }
  void add(string timer) {
    if (!t.count(timer)) 
       t[timer] = 0.0; // initialize
    t[timer] += stop();
  }
  double value(string timer) const {
    if (!t.count(timer)) 
       return 0.0;
    return t.find(timer)->second;
  }
  double total() {
     double end = gettime();
     return end-all;
  }   
  void report(void) {
    const double t_all = total();
    cout << endl;
    cout << "Timing report (process " << myrank << ")" << endl;
    cout << "=============" << endl;
    cout << setw(T_WIDTH) << "All" << ": " << long(t_all) << " s" << endl;
    double t_sum = 0.0;
    for (const auto &i : t) {
       // Only show those that contribute more than 1% of the total time!
       if ( (i.second/t_all) > 0.01) {
	  cout << setw(T_WIDTH) << i.first << ": " << long(i.second) << " s" << endl;
	  if (i.first[0] != '*') 
	     t_sum += i.second;
       }
    }
     cout << setw(T_WIDTH) << "Other" << ": " << long(t_all-t_sum) << " s" << endl;
  }
  ~Timing() { if (reportonexit) report(); }
};

// Higher-level timing code

// Time a section for as long as the object is in scope.
class TimeScope 
{
 private:
   Timing & timer;
   string timer_name;
 public:
   TimeScope(Timing & _timer, const string &_timer_name) :
    timer(_timer), timer_name(_timer_name) {
       timer.start();
    }
   ~TimeScope() {
      timer.add(timer_name);
   }
};

#define TIME(timer_name) TimeScope timer(t, timer_name)
#define TIME_SECTION(timer_name, section)  { TIME(timer_name); section }

// Stores maximal memory usage at various breakpoints.
// This is useful for estimating memory requirements at various
// points of the execution path.

#define MS_WIDTH 12

class MemoryStats {
private:
  map<string, int> maxvals;
  bool reportonexit;
  int peakusage;
public:
  MemoryStats(bool rpexit = true) {
    reportonexit = rpexit;
    peakusage = 0;
  }
  // Intermediate level routine.
  int used() {
     const int memused = memoryused();
     peakusage = max(peakusage, memused);
     return memused; 
  }
  // Sample memory usage at an arbitrarily named "breakpoint".
  int check(string breakpoint) {
    const int memused = used();
    if (maxvals.count(breakpoint) == 0) 
      maxvals[breakpoint] = 0;
    maxvals[breakpoint] = max(maxvals[breakpoint], memused);
    return memused;
  }
  void report(ostream &F = cout) const {
#ifdef HAS_MEMORY_USAGE     
     F << endl;
     F << "Memory usage report (process " << myrank << ")" << endl;
     F << "===================" << endl;
     int topusage = 0; // top usage recorded by check()
     for (const auto &i : maxvals)
	topusage = max(topusage, i.second);
     if (topusage != 0) {
	for (const auto &i : maxvals)
	   F << setw(MS_WIDTH) << i.first << ": " << i.second << " kB" << endl;
     }
     my_assert(topusage <= peakusage);
     F << endl << "Peak usage: " << peakusage << " kB" << endl;
#endif
  }   
  ~MemoryStats() { if (reportonexit) report(); }
};

#define HIGHPREC(val) setw(30) << setprecision(16) << (val) \
                                 << setprecision(COUT_PRECISION)  

#endif
