/*
 "NRG Ljubljana" - Numerical renormalization group for multiple
 impurities and an arbitrary number of channels

 Copyright (C) 2005-2020 Rok Zitko

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

   Contact information:
   Rok Zitko
   F1 - Theoretical physics
   "Jozef Stefan" Institute
   Jamova 39
   SI-1000 Ljubljana
   Slovenia

   rok.zitko@ijs.si
*/

#include "nrg-general.h"
#include "nrg-lib.h" // exposed interfaces for wrapping into a library
#include "portabil.h"
#include "time_mem.h"
#include "debug.h"
#include "misc.h"
#include "openmp.h"

#include "param.cc"
#include "outfield.cc"

// This is included in the library only. Should not be used if cblas library is available.
#ifdef CBLAS_WORKAROUND
 #define ADD_
 #include "cblas_globals.c"
 #include "cblas_dgemm.c"
 #ifdef NRG_COMPLEX
  #include "cblas_zgemm.c"
 #endif
 #include "cblas_xerbla.c"
#endif

// Timing of various parts of the code and memory statistics
namespace time_mem {
Timing tm;
MemoryStats ms;
void timing_report() { tm.report(); }
void memory_report() { ms.report(); }
void memory_time_brief_report() {
#ifdef HAS_MEMORY_USAGE
  cout << "Memory used: " << long(ms.used() / 1024) << " MB "; // NOLINT
#endif
  cout << "Time elapsed: " << prec3(tm.total_in_seconds()) << " s" << endl;
}
} // namespace time_mem

#ifdef NRG_MPI
mpi::environment *mpienv;
mpi::communicator *mpiw;
int myrank() { return mpiw->rank(); } // used in diag.h, time_mem.h
#else
int myrank() { return 0; }
#endif

#include <utility>

// Shared parameters for MPI parallelization.
class sharedparam {
  public:
  // Parameters which have to be known to the slave processes (which only
  // perform diagonalizations).
  dr_value diagroutine{undefined};
  double diagratio{};
  size_t dsyevrlimit{};
  size_t zheevrlimit{};
  bool logall{};
  string log;
  void init();
//  sharedparam() {};
  private:
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int version) {
    ar &diagroutine;
    ar &diagratio;
    ar &dsyevrlimit;
    ar &zheevrlimit;
    ar &logall;
    ar &log;
  }
};

sharedparam sP;

// Symmetry type specification
string sym_string = "";

// Quantum number types defined to enforce type checking
using Number = int;
using Ispin = int;
using Sspin = int;
using Tangmom = int;
using SZspin = int;

// Invariant subspace abstraction (container with quantum numbers)
#include "invar.cc"

// *** Commonly used types ***

#ifdef NRG_REAL
using t_matel = double;                       // type for the matrix elements
using t_eigen = double;                       // type for the eigenvalues
using t_coef = double;                        // type for the Wilson chain coefficients
using t_factor = double;                      // type for various prefactors in recalculations
using t_expv = double;                        // type for expectation values of operators
inline double CONJ_ME(double x) { return x; } // Conjugation of matrix elements: no op
#endif

#ifdef NRG_COMPLEX
using t_matel = cmpl;
using t_eigen = double;
using t_coef = cmpl;
using t_factor = cmpl;
using t_expv = cmpl; // we allow the calculation of expectation values of
                     // non-Hermitian operators!
inline cmpl CONJ_ME(cmpl z) { return conj(z); }
#endif

using t_weight = cmpl; // spectral weight accumulators (complex in general)

using EVEC = ublas::vector<t_eigen>; // Type for arrays of eigenvalues
using STDEVEC = std::vector<t_eigen>;
using CVEC = ublas::vector<t_coef>; // Type for arrays of coefficients for Wilson chains
using DVEC = ublas::vector<double>;
using IVEC = ublas::vector<size_t>;

/* NOTE: Row major is the C array format: A[0][0], A[0][1], A[0][2],
A[1][0], A[1][1], etc. The default in UBLAS is row major, while LAPACK
routines expect column major matrices. Of course, this is of no concern for
symmetric matrices. Default storage type is unbounded_array<T>.

Thus, as always:
  first index - row
  second index - column

 when accessing columns - stride=m.size2()
 when accessing rows - stide=1

 Optimization rule: use stride 1 sequential access where possible.
 ublas default matrix storage is row major (i.e. C-like). The rule is
 "right index the same as inner loop variable".
*/

using Matrix = matrix<t_matel, ublas::row_major>;

#include "numerics.h"

using InvarVec = std::vector<Invar>; // vector of Invars

using Twoinvar = pair<Invar, Invar>;
using MatrixElements = map<Twoinvar, Matrix>;
using DensMatElements = map<Invar, Matrix>;

using ThreeInvar = tuple<Invar, Invar, Invar>; // For 3-leg vertex functions

// Dump irreducible matrix elements between two subspaces
void dump_matrix_elements(const MatrixElements &m, ostream &fout = cout) {
  const double CHOPSMALL = 1e-8;
  const size_t MAXDUMP   = 10;
  for (const auto &i : m) {
    fout << "----" << i.first << "----" << endl;
    for (size_t r1 = 0; r1 < min(i.second.size1(), MAXDUMP); r1++) {
      for (size_t r2 = 0; r2 < min(i.second.size2(), MAXDUMP); r2++) fout << chop((i.second)(r1, r2), CHOPSMALL) << ' ';
      fout << endl;
    }
  }
}

// (Reduced) density matrix
DensMatElements rho;
DensMatElements rhoFDM;

template <typename T> inline pair<T, T> reverse_pair(const pair<T, T> &i) { return make_pair(i.second, i.first); }

// Map of operators matrices
using CustomOp = map<string, MatrixElements>;

// Name of the operator corresponding to a CustomOp map element
string NAME(const CustomOp::value_type &i) { return i.first; }

// Vector containing irreducible matrix elements of f operators.
using OpchChannel = std::vector<MatrixElements>;
// Each channel contains P::perchannel OpchChannel matrices.
using Opch = std::vector<OpchChannel>;

// Object of class IterInfo cotains full information about matrix representations 
// when entering stage N of the NRG iteration.
class IterInfo {
  public:
  Opch opch;     // f operators (channels)
  CustomOp ops;  // singlet operators (even parity)
  CustomOp opsp; // singlet operators (odd parity)
  CustomOp opsg; // singlet operators [global op]
  CustomOp opd;  // doublet operators (spectral functions)
  CustomOp opt;  // triplet operators (dynamical spin susceptibility)
  CustomOp opq;  // quadruplet operators (spectral functions for J=3/2)
  CustomOp opot; // orbital triplet operators

  // For debugging purposes we also store the information about the
  // ancestors of each invariant subspace. Ancestors are those invariant
  // subspaces which are combined with the states on the newly added
  // lattice site(s) to give each new invariant subspace.
  map<Invar, InvarVec> ancestors;

  // Erases all information in the structures of IterInfo object. This is
  // called at the start-up and before the DM-NRG run, so that the data
  // structures are properly reset. opch is erased in read_ireducf().
  void cleanup() {
    ops.clear();
    opsp.clear();
    opsg.clear();
    opd.clear();
    opt.clear();
    opq.clear();
    opot.clear();
    ancestors.clear();
  }
};

IterInfo iterinfo; // NOTE: global object! (directly used in matrix.cc, recalc.cc, nrg-recalc-*

// We need to store the dimensions of the invariant subspaces |r,1>,
// |r,2>, |r,3>, etc. The name "rmax" comes from the maximal value of
// the index "r" which ranges from 1 through rmax.

class Rmaxvals {
  private:
  IVEC values;
  friend ostream &operator<<(ostream &os, const Rmaxvals &rmax);
  void store(IVEC &rmx) {
    my_assert(rmx.size() == P::combs);
    values = rmx;
  }
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &values; }
  public:
  Rmaxvals() = default;
  Rmaxvals(const Rmaxvals &) = default;
  Rmaxvals(Rmaxvals &&) = default;
  Rmaxvals &operator=(const Rmaxvals &) = default;
  Rmaxvals &operator=(Rmaxvals &&) = default;
  ~Rmaxvals() = default;
  size_t rmax(size_t i) const {
    allowed_block_index(i);
    return values[i - 1]; // FOR COMPATIBILITY OFFSET 1!
  } 
  size_t offset(size_t i) const {
    allowed_block_index(i);
    return accumulate(begin(values), begin(values) + (i - 1), 0);
  }
  size_t operator[](size_t i) const { return rmax(i); }
  // The only way to set up the values in Rmaxvals is by calling
  // determine_ranges(), or by using the copy constructor.
  void determine_ranges(const Invar &I, const InvarVec &In);
  // total() returns the total number of states. There is therefore
  // no need to store this value separately.
  size_t total() const { return accumulate(begin(values), end(values), 0); }
};

// Information about the number of states, kept and discarded, rmax,
// and eigenenergies.
class DimSub {
  public:
  size_t kept      = 0;
  size_t discarded = 0;
  size_t total     = 0;
  Rmaxvals rmax;   // substructure of vectors omega
  EVEC eigenvalue; // all eigenvalues
  EVEC absenergy;  // absolute energies (for FDM)
  DimSub() = default;
  DimSub(size_t _kept, size_t _total, const Rmaxvals &_rmax, const EVEC &_eigenvalue, const EVEC &_absenergy) : // NOLINT
     kept(_kept), total(_total), rmax(_rmax), eigenvalue(_eigenvalue), absenergy(_absenergy) {
    my_assert(kept <= total);
    discarded = total - kept;
  }
  DimSub(const DimSub &) = default;
  DimSub(DimSub &&) = default;
  DimSub &operator=(const DimSub &) = default;
  DimSub &operator=(DimSub &&) = default;
  ~DimSub() = default;
};

// Full information about the number of states and matrix dimensions
// Example: dm[N].rmax[I] etc.
using Subs = map<Invar, DimSub>;
using AllSteps = std::vector<Subs>;
AllSteps dm;

// Result of a diagonalisation: eigenvalues and eigenvectors
class Eigen {
  public:
  size_t nr     = 0;   // number of eigenpairs (currently stored)
  size_t rmax   = 0;   // dimensionality of the matrix space
  size_t nrpost = 0;   // number of eigenpairs after truncation
  double shift  = 0.0; // shift of eigenvalues (0 or Egs)
  EVEC value;          // eigenvalues
  EVEC absenergy;      // absolute energies (0 is the absolute ground state of the system)
  EVEC boltzmann;      // Boltzmann factors
  Matrix matrix0;      // eigenvectors in matrix form
  // 'blocks' contains eigenvectors separated according to the invariant
  // subspace from which they originate. This separation is required for
  // using the efficient BLAS routines when performing recalculations of
  // the matrix elements.
  std::vector<Matrix> blocks;
  Eigen() = default;
  Eigen(const Eigen &t) = default;
  Eigen(Eigen &&t) = default;
  // nr - number of eigenpairs, rmax - dimensionality of the matrix space
  Eigen(size_t _nr, size_t _rmax) : nr(_nr), rmax(_rmax) {
    my_assert(rmax >= nr);
    value.resize(nr);
    absenergy.resize(nr);
    matrix0.resize(nr, rmax);
    perform_checks();
  }
  Eigen &operator=(const Eigen &) = default;
  Eigen &operator=(Eigen &&) = default;
  ~Eigen() = default;
  // Various assertion checks; to be called after the diagonalisation routine or reading Eigen objects through MPI or from disk.
  void perform_checks() const;
  // Accessor routine for j-th element of i-th eigenvector.
  inline t_matel &vektor(size_t i, size_t j) { return matrix0(i, j); }
  // Returns the number of eigenpairs CURRENTLY STORED.
  size_t getnr() const { return nr; }
  // Returns the dimensionality of a subspace, i.e. the number of
  // components of each eigenvector.
  size_t getrmax() const { return rmax; }
  size_t getdim() const { return rmax; }
  // Returns the number of eigenpairs after truncation.
  size_t getnrpost() const { return nrpost; }
  // Truncate to nrpost states.
  void truncate_prepare(size_t _nrpost) {
    nrpost = _nrpost;
    my_assert(nrpost <= nr);
  }
  void truncate_perform();
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors
  // form an identity matrix. This is used to represent the spectral
  // decomposition in the eigenbasis itself.
  void diagonal(EVEC &v) {
    nr = rmax = v.size();
    value     = v;
    shift     = 0.0;
    matrix0   = identity_matrix<t_eigen>(nr);
  }
  private:
  friend class boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, const unsigned int version) {
    ar &nr;
    ar &rmax;
    ar &value;
    ar &shift;
    ar &absenergy;
    ar &matrix0;
  }
};

void Eigen::truncate_perform() {
  for (auto &i : blocks) {
    my_assert(nrpost <= i.size1());
    i.resize(nrpost, i.size2());
  }
  value.resize(nrpost); // new, 19.7.2019
  nr = nrpost;
}

void Eigen::perform_checks() const {
  my_assert(value.size() == matrix0.size1());
  my_assert(matrix0.size1() <= matrix0.size2());
  my_assert(nr == value.size()); // new, 19.7.2019
  my_assert(nr == matrix0.size1());
  my_assert(rmax == matrix0.size2());
}

// Full information after diagonalizations.
using DiagInfo = map<Invar, Eigen>;

#define LOOP(diag, var) for (auto &var : diag) // NOLINT
#define LOOP_const(diag, var) for (const auto &var : diag) // NOLINT

Invar INVAR(const DiagInfo::value_type &i) { return i.first; }
// cppcheck-suppress constParameter symbolName=i
Eigen &EIGEN(DiagInfo::value_type &i) { return i.second; } // NOLINT
Eigen const &EIGEN(const DiagInfo::value_type &i) { return i.second; }

// Number of calculated states
size_t NRSTATES(const DiagInfo::value_type &i) { return i.second.getnr(); }

// Dimensionality of the subspace (dim)
size_t RMAX(const DiagInfo::value_type &i) { return i.second.getrmax(); }

ostream &operator<<(ostream &os, const Twoinvar &p) { return os << "(" << p.first << ") (" << p.second << ")"; }

template <typename T> ostream &operator<<(ostream &os, const ublas::vector<T> &v) {
  for (const auto &x : v) os << x << ' ';
  return os;
}

using QSrmax = map<Invar, Rmaxvals>;

class ChainSpectrum;
class BaseSpectrum;

using spCS_t = shared_ptr<ChainSpectrum>;

class SPEC {
  public:
  SPEC() = default;
  SPEC(const SPEC &) = default;
  SPEC(SPEC &&) = default;
  SPEC &operator=(const SPEC &) = default;
  SPEC &operator=(SPEC &&) = default;
  virtual ~SPEC() = default;
  virtual spCS_t make_cs(const BaseSpectrum &) = 0;
  virtual void calc(const Eigen &, const Eigen &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor, spCS_t, const Invar &,
                    const Invar &){};
  virtual void calc_A(const Eigen &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor,
                      spCS_t, const Invar &, const Invar &, const Invar &){};
  virtual void calc_B(const Eigen &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, const Matrix &, const BaseSpectrum &, t_factor,
                      spCS_t, const Invar &, const Invar &, const Invar &){};
  virtual string name() = 0;
  virtual string merge() { return ""; } // what merging rule to use
};

using SPECTYPE = shared_ptr<SPEC>;

// In namespace NRG we store run-time information about the calculation.
namespace NRG {
  // Flag to signal an insufficient number of states computed.
  bool notenough;
  // Diagratio for the following diagonalisation call
  double diagratio;
  // Invariant subspaces to be computed.
  std::vector<Invar> tasks;
#ifdef NRG_REAL
  const string v = "real";
#endif
#ifdef NRG_COMPLEX
  const string v = "complex";
#endif
} // namespace NRG

// Includes which require P:: parameters.
#include "spectral.h"

// The factor that multiplies the eigenvalues of the length-N Wilson
// chain Hamiltonian in order to obtain the energies on the original
// scale. Also named the "reduced bandwidth". Note that STAT::scale =
// SCALE(STAT::N+1) [see function set_N], thus for Ninit=0 and calc0,
// setting N=-1 results in a call SCALE(0). In the actual NRG run,
// the scale is at least SCALE(1). This is important for correct
// handling of rescaling for substeps==true.
double SCALE(int N) // int is correct here, this N may be negative
{
  double scale = 0.0;
  if (string(P::discretization) == "Y")
    // Yoshida,Whitaker,Oliveira PRB 41 9403 Eq. (39)
    scale = 0.5 * (1. + 1. / P::Lambda); // NOLINT
  if (string(P::discretization) == "C" || string(P::discretization) == "Z")
    // Campo, Oliveira PRB 72 104432, Eq. (46) [+ Lanczos]
    scale = (1.0 - 1. / P::Lambda) / log(P::Lambda); // NOLINT
  if (!P::substeps)
    scale *= pow(P::Lambda, -(N - 1) / 2. + 1 - P::z); // NOLINT
  else
    scale *= pow(P::Lambda, -N / (2. * P::channels) + 3 / 2. - P::z); // NOLINT
  my_assert(scale != 0.0);        // yes, != is intentional here.
  scale = scale * P::bandrescale; // RESCALE
  return scale;
}

// M ranges 0..channels-1
pair<size_t, size_t> get_Ntrue_M(size_t N) {
  size_t i = N / P::channels;
  return make_pair(i, N - i * P::channels);
}

// Compensate for different definition of SCALE in initial.m and C++
// code in case of substeps==true.
double scale_fix(size_t N) {
  size_t Ntrue, M;
  tie(Ntrue, M) = get_Ntrue_M(N);
  my_assert(N == Ntrue * P::channels + M);
  size_t N_at_end_of_full_step     = Ntrue * P::channels + P::channels - 1; // M=0,...,channels-1
  double scale_now                 = SCALE(N + 1); // NOLINT
  double scale_at_end_of_full_step = SCALE(N_at_end_of_full_step + 1); // NOLINT
  return scale_now / scale_at_end_of_full_step;
}

enum class RUNTYPE { NRG, DMNRG };

// Map from double to double. This is used for the "fixeps" trick as a data
// structure which holds the transformation rules from original eigenvalues
// to fixed eigenvalues.
using mapdd = unordered_map<t_eigen, t_eigen>;

// Pair of energy and multiplicity. This is used to compute the true
// grand-canonical partition function STAT::ZZ.
using energy_mult_type = pair<t_eigen, int>;
using excitation_list = list<energy_mult_type>;

// Namespace for storing various statistical quantities calculated
// during iteration.
namespace STAT {
  RUNTYPE runtype; // NRG vs. DM-NRG run
  size_t N;        // iteration step, N>=0
  double scale;    // scale factor
  double Teff;     // effective temperature
  double scT;      // scT = scale/P::T, this combination appears in the
                   // exponents (Boltzmann weights)
  /* After the diagonalization, the energy scale is reduced by a factor of
   \Lambda^{1/2}. True energies are now eigenvalues multiplied by this new
   scale. The effective temperature is T_(N+1)=scale/P::betabar. Note that
   E/T = E_N * scale / (scale / betabar) = betabar E_N. P::betabar
   determines the factor c in the expression T_N = c Lambda^-(N-1)/2. */
  void set_N(int newN) {
    N = (newN >= 0 ? newN : 0);                                          // N is used as an array index, must be non-negative
    // The following is evaluated with newN, which may be negative!
    scale = SCALE(newN + 1); // Current energy scale in units of bandwidth D
    Teff  = scale / P::betabar;
    scT   = scale / P::T; // Boltzmann weight
  }
  // Ground-state subspace and energy at the current NRG iteration
  Invar ground;
  t_eigen Egs;
  // All energies, sorted in increasing order
  STDEVEC energies;
  // Mapping old_energy -> new_energy for curing small splittings due to
  // the floating point roundoff errors.
  mapdd cluster_mapping;
  double Zft;   // grand-canonical partition function (at shell n)
  double Zgt;   // grand-canonical partition function for computing G(T)
  double Zchit; // grand-canonical partition function for computing chi(T)
  double E;     // energy times beta
  double E2;    // (energy times beta) squared
  double C;     // heat capcity (in units of 1/k_B)
  double F;     // free energy times beta
  double S;     // entropy: (beta E)-(beta F)
  // Expectation values
  using t_mapexpv = map<string, t_expv>;
  t_mapexpv expv;    // expectation values of custom operators
  t_mapexpv fdmexpv; // Expectation values computed using the FDM algorithm
                     // Total energy of the ground state. This is the sum of all the zero state
  // energies for all the iteration steps.
  t_eigen totalenergy;
  // This is the value of the variable "totalenergy" at the end of the
  // iteration. This is different from 'Egs'.
  t_eigen GSenergy;
  // Containers related to the FDM-NRG approach.
  // Consult A. Weichselbaum, J. von Delft, PRL 99, 076402 (2007).
  std::vector<double> ZnD;      // Z_n^D=\sum_s^D exp(-beta E^n_s),
                                // sum over **discarded** states at shell n
                                // cf. Eq. (8) and below in WvD paper
  std::vector<double> wn;       // Weights w_n. They sum to 1.
  std::vector<double> wnfactor; // wn/ZnD
  double ZZ = 0;                // grand-canonical partition function (full-shell),
                                // evaluated at the temperature P::T
} // namespace STAT

// Return true if this is the first step of the NRG iteration
bool FIRST_ITERATION(size_t N = STAT::N) { return N == P::Ninit; }

// Return true if this is the last step of the NRG iteration
bool LAST_ITERATION(size_t N = STAT::N) {
  return (N + 1) == P::Nmax || (P::ZBW && N == P::Ninit); // special case!
}

// NOTE: for zero-bandwidth calculations, Ninit=0 and Nmax=0, so that
// FIRST_ITERATION(0) == true and LAST_ITERATION(0) == true.
// The same is true for Nmax=1.

// Test for first NRG run [diagonalisations] with (nrgrun)
class runtypenrg {
  public:
  operator bool() { return STAT::runtype == RUNTYPE::NRG; }
} nrgrun;

// Test for second NRG run [density matrix computations] with (dmnrgrun)
class runtypedmnrg {
  public:
  operator bool() { return STAT::runtype == RUNTYPE::DMNRG; }
} dmnrgrun;

ostream &operator<<(ostream &os, const Rmaxvals &rmax) {
  for (const auto &x : rmax.values) os << x << ' ';
  return os;
}

// Data structures containing NRG eigenvalues and eigenvectors

DiagInfo diag;     // holds eigenvalues and unitary matrices for all subspaces
DiagInfo diagprev; // holds eigenenegies (from previous/current iteration and initial ones)
QSrmax qsrmax;     // holds the number of states for |r,1>, |r,2>,...

// Returns true if option 'c' is selected for logging
bool logletter(char c) { return (sP.logall ? true : sP.log.find(c) != string::npos); }

// Index 'n' of the last site in the existing chain, f_n (at iteration 'N').
// The site being added is f_{n+1}.
// This is the value that we use in building the matrix,
// cf. nrg-makematrix-ISO.cc
int getnn() { return STAT::N; }

// Energy scale at the last NRG iteration
double LAST_STEP_SCALE() { return SCALE(P::Nmax); }

void validate_parameters();

map<string, string> parsed_params;

// Parse command line and parameter file
// Warning: there's no error checking for parameters.
// Check the parameter dump to see if everything is OK.
void read_parameters() {
  parsed_params.clear(); // XXX: remove after building a nrg class
  parser(parsed_params, "param");
  for (const auto &i : allparams) {
    const string keyword = i->getkeyword();
    if (parsed_params.count(keyword) == 1) {
      i->setvalue_str(parsed_params[keyword]);
      parsed_params.erase(keyword);
    }
  }
  if (parsed_params.size()) {
    cout << "Unused settings: " << endl;
    for (const auto &i : parsed_params) cout << " " << i.first << "=" << i.second << endl;
    cout << endl;
  }
}

void dump_parameters() {
  allparams.sort([](auto a, auto b) { return a->getkeyword() < b->getkeyword(); });
  for (const auto &i : allparams) i->dump();
}

void remove_workdir() { remove(P::workdir.c_str()); }

const string default_workdir = ".";

void create_workdir(const string &workdir) {
  const string workdir_template = workdir + "/XXXXXX";
  size_t len = workdir_template.length()+1;
  auto x = make_unique<char[]>(len); // NOLINT
  strncpy(x.get(), workdir_template.c_str(), len);
  if (char *w = mkdtemp(x.get())) // create a unique directory
    P::workdir = w;
  else
    P::workdir = default_workdir;
  cout << "workdir=" << P::workdir << endl << endl;
  atexit(remove_workdir);
}

void set_workdir(const string &workdir_) {
  string workdir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) workdir = env_w;
  if (!workdir_.empty()) workdir = workdir_;
  create_workdir(workdir);
}

void set_workdir(int argc, char **argv) {
  std::vector<string> args(argv+1, argv+argc); // NOLINT
  string workdir = default_workdir;
  if (const char *env_w = std::getenv("NRG_WORKDIR")) workdir = env_w;
  if (args.size() == 2 && args[0] == "-w") workdir = args[1];
  create_workdir(workdir);
}

// This class holds table of generalized xi/zeta/etc. coefficients
class coef_table {
  private:
  CVEC t;

  public:
  // Read values from a stream f
  void read_values(ifstream &f) {
    size_t len;
    f >> len; // get length (= last index n still included)
    read_vector(f, t, len + 1);
  }
  t_coef coef(size_t n) const {
    my_assert(n < t.size());
    return t[n];
  }
  // Returns the index of the last coefficient still included in the table.
  size_t max() const {
    my_assert(t.size() >= 1);
    return t.size() - 1;
  }
  void setvalue(size_t n, t_coef val) {
    if (n + 1 > t.size()) t.resize(n + 1);
    t[n] = val;
  }
};

// NOTE: One table of discretization coefficients for each channel
class set_of_tables {
  private:
  std::vector<coef_table> tabs;

  public:
  size_t nr_tabs() const { return tabs.size(); }
  void read(ifstream &fdata) {
    tabs.resize(P::coefchannels);
    for (auto &i : tabs) i.read_values(fdata);
  }
  t_coef operator()(size_t N, size_t alpha) const {
    allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size());
    return tabs[alpha].coef(N);
  }
  size_t max(size_t alpha) const {
    allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size());
    return tabs[alpha].max();
  }
  void setvalue(size_t N, size_t alpha, t_coef val) {
    allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size() && N <= P::Nmax);
    tabs[alpha].setvalue(N, val);
  }
};

set_of_tables xi;   // f^dag_N f_N+1 terms
set_of_tables zeta; // f^dag_N f_N terms

// Support for spin-polarized conduction bands. See also P::polarized.
// Hack: the total number of channels is doubled, the index runs from
// 0...2*P::channels-1. Numbers 0...P::channels-1 correspond to spin up,
// while P::channels...2*P::channels-1 correspond to spin down.
// Compare P::channels and P::coefchannels (which reflects the same
// convention in initial.m, i.e. CHANNELS vs. COEFCHANNELS).
t_coef xiUP(size_t N, size_t ch) { return xi(N, ch); }
t_coef xiDOWN(size_t N, size_t ch) { return xi(N, ch + P::channels); }
t_coef zetaUP(size_t N, size_t ch) { return zeta(N, ch); }
t_coef zetaDOWN(size_t N, size_t ch) { return zeta(N, ch + P::channels); }

// Support for conduction bands with full 2x2 matrix structure, a
// generalization of P::polarized. The total number of "channels" is
// here multiplied by 4, i.e., the index runs from 0 to
// 4*P::channels-1. Numbers 2*P::channels...3*P::channels-1 correspond
// to UP/DO, 3*P::channels...4*P::channels-1 correspond to DO/UP.
t_coef xiUPDO(size_t N, size_t ch) { return xi(N, ch + 2 * P::channels); }
t_coef xiDOUP(size_t N, size_t ch) { return xi(N, ch + 3 * P::channels); }
t_coef zetaUPDO(size_t N, size_t ch) { return zeta(N, ch + 2 * P::channels); }
t_coef zetaDOUP(size_t N, size_t ch) { return zeta(N, ch + 3 * P::channels); }

// Support for channel-mixing Wilson chains
set_of_tables xiR;
set_of_tables zetaR;
// Support for superconducting bands.
set_of_tables delta; // f^dag_up,N f^dag_down,N terms
set_of_tables kappa; // f^dag_N f^dag_down,N+1 terms

set_of_tables ep, em;   // e_n coefficients
set_of_tables u0p, u0m; // u_{0,m} coefficients

#include "tridiag.h"
#include "diag.h"
#include "symmetry.cc"
#include "matrix.cc"
#include "recalc.cc"

// Select which symmetries to compile in.
#ifdef NRG_SYM_BASIC
#include "sym-QS.cc"
#include "sym-QSZ.cc"
#endif
#ifdef NRG_SYM_MORE
#include "sym-ISO.cc"
#include "sym-ISOSZ.cc"
#include "sym-SPSU2.cc"
#include "sym-SPU1.cc"
#endif
#ifdef NRG_SYM_ALL
#include "sym-DBLSU2.cc"
#include "sym-DBLISOSZ.cc"
#include "sym-ISOLR.cc"
#include "sym-ISOSZLR.cc"
#include "sym-NONE.cc"
#include "sym-P.cc"
#include "sym-PP.cc"
#include "sym-SL.cc"
#include "sym-SL3.cc"
#include "sym-SPSU2LR.cc"
#include "sym-SPSU2T.cc"
#include "sym-SPSU2C3.cc"
#include "sym-SPU1LR.cc"
#include "sym-SU2.cc"
#include "sym-QSLR.cc"
#include "sym-QSC3.cc"
#include "sym-QST.cc"
#include "sym-QSTZ.cc"
#include "sym-QSZTZ.cc"
#include "sym-QSZLR.cc"
#include "sym-QJ.cc"
#include "sym-U1.cc"
#endif

#include "read-input.cc"

/**** Calculation of traces ****/

// Used in calculate_TD().
CONSTFNC double calculate_Z(const DiagInfo::value_type &is, double rescale_factor) {
  double sumZ = 0;
  for (const auto &x : EIGEN(is).value) sumZ += exp(-rescale_factor * x);
  sumZ *= mult(INVAR(is));
  return assert_isfinite(sumZ);
}

// Formated output for the expectation values
template <typename T> string output_val(const T &x, size_t prec = std::numeric_limits<double>::max_digits10) {
  ostringstream F;
  F << setprecision(prec) << x;
  return F.str();
}

// Specialization for complex values, the output format is X+IY or
// X-IY, where X and Y are real and imaginary part, respectively. The
// imaginary part is only shown where its value relative to the real
// part is sufficiently large. No space is used in the outputted
// string in order to simplify parsing. This behavior can be turned
// off using P::noimag, which is the default.
const double OUTPUT_IMAG_EPS = 1.0e-13;
string output_val(const cmpl &val) {
  ostringstream F;
  if (P::noimag || abs(val.imag()) < abs(val.real()) * OUTPUT_IMAG_EPS) {
    F << val.real();
  } else {
    F << val.real();
    if (val.imag() > 0.0)
      F << "+I" << val.imag();
    else
      F << "-I" << -val.imag();
  }
  return F.str();
}

template <typename T> void formatted_output(ostream &F, T x) {
  // Important: setw first, setprecision second
  F << setw(P::width_custom) << setprecision(P::prec_custom) << x << ' ';
}

void formatted_output(ostream &F, cmpl val) {
  ostringstream str;
  // This sets precision for both real and imaginary parts.
  str << setprecision(P::prec_custom);
  if (P::noimag || abs(val.imag()) < abs(val.real()) * OUTPUT_IMAG_EPS) {
    str << val.real();
  } else {
    str << val.real();
    if (val.imag() > 0.0)
      str << "+I" << val.imag();
    else
      str << "-I" << -val.imag();
  }
  // The width for the whole X+IY string.
  F << setw(P::width_custom) << str.str() << ' ';
}

/* "Trace" of a singlet operator: actually this is the statistical average
 with respect to exp(-beta H), but without the 1/Z factor. I.e.
 \sum_n <n|exp(-beta H) O|n>, for a given singlet operator 'O'.  As a
 side effect, dump the matrix elements to stream F if so requested
 (parameter "dumpdiagonal").  */
CONSTFNC t_expv calc_trace_singlet(const DiagInfo &diag, const MatrixElements &n, ostream &F = cout) {
  matel_bucket tr; // note: t_matel = t_expv
  LOOP_const(diag, is) {
    const Twoinvar I = make_pair(INVAR(is), INVAR(is));
    if (n.count(I) != 1) exit1("ERROR: calc_trace_singlet() I=(" << I << ") cnt=" << n.count(I));
    const Matrix &nI = n.find(I)->second;
    const size_t dim = NRSTATES(is);
    my_assert(dim == nI.size2());
    matel_bucket sum;
    for (size_t r = 0; r < dim; r++) sum += exp(-P::betabar * EIGEN(is).value(r)) * nI(r, r);
    if (P::dumpdiagonal != 0) {
      F << INVAR(is) << ": ";
      for (size_t r = 0; r < dim; r++)
        if (r < P::dumpdiagonal) F << nI(r, r) << ' ';
      F << endl;
    }
    tr += t_matel(mult(INVAR(is))) * t_matel(sum);
  }
  return tr;
}

CONSTFNC t_expv calc_trace_fdm_kept(const DiagInfo &diag, const MatrixElements &n) {
  matel_bucket tr;
  LOOP_const(diag, is) {
    const Invar I     = INVAR(is);
    const Twoinvar II = make_pair(INVAR(is), INVAR(is));
    my_assert(n.count(II) == 1);
    const Matrix &nI  = n.find(II)->second;
    const size_t ret  = dm[STAT::N][I].kept;
    matel_bucket sum;
    for (size_t i = 0; i < ret; i++) // over kept states ONLY
      for (size_t j = 0; j < ret; j++) sum += rhoFDM[I](i, j) * nI(j, i);
    tr += t_matel(mult(INVAR(is))) * t_matel(sum);
  }
  return tr;
}

#include "bins.h"

/* Object of class 'ChainSpectrum' will contain information about the
 spectral density calculated at a given stage of the NRG run, i.e. for a
 finite Wilson chain.  We then merge it in an object of class 'Spectrum'
 which holds the spectral information for the entire run (i.e. the physical
 spectral density). */

class ChainSpectrum {
  public:
  virtual void add(double energy, t_weight weight) = 0;
  ChainSpectrum() = default;
  ChainSpectrum(const ChainSpectrum &) = default;
  ChainSpectrum(ChainSpectrum &&) = default;
  ChainSpectrum &operator=(const ChainSpectrum &) = default;
  ChainSpectrum &operator=(ChainSpectrum &&) = default;
  virtual ~ChainSpectrum() = default; // required, because there are virtual members
};

class ChainSpectrumBinning : public ChainSpectrum {
  private:
  Bins spos, sneg;
  public:
  void add(double energy, t_weight weight) override {
    if (energy >= 0.0)
      spos.add(energy, weight);
    else
      sneg.add(-energy, weight);
  }
  t_weight total_weight() const { return spos.total_weight() + sneg.total_weight(); }
  friend class SpectrumRealFreq;
};

class ChainSpectrumTemp : public ChainSpectrum {
  private:
  Temp v;
  public:
  void add(double T, t_weight value) override { v.add_value(T, value); }
  friend class SpectrumTemp;
};

#include "matsubara.h"
#include "matsubara2.h"

class ChainSpectrumMatsubara : public ChainSpectrum {
  private:
  Matsubara m;
  public:
  ChainSpectrumMatsubara() = delete;
  explicit ChainSpectrumMatsubara(matstype _mt) : m(P::mats, _mt){};
  void add(size_t n, t_weight w) { m.add(n, w); }
  void add(double energy, t_weight w) override { my_assert_not_reached(); }
  t_weight total_weight() const { return m.total_weight(); }
  friend class SpectrumMatsubara;
};

class ChainSpectrumMatsubara2 : public ChainSpectrum {
  private:
  Matsubara2 m;
  public:
  ChainSpectrumMatsubara2() = delete;
  explicit ChainSpectrumMatsubara2(matstype _mt) : m(P::mats, _mt){};
  void add(size_t i, size_t j, t_weight w) { m.add(i, j, w); }
  void add(double energy, t_weight w) override { my_assert_not_reached(); }
  t_weight total_weight() const { return m.total_weight(); }
  friend class SpectrumMatsubara2;
};

// Object of class spectrum will contain everything that we know about a
// spectral density.
class Spectrum {
  public:
  string opname, filename;
  SPECTYPE spectype;
  Spectrum(const string &_opname, const string &_filename, SPECTYPE _spectype) : opname(_opname), filename(_filename), spectype(_spectype){}; // NOLINT
  Spectrum(const Spectrum &) = default;
  Spectrum(Spectrum &&) = default;
  Spectrum &operator=(const Spectrum &) = default;
  Spectrum &operator=(Spectrum &&) = default;
  virtual ~Spectrum()= default; // required (the destructor saves the results to a file)
  virtual void merge(spCS_t) = 0; // called from spec.cc as the very last step
  string name() { return opname; }
};

#include "spectrumrealfreq.cc"

// G(T) type of results, i.e. not a real spectrum
class SpectrumTemp : public Spectrum {
  private:
  std::vector<pair<double, t_weight>> results;
  public:
  SpectrumTemp(const string &_opname, const string &_filename, SPECTYPE _spectype) : Spectrum(_opname, _filename, _spectype) {}
  void merge(spCS_t) override;
  SpectrumTemp(const SpectrumTemp &) = default;
  SpectrumTemp(SpectrumTemp &&) = default;
  SpectrumTemp &operator=(const SpectrumTemp &) = default;
  SpectrumTemp &operator=(SpectrumTemp &&) = default;
  ~SpectrumTemp() override;
};

void SpectrumTemp::merge(spCS_t cs) {
  auto t = dynamic_pointer_cast<ChainSpectrumTemp>(cs);
  copy(begin(t->v), end(t->v), back_inserter(results));
}

SpectrumTemp::~SpectrumTemp() {
  string fn = filename + ".dat";
  cout << "Spectrum: " << opname << " " << spectype->name() << " -> " << fn << endl;
  Spikes d(results);
  sort(begin(d), end(d), sortfirst());
  ofstream Fd = safeopen(fn);
  save_densfunc(Fd, d, P::reim);
}

// This container actually holds the GF on the Matsubara axis, not a
// spectral function.
class SpectrumMatsubara : public Spectrum {
  private:
  Matsubara results;
  public:
  SpectrumMatsubara(const string &_opname, const string &_filename, SPECTYPE _spectype, matstype _mt)
     : Spectrum(_opname, _filename, _spectype), results(P::mats, _mt) {}
  void merge(spCS_t) override;
  SpectrumMatsubara(const SpectrumMatsubara &) = default;
  SpectrumMatsubara(SpectrumMatsubara &&) = default;
  SpectrumMatsubara &operator=(const SpectrumMatsubara &) = default;
  SpectrumMatsubara &operator=(SpectrumMatsubara &&) = default;
  ~SpectrumMatsubara() override;
};

void SpectrumMatsubara::merge(spCS_t cs) {
  auto t = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
  nrglog('*', "weight=" << t->total_weight()); // useful for debugging
  for (size_t n = 0; n < P::mats; n++) results.v[n].second += t->m.v[n].second;
}

SpectrumMatsubara::~SpectrumMatsubara() {
  cout << "Spectrum: " << opname << " " << spectype->name() << endl;
  ofstream Fd = safeopen(filename + ".dat");
  results.save(Fd);
}

class SpectrumMatsubara2 : public Spectrum {
  private:
  Matsubara2 results;
  public:
  SpectrumMatsubara2(const string &_opname, const string &_filename, SPECTYPE _spectype, matstype _mt)
     : Spectrum(_opname, _filename, _spectype), results(P::mats, _mt) {}
  void merge(spCS_t) override;
  SpectrumMatsubara2(const SpectrumMatsubara2 &) = default;
  SpectrumMatsubara2(SpectrumMatsubara2 &&) = default;
  SpectrumMatsubara2 &operator=(const SpectrumMatsubara2 &) = default;
  SpectrumMatsubara2 &operator=(SpectrumMatsubara2 &&) = default;
  ~SpectrumMatsubara2() override;
};

void SpectrumMatsubara2::merge(spCS_t cs) {
  auto t = dynamic_pointer_cast<ChainSpectrumMatsubara2>(cs);
  nrglog('*', "weight=" << t->total_weight());
  for (size_t m = 0; m < P::mats; m++)
    for (size_t n = 0; n < P::mats; n++) results.v(m, n) += t->m.v(m, n);
}

SpectrumMatsubara2::~SpectrumMatsubara2() {
  cout << "Spectrum: " << opname << " " << spectype->name() << endl;
  ofstream Fd = safeopen(filename + ".dat");
  results.save(Fd);
}

// This is mathematical trace, i.e. the sum of the diagonal elements.
CONSTFNC double trace(const DensMatElements &m) {
  double tr = 0.0;
  for (const auto &i : m) tr += mult(i.first) * trace_real_nochecks(i.second);
  return tr;
}

// Check if the trace of the density matrix equals 'ref_value'.
void check_trace_rho(const DensMatElements &m, double ref_value = 1.0) {
  const double tr = trace(m);
  if (!num_equal(trace(m), ref_value)) exit1("check_trace_rho() failed: " << tr << " instead of " << ref_value);
}

enum class axis { RealFreq, Temp, Matsubara, Matsubara2 };

string axisstring(axis a) {
  switch (a) {
    case axis::RealFreq: return "RealFreq";
    case axis::Temp: return "Temp";
    case axis::Matsubara: return "Matsubara";
    case axis::Matsubara2: return "Matsubara,Matsubara";
    default: my_assert_not_reached();
  }
}

ostream &operator<<(ostream &os, const axis a) { return os << axisstring(a); }

/* class BaseSpectrum contains all information about calculating the
 spectrum: pointers to the operator data and miscelaneous data, such
 as the spectrum type. Functions calc_specdens() et al. receive an
 object of this type as input. */
class BaseSpectrum {
  public:
  string name;
  string prefix; // "dens", "corr", etc.
  size_t nr;     // number of operators
  const MatrixElements &op1, &op2, &op3;
  shared_ptr<Spectrum> spec;
  SPECTYPE spectype{}; // SPEC_FT, ...
  axis a;            // axis::RealFreq, axis::Temp, axis::Matsubara, etc.
  matstype mt;       // matstype::bosonic, matstype::fermionic, etc.
  int spin{};          // -1 or +1, or 0 where irrelevant
  string fullname() const {
    string s = name + " " + prefix + " " + spectype->name() + " " + axisstring(a);
    if (a != axis::RealFreq && a != axis::Temp) s += " " + matstypestring(mt);
    return s;
  }
  void about() { cout << "Spectrum: " << fullname() << endl; }
  BaseSpectrum(const MatrixElements &_op1, const MatrixElements &_op2) : 
     op1(_op1), op2(_op2), op3(_op2), a(axis::RealFreq), mt(matstype::fermionic) { nr = 2; } // op3 initialization is a hack
  BaseSpectrum(const MatrixElements &_op1, const MatrixElements &_op2, const MatrixElements &_op3) : 
     op1(_op1), op2(_op2), op3(_op3), a(axis::RealFreq), mt(matstype::fermionic) { nr = 3; }
};

class speclist;
using lsl = list<speclist *>;
lsl allspectra; // list of list of spectra

class speclist {
  private:
  list<BaseSpectrum> spectra;

  public:
  speclist() { allspectra.push_back(this); }
  auto begin() { return spectra.begin(); }
  auto end() { return spectra.end(); }
  void push_back(BaseSpectrum &bs) { spectra.push_back(bs); }
  // Broaden spectra, close spectral files and deallocate all data storage (in destructor!)
  void clear() { spectra.clear(); }
  void about() {
    for (auto &i : spectra) i.about();
  }
};

speclist spectraD, spectraS, spectraT, spectraQ, spectraGT, spectraI1T, spectraI2T, spectraK, spectraCHIT, spectraC, spectraOT, spectraV3;

/**** CALCULATION OF SPECTRAL FUNCTIONS ****/

auto CorrelatorFactorFnc = [](const Invar &Ip, const Invar &I1) { return mult(I1); };
auto SpecdensFactorFnc = [](const Invar &Ip, const Invar &I1) { return Sym->specdens_factor(Ip, I1); };
auto SpecdensquadFactorFnc = [](const Invar &Ip, const Invar &I1) { return Sym->specdensquad_factor(Ip, I1); };
auto SpinSuscFactorFnc = [](const Invar &Ip, const Invar &I1) { return Sym->dynamicsusceptibility_factor(Ip, I1); };
auto OrbSuscFactorFnc = [](const Invar &Ip, const Invar &I1) { return Sym->dynamic_orb_susceptibility_factor(Ip, I1); };
auto TrivialCheckSpinFnc = [](const Invar &, const Invar &, int) { return true; };
auto SpecdensCheckSpinFnc = [](const Invar &I1, const Invar &Ip, int SPIN) { return Sym->check_SPIN(I1, Ip, SPIN); };

void doublet_check_norm(const CustomOp::value_type &op, const DiagInfo &diag, int SPIN) {
  weight_bucket sum;
  LOOP_const(diag, isp) {
    const Invar Ip = INVAR(isp);
    LOOP_const(diag, is1) {
      const Invar I1    = INVAR(is1);
      const Twoinvar II = make_pair(I1, Ip);
      auto it = op.second.find(II);
      if (it != op.second.end()) {
        if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
        t_factor spinfactor = Sym->specdens_factor(Ip, I1);
        auto mat = it->second;
        for (size_t r1 = 0; r1 < mat.size1(); r1++)
          for (size_t rp = 0; rp < mat.size2(); rp++) sum += spinfactor * sqr(abs(mat(r1, rp)));
      }
    }
  }
  const double result = 2.0 * cmpl(sum).real();
  // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
  cout << "check_norm[" << NAME(op) << "]=" << result << endl;
}

void quadruplet_check_norm(const CustomOp::value_type &op, const DiagInfo &diag, int SPIN) {
  weight_bucket sum;
  LOOP_const(diag, isp) {
    const Invar Ip = INVAR(isp);
    LOOP_const(diag, is1) {
      const Invar I1    = INVAR(is1);
      const Twoinvar II = make_pair(I1, Ip);
      auto it = op.second.find(II);
      if (it != op.second.end()) {
        if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
        t_factor spinfactor = Sym->specdensquad_factor(Ip, I1);
        auto mat = it->second;
        for (size_t r1 = 0; r1 < mat.size1(); r1++)
          for (size_t rp = 0; rp < mat.size2(); rp++) sum += spinfactor * sqr(abs(mat(r1, rp)));
      }
    }
  }
  const double result = 2.0 * cmpl(sum).real();
  // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
  cout << "check_norm[" << NAME(op) << "]=" << result << endl;
}

#include "spec.cc"

InvarVec input_subspaces() { return In; }

#include "dmnrg.h"

// **** Helper functions for the NRG RUN ****

void copy_sort_energies(const DiagInfo &diag, STDEVEC &energies) {
  energies.clear();
  LOOP_const(diag, is) energies.insert(end(energies), begin(EIGEN(is).value), end(EIGEN(is).value));
  sort(begin(energies), end(energies));
}

#include "splitting.cc"

// Number of eigenstates in subspace I at the previous (!) iteration.
size_t size_subspace_prev(const Invar &I) { return diagprev[I].getnr(); }

// Generate an info string for headers
void infostring() {
  string info = " ***** [" + (string)(nrgrun ? "NRG" : "DM") + "] " + "Iteration " + to_string(STAT::N + 1) + "/" + to_string(P::Nmax) + " (scale "
     + to_string(STAT::scale) + ")" + " ***** ";
  if (P::substeps) {
    size_t Ntrue, M;
    tie(Ntrue, M) = get_Ntrue_M(STAT::N);
    info += " step " + to_string(Ntrue + 1) + " substep " + to_string(M + 1);
  }
  cout << endl << info << endl;
}

// TODO: generalize to all symmetry types! (i.e. add additional
// coefficients of interest). This should ideally be a call to a
// method in class Symmetry.
void show_coefficients() {
  cout << setprecision(std::numeric_limits<double>::max_digits10);
  if (!P::substeps) {
    using namespace STAT;
    for (size_t i = 0; i < P::coefchannels; i++) {
      double scale = SCALE(static_cast<int>(N)+1);
      cout << "[" << i + 1 << "]"
           << " xi(" << N << ")=" << xi(N, i) << " xi_scaled(" << N << ")=" << xi(N, i)/scale << " zeta(" << N + 1 << ")=" << zeta(N + 1, i)
           << endl;
    }
  } else {
    size_t Ntrue, M;
    tie(Ntrue, M) = get_Ntrue_M(STAT::N);
    for (size_t i = 0; i < P::coeffactor; i++) {
      size_t index = M + P::channels * i;
      cout << "[" << index << "]"
           << " xi(" << Ntrue << ")=" << xi(Ntrue, index) << " zeta(" << Ntrue + 1 << ")=" << zeta(Ntrue + 1, index) << endl;
    }
  }
  Sym->show_coefficients();
}

// Information about ancestor subspaces
void logancestors(const Invar &I, const InvarVec &In, const Rmaxvals &rmx) {
  if (logletter('s')) {
    cout << "Ancestors of (" << I << "): ";
    for (size_t i = 1; i <= P::combs; i++)
      if (rmx[i] > 0) cout << "(" << In[i] << ", dim=" << rmx[i] << ") ";
    cout << endl;
  }
}

// Find the ground state in the current NRG shell.
void find_groundstate(const DiagInfo &diag) {
  STAT::Egs = DBL_MAX;
  LOOP_const(diag, is) {
    my_assert(EIGEN(is).value.size() > 0);
    t_eigen Emin = EIGEN(is).value(0); // Eigenvalues are sorted
    Emin += EIGEN(is).shift;           // if something was subtracted, we add it back
    if (Emin < STAT::Egs) {
      STAT::Egs    = Emin;
      STAT::ground = INVAR(is);
    }
  }
}

// Dump all energies in diag to a file
void dumptofile(const DiagInfo &diag, ostream &file) {
  file << endl << "===== Iteration number: " << STAT::N << endl;
  LOOP_const(diag, is) {
    file << "Subspace: " << is.first << endl;
    for (const auto &x : EIGEN(is).value) file << x << ' ';
    file << endl;
  }
}

/*** Subtract ground state energy ***/
void subtract_groundstate_energy(DiagInfo &diag) {
  if (nrgrun) // should be done only once!
    LOOP(diag, is) {
      for (auto &x : EIGEN(is).value) x -= STAT::Egs;
      EIGEN(is).shift = STAT::Egs;
      my_assert(EIGEN(is).value[0] >= 0.0);
    }
}

// Calculate the (on-shell) statistical sum using the excitation
// energies relative to the current step energy scale.
// Used in nrg_measure_singlet().
double calc_Z(const DiagInfo &diag) {
  bucket Z;
  LOOP_const(diag, is) for (const auto &x : EIGEN(is).value) Z += mult(INVAR(is)) * exp(-P::betabar * x);
  return Z;
}

// Function newcombination_allowed() checks if states |q,ss,i>_{N+1} can be
// formed for a given spin. Some combinations might be forbidden even when
// the corresponding |q',ss'> subspaces exist at iteration step N (Invar In).
bool newcombination_allowed(size_t i, const Invar &I, const Invar &In) { return Sym->triangle_inequality(I, In, QN[i]); }

// Determine the ranges of index r
void Rmaxvals::determine_ranges(const Invar &I, const InvarVec &In) {
  IVEC rmx(P::combs);
  rmx.clear();
  for (size_t i = 0; i < P::combs; i++)
    if (newcombination_allowed(i + 1, I, In[i + 1])) rmx[i] = size_subspace_prev(In[i + 1]);
  store(rmx);
}

// *********************************** NRG RUN **********************************

ofstream Fenergies;  // all energies (different file for NRG and for DMNRG)
ofstream Ftd;        // magnetic and charge susceptibility
ofstream Fannotated; // annotated eigenvalue spectrum

// Construct the suffix of the filename for spectral density files: 'A_?-A_?'.
// If SPIN == 1 or SPIN == -1, '-u' or '-d' is appended to the string.
string SDNAME(const string &a, const string &b, int SPIN = 0) {
  string sdname = a + "-" + b;
  if (SPIN == 1) sdname += "-u";
  if (SPIN == -1) sdname += "-d";
  return sdname;
}

// Formatted output of the computed expectation values
class ExpvOutput {
  private:
  ofstream F;             // output stream
  map<string, t_expv> &m; // reference to the name->value mapping
  list<string> fields;    // list of fields to be output (may be a subset
                          // of the fields actually present in m)
  // Consecutive numbers for the columns
  void field_numbers() {
    F << "#";
    formatted_output(F, 1);
    for (size_t ctr = 1; ctr <= fields.size(); ctr++) formatted_output(F, 1 + ctr);
    F << endl;
  }
  // Label and field names. Label is the first column (typically the
  // temperature).
  void field_names(string labelname = "T") {
    F << "#";
    formatted_output(F, labelname);
    for (const auto &op : fields) formatted_output(F, op);
    F << endl;
  }
  public:
  // Output the current values for the label and for all the fields
  void field_values(double labelvalue = STAT::Teff) {
    F << ' ';
    formatted_output(F, labelvalue);
    for (const auto &op : fields) formatted_output(F, m[op]);
    F << endl;
  }
   ExpvOutput(const string &fn, map<string, t_expv> &_m, list<string> _fields) : m(_m), fields(std::move(_fields)) {
    F.open(fn);
    field_numbers();
    field_names();
  }
};

unique_ptr<ExpvOutput> custom, customfdm;

// Prepare "td" for output: write header with parameters and
// a line with field names.
void open_Ftd(ofstream &Ftd) {
  Ftd.open(FN_TD);
  outfield::width = P::width_td;
  outfield::prec  = P::prec_td;
  TD::save_header(Ftd);
}

void open_files(speclist &sl, BaseSpectrum &spec, SPECTYPE spectype, axis a) {
  const string fn = spec.prefix + "_" + spectype->name() + "_dens_" + spec.name; // no suffix (.dat vs. .bin)
  switch (a) {
    case axis::RealFreq: spec.spec = make_shared<SpectrumRealFreq>(spec.name, fn, spectype); break;
    case axis::Temp: spec.spec = make_shared<SpectrumTemp>(spec.name, fn, spectype); break;
    case axis::Matsubara: spec.spec = make_shared<SpectrumMatsubara>(spec.name, fn, spectype, spec.mt); break;
    case axis::Matsubara2: spec.spec = make_shared<SpectrumMatsubara2>(spec.name, fn, spectype, spec.mt); break;
    default: my_assert_not_reached();
  }
  spec.spectype = spectype;
  spec.a        = a;
  nrglog('c', "Spectrum " << spec.fullname() << " -> " << fn);
  sl.push_back(spec);
}

// open_files_spec() opens the output files and establishes the data structures
// for storing spectral information.
void open_files_spec(speclist &sl, BaseSpectrum &spec) {
  if (spec.prefix == "gt") {
    if (nrgrun) open_files(sl, spec, make_shared<SPEC_GT>(), axis::Temp);
    return;
  }
  if (spec.prefix == "i1t") {
    if (nrgrun) open_files(sl, spec, make_shared<SPEC_I1T>(), axis::Temp);
    return;
  }
  if (spec.prefix == "i2t") {
    if (nrgrun) open_files(sl, spec, make_shared<SPEC_I2T>(), axis::Temp);
    return;
  }
  if (spec.prefix == "chit") {
    if (nrgrun) open_files(sl, spec, make_shared<SPEC_CHIT>(), axis::Temp);
    return;
  }
  // If we did not return from this funciton by this point, what we
  // are computing is the spectral function. There are several
  // possibilities in this case, all of which may be enabled at the
  // same time.
  if (nrgrun && P::finite) open_files(sl, spec, make_shared<SPEC_FT>(), axis::RealFreq);
  if (nrgrun && P::finitemats) open_files(sl, spec, make_shared<SPEC_FTmats>(), axis::Matsubara);
  if (dmnrgrun && P::dmnrg) open_files(sl, spec, make_shared<SPEC_DMNRG>(), axis::RealFreq);
  if (dmnrgrun && P::dmnrgmats) open_files(sl, spec, make_shared<SPEC_DMNRGmats>(), axis::Matsubara);
  if (dmnrgrun && P::cfs) open_files(sl, spec, make_shared<SPEC_CFS>(), axis::RealFreq);
  if (dmnrgrun && P::cfsgt) open_files(sl, spec, make_shared<SPEC_CFSgt>(), axis::RealFreq);
  if (dmnrgrun && P::cfsls) open_files(sl, spec, make_shared<SPEC_CFSls>(), axis::RealFreq);
  if (dmnrgrun && P::fdm) open_files(sl, spec, make_shared<SPEC_FDM>(), axis::RealFreq);
  if (dmnrgrun && P::fdmgt) open_files(sl, spec, make_shared<SPEC_FDMgt>(), axis::RealFreq);
  if (dmnrgrun && P::fdmls) open_files(sl, spec, make_shared<SPEC_FDMls>(), axis::RealFreq);
  if (dmnrgrun && P::fdmmats) open_files(sl, spec, make_shared<SPEC_FDMmats>(), axis::Matsubara);
}

void open_files_spec3(speclist &sl, BaseSpectrum &spec) {
  if (dmnrgrun && P::fdm && P::v3mm) // both options, fdm and v3mm
    open_files(sl, spec, make_shared<SPEC_FDM_v3mm>(), axis::Matsubara2);
}

namespace oprecalc {
  /* The following lists hold the names of operators which need to be
   recomputed. The default behavior is to recompute all the operators
   that are required to calculate the requested spectral densities, see
   function open_files(). In addition, singlet operators are always
   recomputed in the first NRG run, so that we can calculate the
   expectation values. In addition, if fdmexpv=true, the singlet operators
   are also recomputed in the second run if fdmexpvn=-1. */
  set<string> s, p, g, d, v, t, q, ot;

  void clear() {
    s.clear();
    p.clear();
    g.clear();
    d.clear();
    v.clear();
    t.clear();
    q.clear();
    ot.clear();
  }

  void report(ostream &F, const string &name, const set<string> &x) {
    F << name << "=[";
    for (const auto &i : x) F << i << ' ';
    F << "]" << endl;
  }

  void report(ostream &F = cout) {
    F << "Computing the following operators:" << endl;
    report(F, "s", s);
    report(F, "p", p);
    report(F, "g", g);
    report(F, "d", d);
    report(F, "v", v);
    report(F, "t", t);
    report(F, "q", q);
    report(F, "ot", ot);
  }

  bool do_s(const string &name) {
    if (nrgrun) return true;                               // for expectation values
    if (P::fdmexpv && STAT::N <= P::fdmexpvn) return true; // Calculate <O> using FDM algorithm
    return s.count(name);
  }

  bool do_g(const string &name) {
    if (nrgrun) return true;                               // for expectation values
    if (P::fdmexpv && STAT::N <= P::fdmexpvn) return true; // Calculate <O> using FDM algorithm
    return g.count(name);
  }

  bool do_p(const string &name) { return p.count(name); }
  bool do_d(const string &name) { return d.count(name); }
  bool do_v(const string &name) { return v.count(name); }
  bool do_t(const string &name) { return t.count(name); }
  bool do_q(const string &name) { return q.count(name); }
  bool do_ot(const string &name) { return ot.count(name); }

#define LOOPOVER(set1, set2, job) for (const auto &op1 : set1) for (const auto &op2 : set2) { job; } // NOLINT
#define LOOPOVER3(set1, set2, set3, job) for (const auto &op1 : set1) for (const auto &op2 : set2) for (const auto &op3 : set3) { job; } // NOLINT

  void OPENSPEC(const CustomOp::value_type &op1, const CustomOp::value_type &op2, const string_token &stringtoken, speclist &spectra, const string &prefix,
                set<string> &rec1, set<string> &rec2, matstype mt, int Spin = 0) {
    const string spname = SDNAME(NAME(op1), NAME(op2), Spin);
    if (stringtoken.find(spname)) {
      BaseSpectrum spec(op1.second, op2.second);
      spec.name   = spname;
      spec.prefix = prefix;
      spec.mt     = mt;
      spec.spin   = Spin;
      open_files_spec(spectra, spec);
      rec1.insert(NAME(op1));
      rec2.insert(NAME(op2));
    }
  }

  void OPENSPEC3(const CustomOp::value_type &op1, const CustomOp::value_type &op2, const CustomOp::value_type &op3, const string_token &stringtoken,
                 speclist &spectra, const string &prefix, set<string> &rec1, set<string> &rec2, set<string> &rec3, matstype mt) {
    const string spname = NAME(op1) + "-" + NAME(op2) + "-" + NAME(op3);
    if (stringtoken.find(spname)) {
      BaseSpectrum spec(op1.second, op2.second, op3.second);
      spec.name   = spname;
      spec.prefix = prefix;
      spec.mt     = mt;
      open_files_spec3(spectra, spec);
      rec1.insert(NAME(op1));
      rec2.insert(NAME(op2));
      rec3.insert(NAME(op3));
    }
  }

  // Reset lists of operators which need to be iterated.
  void reset_operator_lists_and_open_spectrum_files(const IterInfo &a) {
    oprecalc::clear();
    // Correlators (singlet operators of all kinds).
    string_token sts(P::specs);
    LOOPOVER(a.ops, a.ops, OPENSPEC(op1, op2, sts, spectraS, "corr", s, s, matstype::bosonic));
    LOOPOVER(a.opsp, a.opsp, OPENSPEC(op1, op2, sts, spectraS, "corr", p, p, matstype::bosonic));
    LOOPOVER(a.opsg, a.opsg, OPENSPEC(op1, op2, sts, spectraS, "corr", g, g, matstype::bosonic));
    LOOPOVER(a.ops, a.opsg, OPENSPEC(op1, op2, sts, spectraS, "corr", s, g, matstype::bosonic));
    LOOPOVER(a.opsg, a.ops, OPENSPEC(op1, op2, sts, spectraS, "corr", g, s, matstype::bosonic));
    // Global susceptibilities (global singlet operators).
    string_token stchit(P::specchit);
    LOOPOVER(a.ops, a.ops, OPENSPEC(op1, op2, stchit, spectraCHIT, "chit", s, s, matstype::bosonic));
    LOOPOVER(a.ops, a.opsg, OPENSPEC(op1, op2, stchit, spectraCHIT, "chit", s, g, matstype::bosonic));
    LOOPOVER(a.opsg, a.ops, OPENSPEC(op1, op2, stchit, spectraCHIT, "chit", g, s, matstype::bosonic));
    LOOPOVER(a.opsg, a.opsg, OPENSPEC(op1, op2, stchit, spectraCHIT, "chit", g, g, matstype::bosonic));
    // Dynamic spin susceptibilities (triplet operators).
    string_token stt(P::spect);
    LOOPOVER(a.opt, a.opt, OPENSPEC(op1, op2, stt, spectraT, "spin", t, t, matstype::bosonic));
    string_token stot(P::specot);
    LOOPOVER(a.opot, a.opot, OPENSPEC(op1, op2, stot, spectraOT, "orbspin", ot, ot, matstype::bosonic));
    const int varmin = (Sym->isfield() ? -1 : 0);
    const int varmax = (Sym->isfield() ? +1 : 0);
    // Spectral functions (doublet operators).
    string_token std(P::specd);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      LOOPOVER(a.opd, a.opd, OPENSPEC(op1, op2, std, spectraD, "spec", d, d, matstype::fermionic, SPIN));
    string_token stgt(P::specgt);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      LOOPOVER(a.opd, a.opd, OPENSPEC(op1, op2, stgt, spectraGT, "gt", d, d, matstype::fermionic, SPIN));
    string_token sti1t(P::speci1t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      LOOPOVER(a.opd, a.opd, OPENSPEC(op1, op2, sti1t, spectraI1T, "i1t", d, d, matstype::fermionic, SPIN));
    string_token sti2t(P::speci2t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      LOOPOVER(a.opd, a.opd, OPENSPEC(op1, op2, sti2t, spectraI2T, "i2t", d, d, matstype::fermionic, SPIN));
    // Spectral functions (quadruplet operators).
    string_token stq(P::specq);
    LOOPOVER(a.opq, a.opq, OPENSPEC(op1, op2, stq, spectraQ, "specq", q, q, matstype::fermionic));
    // Vertex functions
    string_token stv3(P::specv3);
    LOOPOVER3(a.opd, a.opd, a.ops, OPENSPEC3(op1, op2, op3, stv3, spectraV3, "specv3", d, d, s, matstype::fb));
    oprecalc::report();
    cout << endl << "Computing the following spectra:" << endl;
    for (const auto &i : allspectra) i->about();
  }
} // namespace oprecalc

// Open output files and write headers
void open_output_files(const IterInfo &iterinfo) {
  nrglog('@', "@ open_output_files()");
  // We dump all energies to separate files for NRG and DM-NRG runs.
  // This is a very convenient way to check if both runs produce the
  // same results.
  if (P::dumpenergies) Fenergies.open((nrgrun ? FN_ENERGIES_NRG : FN_ENERGIES_DMNRG));
  if (nrgrun) {
    open_Ftd(Ftd);
    if (P::dumpannotated) Fannotated.open(FN_ANNOTATED);
  }
  list<string> ops;
  for (const auto &op : iterinfo.ops) ops.push_back(NAME(op));
  for (const auto &op : iterinfo.opsg) ops.push_back(NAME(op));
  if (nrgrun)
    custom = make_unique<ExpvOutput>(FN_CUSTOM, STAT::expv, ops);
  else if (dmnrgrun && P::fdmexpv) 
    customfdm = make_unique<ExpvOutput>(FN_CUSTOMFDM, STAT::fdmexpv, ops);
  oprecalc::reset_operator_lists_and_open_spectrum_files(iterinfo);
}

// Close files. This has to be called explicitly, because there can be two
// separate runs (NRG and DMNRG). No error checking is done, nor do we test
// if the files are actually open.
void close_output_files() {
  if (nrgrun) {
    Fenergies.close();
    Ftd.close();
    Fannotated.close();
    custom.reset();
  }
  if (dmnrgrun && P::fdmexpv) customfdm.reset();
}

// DM-NRG: initialization of the density matrix -----------------------------

// Calculate grand canonical partition function at current NRG energy
// shell. This is not the same as the true partition function of the
// full problem! Instead this is the Z_N that is used to initialize
// the density matrix, i.e. rho = 1/Z_N \sum_{l} exp{-beta E_l} |l;N>
// <l;N|.  calc_grand_canonical_Z() is also used to calculate
// STAT::Zft, that is used to compute the spectral function with the
// conventional approach, as well as STAT::Zgt for G(T) calculations,
// STAT::Zchit for chi(T) calculations.
double calc_grand_canonical_Z(const DiagInfo &diag, double temperature = P::T) {
  bucket ZZ;
  LOOP_const(diag, is) for (const auto &x : EIGEN(is).value) ZZ += mult(INVAR(is)) * exp(-x * (STAT::scale / temperature));
  my_assert(ZZ >= 1.0);
  return ZZ;
}

// Calculate rho_N, the density matrix at the last NRG iteration. It is
// normalized to 1. Note: in CFS approach, we consider all states in the
// last iteration to be "discarded".
// For the details on the full Fock space approach see:
// F. B. Anders, A. Schiller, Phys. Rev. Lett. 95, 196801 (2005).
// F. B. Anders, A. Schiller, Phys. Rev. B 74, 245113 (2006).
// R. Peters, Th. Pruschke, F. B. Anders, Phys. Rev. B 74, 245114 (2006).
void init_rho(const DiagInfo &diag, DensMatElements &rho) {
  nrglog('@', "@ init_rho()");
  const double ZZ = calc_grand_canonical_Z(diag);
  rho.clear();
  LOOP_const(diag, is) {
    // The size of the matrix is determined by NRSTATES().
    const size_t dim = NRSTATES(is);
    const Invar I    = INVAR(is);
    rho[I]           = Matrix(dim, dim);
    rho[I].clear();
    for (size_t i = 0; i < dim; i++) rho[I](i, i) = exp(-EIGEN(is).value(i) * STAT::scale / P::T) / ZZ;
  }
  const double Tr = trace(rho);
  my_assert(num_equal(Tr, 1.0, 1e-8)); // NOLINT
}

// Calculate the shell-N statistical sum as used in the FDM algorithm.
double calc_ZnD_one(size_t N, double T) {
  double ZZ = 0;
  for (const auto &j : dm[N]) {
    // Determine which states count as discarded in the FDM sense
    const size_t min = (LAST_ITERATION(N) ? 0 : j.second.kept);
    const size_t max = j.second.total;
    for (size_t i = min; i < max; i++) {
      const double Eabs = j.second.absenergy[i] - STAT::GSenergy;
      my_assert(Eabs >= 0.0);
      const double betaE = Eabs / T;
      ZZ += mult(j.first) * exp(-betaE);
    }
  }
  nrglog('w', "ZnD[" << N << "]=" << HIGHPREC(ZZ));
  return ZZ;
}

/*** Truncation ***/

// Determine the number of states to be retained.
// Returns Emax - the highest energy to still be retained.
t_eigen highest_retained_energy(const STDEVEC &energies) {
  my_assert(energies.front() == 0.0); // check for the subtraction of Egs
  const size_t totalnumber = energies.size();
  size_t nrkeep;
  if (P::keepenergy <= 0.0) {
    nrkeep = P::keep;
  } else {
    double keepenergy = P::keepenergy;
    // We add 1 for historical reasons. We thus keep states with E<=Emax,
    // and one additional state which has E>Emax.
    nrkeep = 1 + count_if(begin(energies), end(energies), [=](double e) { return e <= keepenergy; });
    nrkeep = CLIP(nrkeep, size_t(P::keepmin), size_t(P::keep));
  }
  // Check for near degeneracy and ensure that the truncation occurs in a
  // "gap" between nearly-degenerate clusters of eigenvalues.
  if (P::safeguard > 0.0) {
    size_t cnt_extra = 0;
    while (nrkeep < totalnumber && (energies[nrkeep] - energies[nrkeep - 1]) <= P::safeguard && cnt_extra < P::safeguardmax) {
      nrkeep++;
      cnt_extra++;
    }
    if (cnt_extra) debug("Safeguard: keep additional " << cnt_extra << " states");
  }
  nrkeep = CLIP(nrkeep, size_t(1), totalnumber);
  return energies[nrkeep - 1];
}

// This function doesn't really perform the truncation in the sense
// of deleting matrix elements. It just computes the number of states
// to keep for each invariant subspace. Truncation is performed by
// nrg_trucate_perform().
void nrg_truncate_prepare(DiagInfo &diag) {
  nrglog('@', "@ nrg_truncate_prepare()");
  t_eigen Emax      = highest_retained_energy(STAT::energies);
  size_t nrkept     = 0; // counter of states actually retained
  size_t nrkeptmult = 0; // counter of states actually retained,
                         // taking into account the multiplicity
  LOOP(diag, is) {
    const Invar I             = INVAR(is);
    const size_t multiplicity = mult(INVAR(is));
    const size_t nrc          = NRSTATES(is); // number of (calculated) states in the subspace
    my_assert(nrc > 0);
    // Count the number of elements to keep
    size_t count = count_if(begin(EIGEN(is).value), end(EIGEN(is).value), [=](double e) { return e <= Emax; });
    my_assert(count <= nrc);
    // In CFS=true, we keep all states in the last iteration.
    if (cfs_flags() && LAST_ITERATION() && !P::lastalloverride) count = nrc;
    // Exception: in the last iteration of the density-matrix run, we
    // may override strategy=kept and keep all (computed) states if
    // lastall=true. This is only relevant if CFS=false and
    // DMNRG=true.
    if (dmnrgrun && LAST_ITERATION() && P::lastall) count = nrc;
    if (count == nrc && EIGEN(is).value(nrc - 1) != Emax) {
      /* We determined that all calculated states in this invariant
        subspace will be retained (and that perhaps additional should
        have been computed). RMAX() gives actual dimensionality of
        the invariant subspace (i.e. the maximal number of
        eigenpairs), while the variable nrc given by NRSTATES() is
        the number of actually calculated states, which is only equal
        to the dimensionality if in the computational scheme used all
        the states are retained. */
      const size_t truedim = RMAX(is);
      if (nrc < truedim) NRG::notenough = true;
    }
    diag[I].truncate_prepare(count);
    nrkept += count;
    nrkeptmult += count * multiplicity;
  }
  nrgdump3(Emax, nrkept, nrkeptmult) << endl;
}

// Calculate partial statistical sums, ZnD, and the grand canonical Z
// (STAT::ZZ). calc_ZnD() must be called before the second NRG run.
void calc_ZnD() {
  // Partial statistical sums, computed with respect to absolute excitation
  // energies! This implies that the absolute energies have to be used
  // in the calc_generic_FDM() function, since 1/ZnD are used as prefactors.
  STAT::ZnD = std::vector<double>(P::Nlen, 0.0);
  for (size_t N = P::Ninit; N < P::Nlen; N++) // here size_t, because Ninit>=0
    STAT::ZnD[N] = calc_ZnD_one(N, P::T);
  // Note: for ZBW, Nlen=Nmax+1. For Ninit=Nmax=0, index 0 will thus be included here.
  // Grand-canonical partition function (over all shells!!), i.e. \sum
  // \exp(-beta E) for all (!!) states in the complete Fock space, where
  // E is the absolute energy.
  bucket Z;
  for (size_t N = P::Ninit; N < P::Nlen; N++) Z += pow(double(P::combs), int(P::Nlen - N - 1)) * STAT::ZnD[N];
  STAT::ZZ = Z;
  cout << "Grand canonical Z= " << HIGHPREC(Z) << endl;
  cout << "Zft (last shell)=  " << HIGHPREC(STAT::Zft) << endl;
  // Weights for specific shells, w_n
  STAT::wn = std::vector<double>(P::Nlen, 0.0);
  STAT::wnfactor = std::vector<double>(P::Nlen, 0.0);
  bucket sumwn;
  for (size_t N = P::Ninit; N < P::Nlen; N++) {
    // This is w_n defined after Eq. (8) in the WvD paper.
    const double wn = pow(double(P::combs), int(P::Nlen - N - 1)) * STAT::ZnD[N] / Z;
    STAT::wn[N] = wn;
    sumwn += wn;
    const double w = pow(double(P::combs), int(P::Nlen - N - 1)) / Z;
    STAT::wnfactor[N] = w; // These ratios enter the terms for the spectral function.
  }
  for (size_t N = P::Ninit; N < P::Nlen; N++) 
    nrglog('w', "wn[" << N << "]=" << HIGHPREC(STAT::wn[N]));
  nrglog('w', "sumwn=" << sumwn << " sumwn-1=" << sumwn - 1.0);
  for (size_t N = P::Ninit; N < P::Nlen; N++)
    nrglog('w', "wnfactor[" << N << "]=" << HIGHPREC(STAT::wnfactor[N]));
  // Check the sum-rule.
  my_assert(num_equal(sumwn, 1.0));
}

// Actually truncate matrices. Drop subspaces with no states kept.
void nrg_truncate_perform(DiagInfo &diag) {
  for (auto &i : diag) EIGEN(i).truncate_perform(); // Truncate subspace to appropriate size
}

// scaled = true -> output scaled energies (i.e. do not multiply
// by the rescale factor)
inline t_eigen scaled_energy(t_eigen e, bool scaled = true, bool absolute = false) {
  return e * (scaled ? 1.0 : STAT::scale) + (absolute ? STAT::totalenergy : 0.0);
}

/* Store (eigenvalue/quantum numbers) pairs at given NRG iteration.
 Eigenenergies are useful for determining RG flows and checking for
 possible phase transitions between various different ground states. */
void dump_annotated(const DiagInfo &diag, bool scaled = true, bool absolute = false) {
  std::vector<pair<t_eigen, Invar>> seznam;
  for (const auto &is : diag)
    for (const auto e : EIGEN(is).value) seznam.emplace_back(e, INVAR(is));
  sort(begin(seznam), end(seznam));
  size_t len = min(seznam.size(), size_t(P::dumpannotated));
  // If states are clustered, we dump the full cluster
  for (size_t i = len; i < seznam.size() - 1; i++) {
    // If the next state has an energy within P::grouptol, add it to
    // the list.
    if (my_fcmp(seznam[i].first, seznam[i - 1].first, P::grouptol) == 0)
      len++;
    else
      break;
  }
  my_assert(len <= seznam.size());
  Fannotated << setprecision(P::dumpprecision);
  if (P::dumpgroups) {
    // Group by degeneracies
    for (size_t i = 0; i < len;) { // i increased in the while loop below
      Fannotated << scaled_energy(seznam[i].first, scaled, absolute);
      std::vector<string> QNstrings;
      size_t total_degeneracy = 0; // Total number of levels (incl multiplicity)
      const size_t i0         = i;
      while (i < len && my_fcmp(seznam[i].first, seznam[i0].first, P::grouptol) == 0) {
        QNstrings.push_back(to_string(seznam[i].second));
        total_degeneracy += mult(seznam[i].second);
        i++;
      }
      sort(begin(QNstrings), end(QNstrings));
      for (const auto &j : QNstrings) Fannotated << " (" << j << ")";
      Fannotated << " [" << total_degeneracy << "]" << endl;
    }
  } else {
    for (const auto &i : seznam) Fannotated << scaled_energy(i.first, scaled, absolute) << " " << i.second << endl;
  }
  // Consecutive iterations are separated by an empty line
  Fannotated << endl;
}

/* recalc_singlet() recalculates irreducible matrix elements of a
 singlet operator (nold) and stores them in a new matrix (nnew).
 This implementation is generic for all the symmetry types! */
void recalc_singlet(const DiagInfo &diag, const MatrixElements &nold, MatrixElements &nnew, int parity) {
  std::vector<Recalc> recalc_table(P::combs);
  const InvarVec input = input_subspaces();
  if (Sym->islr())
    my_assert(parity == 1 || parity == -1);
  else
    my_assert(parity == 1);
  LOOP_const(diag, is1) {
    const Invar I1 = INVAR(is1);
    Invar Ip       = I1;
    if (parity == -1) Ip.InvertParity();
    for (size_t i = 1; i <= P::combs; i++) {
      Recalc r;
      r.i1 = r.ip = i;
      r.factor    = 1.0;
      Invar ancI  = I1;
      ancI.combine(input[i]);
      r.IN1 = r.INp = ancI;
      if (parity == -1) r.INp.InvertParity();
      recalc_table[i - 1] = r; // mind the -1 shift!
    }
    recalc_general(diag, nold, nnew, Ip, I1, &recalc_table[0], P::combs, Sym->InvarSinglet);
  } // loop over is
}

// Wrapper routine for recalculations. Called from nrg_recalculate_operators().
template <class RecalcFnc>
void recalc_common(RecalcFnc recalc_fnc, DiagInfo &dg, CustomOp::value_type &op, const string &tip, bool (*testfn)(const string &)) {
  if (testfn(NAME(op))) {
    TIME("recalc " + tip);
    nrglog('0', "Recalculate " << tip << " " << NAME(op));
    MatrixElements opstore;
    opstore.swap(op.second);
    op.second.clear();
    recalc_fnc(dg, opstore, op.second);
    if (tip == "g") Sym->recalc_global(dg, NAME(op), op.second);
  } else
    op.second.clear(); // save memory!
}

/* We trim the matrices containing the irreducible matrix elements of the
 operators to the sizes that are actually required in the next iterations.
 This saves memory and leads to better cache usage in recalc_general()
 recalculations. Note: this is only needed for strategy=all; copying is
 avoided for strategy=kept. */
void nrg_trim_matel(DiagInfo &dg, MatrixElements &op) {
  for (auto &m : op) {
    const Invar &I1 = m.first.first;
    const Invar &I2 = m.first.second;
    // Current matrix dimensions
    const auto size1 = m.second.size1();
    const auto size2 = m.second.size2();
    if (size1 == 0 || size2 == 0) continue;
    // Target matrix dimensions
    const auto nr1 = dg[I1].getnr();
    const auto nr2 = dg[I2].getnr();
    my_assert(nr1 <= size1 && nr2 <= size2);
    if (nr1 == size1 && nr2 == size2) // Trimming not necessary!!
      continue;
    matrix_range<Matrix> m2(m.second, range(0, nr1), range(0, nr2));
    Matrix m2new = m2;
    m.second.swap(m2new);
  }
}

void nrg_trim_op(DiagInfo &dg, CustomOp &allops) {
  for (auto &op : allops) nrg_trim_matel(dg, op.second);
}

void nrg_trim_matrices(DiagInfo &dg, IterInfo &a) {
  nrg_trim_op(dg, a.ops);
  nrg_trim_op(dg, a.opsp);
  nrg_trim_op(dg, a.opsg);
  nrg_trim_op(dg, a.opd);
  nrg_trim_op(dg, a.opt);
  nrg_trim_op(dg, a.opot);
  nrg_trim_op(dg, a.opq);
}

void nrg_clear_eigenvectors(DiagInfo &diag) {
  LOOP(diag, i)
  for (auto &j : i.second.blocks) j = Matrix(0, 0);
}

// Z_S is the appropriate statistical sum
void nrg_measure_singlet1(const DiagInfo &dg, const CustomOp::value_type &op, double Z_S) {
  const string opname = NAME(op);
  const t_expv expv = calc_trace_singlet(dg, op.second) / Z_S;
  STAT::expv[opname] = expv;
  cout << "<" << opname << ">=" << output_val(expv) << endl;
}

// Expectation values using FDM algorithm
void nrg_measure_singlet1_fdm(const DiagInfo &dg, const CustomOp::value_type &op) {
  const string opname   = NAME(op);
  const t_expv expv   = calc_trace_fdm_kept(dg, op.second);
  STAT::fdmexpv[opname] = expv;
  cout << "<" << opname << ">_fdm=" << output_val(expv) << endl;
}

// Measure thermodynamic expectation values of singlet operators
void nrg_measure_singlet(const DiagInfo &dg, const IterInfo &a) {
  nrglog('@', "@ nrg_measure_singlet()");
  const double Z_S = calc_Z(dg);
  for (const auto &op : a.ops) nrg_measure_singlet1(dg, op, Z_S);
  for (const auto &op : a.opsg) nrg_measure_singlet1(dg, op, Z_S);
  custom->field_values();
}

void nrg_measure_singlet_fdm(const DiagInfo &dg, const IterInfo &a) {
  if (STAT::N != P::fdmexpvn) return;
  nrglog('@', "@ nrg_measure_singlet_fdm()");
  for (const auto &op : a.ops) nrg_measure_singlet1_fdm(dg, op);
  for (const auto &op : a.opsg) nrg_measure_singlet1_fdm(dg, op);
  customfdm->field_values(P::T);
}

void check_operator_sumrules(const DiagInfo &diag, const IterInfo &a) {
  nrglog('@', "@ check_operator_sumrules()");
  // We check sum rules wrt some given spin (+1/2, by convention).
  // For non-spin-polarized calculations, this is irrelevant (0).
  const int SPIN = (Sym->isfield() ? 1 : 0);
  for (const auto &op : a.opd) doublet_check_norm(op, diag, SPIN);
  for (const auto &op : a.opq) quadruplet_check_norm(op, diag, 0);
}

// Recalculate operator matrix representations
ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
void nrg_recalculate_operators(DiagInfo &dg, IterInfo &a) { // XXX: DiagInfo should be const, but all recalc* files would need to be fixed first
  nrglog('@', "@ nrg_recalculate_operators()");
  for (auto &op : a.ops) recalc_common([](const auto &a, const auto &b, auto &c) { recalc_singlet(a, b, c, 1); }, dg, op, "s", oprecalc::do_s);
  for (auto &op : a.opsp) recalc_common([](const auto &a, const auto &b, auto &c) { recalc_singlet(a, b, c, -1); }, dg, op, "p", oprecalc::do_p);
  for (auto &op : a.opsg) recalc_common([](const auto &a, const auto &b, auto &c) { recalc_singlet(a, b, c, 1); }, dg, op, "g", oprecalc::do_g);
  for (auto &op : a.opd) recalc_common([](auto &a, auto &b, auto &c) { Sym->recalc_doublet(a, b, c); }, dg, op, "d", oprecalc::do_d);
  for (auto &op : a.opt) recalc_common([](auto &a, auto &b, auto &c) { Sym->recalc_triplet(a, b, c); }, dg, op, "t", oprecalc::do_t);
  for (auto &op : a.opot) recalc_common([](auto &a, auto &b, auto &c) { Sym->recalc_orb_triplet(a, b, c); }, dg, op, "ot", oprecalc::do_ot);
  for (auto &op : a.opq) recalc_common([](auto &a, auto &b, auto &c) { Sym->recalc_quadruplet(a, b, c); }, dg, op, "q", oprecalc::do_q);
}

// Calculate spectral densities
void nrg_spectral_densities(const DiagInfo &diag) {
  nrglog('@', "@ nrg_spectral_densities()");
  TIME("spec");
  for (auto &i : spectraS) calc_generic(i, diag, CorrelatorFactorFnc, TrivialCheckSpinFnc);
  for (auto &i : spectraCHIT) calc_generic(i, diag, CorrelatorFactorFnc, TrivialCheckSpinFnc);
  for (auto &i : spectraD) calc_generic(i, diag, SpecdensFactorFnc, SpecdensCheckSpinFnc);
  for (auto &i : spectraT) calc_generic(i, diag, SpinSuscFactorFnc, TrivialCheckSpinFnc);
  for (auto &i : spectraOT) calc_generic(i, diag, OrbSuscFactorFnc, TrivialCheckSpinFnc);
  for (auto &i : spectraQ) calc_generic(i, diag, SpecdensquadFactorFnc, TrivialCheckSpinFnc);
  for (auto &i : spectraGT) calc_generic(i, diag, SpecdensFactorFnc, SpecdensCheckSpinFnc);
  for (auto &i : spectraI1T) calc_generic(i, diag, SpecdensFactorFnc, SpecdensCheckSpinFnc);
  for (auto &i : spectraI2T) calc_generic(i, diag, SpecdensFactorFnc, SpecdensCheckSpinFnc);
  // no CheckSpinFnc!! One must use A_u_d, etc. objects for sym=QSZ.
  for (auto &i : spectraV3) calc_generic3(i, diag, SpecdensFactorFnc);
}

/* We calculate thermodynamic quantities before truncation to make better
 use of the available states. Here we compute quantities which are defined
 for all symmetry types. Other calculations are performed by calculate_TD
 member functions defined in symmetry.cc. */
void nrg_calculate_TD(DiagInfo &diag, double additional_factor = 1.0) {
  nrglog('@', "@ nrg_calculate_TD()");
  TD::T = STAT::Teff;
  // Rescale factor for energies. The energies are expressed in
  // units of omega_N, thus we need to appropriately rescale them to
  // calculate the Boltzmann weights at the temperature scale Teff
  // (Teff=scale/betabar).
  double rescale_factor = P::betabar * additional_factor;
  bucket Z, E, E2; // Statistical sum, Tr[beta H], Tr[(beta H)^2]
  LOOP(diag, is) {
    bucket sumZ, sumE, sumE2;
    for (const auto &x : EIGEN(is).value) {
      const double betaE = rescale_factor * x;
      const double expo  = exp(-betaE);
      sumE += betaE * expo;
      sumE2 += sqr(betaE) * expo;
      sumZ += expo;
    }
    // Take into account the multiplicity
    const int multip = mult(INVAR(is));
    Z += multip * sumZ;
    E += multip * sumE;
    E2 += multip * sumE2;
  }
  STAT::Z = Z;
  // Energy, energy squared, heat capacity, free energy, entropy
  TD::E = STAT::E = E / Z;                   // beta <H>
  TD::E2 = STAT::E2 = E2 / Z;                // beta^2 <H^2>
  TD::C = STAT::C = STAT::E2 - sqr(STAT::E); // C/k_B=beta^2(<H^2>-<H>^2)
  TD::F = STAT::F = -log(Z);                 // F/(k_B T)=-ln(Z)
  TD::S = STAT::S = STAT::E - STAT::F;       // S/k_B=beta<H>+ln(Z)
  // Call quantum number specific calculation routine
  Sym->calculate_TD(diag, rescale_factor);
  if (nrgrun) TD::save_TD_quantities(Ftd);
}

void nrg_calculate_spectral_and_expv(const DiagInfo &diag, const IterInfo &iterinfo) {
  nrglog('@', "@ nrg_calculate_spectral_and_expv()");
  // Zft is used in the spectral function calculations using the
  // conventional approach. We calculate it here, in order to avoid
  // recalculations later on.
  STAT::Zft = calc_grand_canonical_Z(diag);
  if (string(P::specgt) != "" || string(P::speci1t) != "" || string(P::speci2t) != "") {
    double GTtemperature = P::gtp * STAT::scale;
    STAT::Zgt            = calc_grand_canonical_Z(diag, GTtemperature);
  }
  if (string(P::specchit) != "") {
    double CHITtemperature = P::chitp * STAT::scale;
    STAT::Zchit            = calc_grand_canonical_Z(diag, CHITtemperature);
  }
  if (dmnrgrun && !LAST_ITERATION()) {
    loadRho(STAT::N, FN_RHO, rho);
    if (P::checkrho) check_trace_rho(rho); // Check if Tr[rho]=1, i.e. the normalization
    if (P::fdm) loadRho(STAT::N, FN_RHOFDM, rhoFDM);
  }
  nrg_spectral_densities(diag);
  if (nrgrun) nrg_measure_singlet(diag, iterinfo);
  if (dmnrgrun && P::fdmexpv) nrg_measure_singlet_fdm(diag, iterinfo);
}

// Perform calculations of physical quantities. Called prior to NRG
// iteration (if calc0=true) and after each NRG step.
void nrg_perform_measurements(DiagInfo &diag) {
  nrglog('@', "@ nrg_perform_measurements()");
  nrg_calculate_TD(diag);
  if (nrgrun && P::dumpannotated) dump_annotated(diag, P::dumpscaled, P::dumpabs);
}

/* Make a list of subspaces for the new iteration. Generic implementation:
 use the quantum number differences in array In[] (obtained by a call to
 function input_subspaces), make a list of all possible subspaces and
 remove duplicates. */
void nrg_make_subspaces_list(list<Invar> &subspaces) {
  for (const auto &ii : diagprev)
    if (NRSTATES(ii)) {
      const Invar I = INVAR(ii);
      InvarVec input = input_subspaces(); // make a new copy of subspaces list
      for (size_t i = 1; i <= P::combs; i++) {
        input[i].inverse(); // IMPORTANT!
        input[i].combine(I);
        if (Sym->Invar_allowed(input[i])) subspaces.push_back(input[i]);
      }
    }
  subspaces.sort();
  subspaces.unique();
}

/* Define recalculation strategy
 all: Recalculate using all vectors
 kept: Recalculate using vectors kept after truncation
 VERY IMPORTANT: Override in the case of CFS (in the second run) */
bool do_recalc_kept() { return string(P::strategy) == "kept" && !(cfs_flags() && dmnrgrun) && !P::ZBW; }

bool do_recalc_all() { return !do_recalc_kept() && !P::ZBW; }

bool do_no_recalc() { return P::ZBW; }

t_eigen EigenvaluePrev(const Invar &I, size_t r) { return diagprev[I].value(r); }

Matrix nrg_prepare_task_for_diag(const Invar &I) {
  nrglog('@', "@ nrg_prepare_task_for_diag()");
  const InvarVec &anc = iterinfo.ancestors[I];
  const Rmaxvals &qq  = qsrmax[I];
  const size_t dim    = qq.total();
  nrglog('i', endl << "Subspace (" << I << ") dim=" << dim); // skip a line
  logancestors(I, anc, qq);
  Matrix h(dim, dim);
  h.clear();
  double scalefactor = (!P::substeps ? sqrt(P::Lambda) : pow(P::Lambda, 1. / (2. * P::channels))); // NOLINT
  // H_{N+1}=\lambda^{1/2} H_N+\xi_N (hopping terms)
  for (size_t i = 1; i <= P::combs; i++) {
    const size_t offset = qq.offset(i);
    for (size_t r = 0; r < qq.rmax(i); r++) h(offset + r, offset + r) = scalefactor * EigenvaluePrev(anc[i], r);
  }
  // Symmetry-type-specific matrix initialization steps.
  Sym->makematrix(h, qq, I, anc);
  if (logletter('m')) dump_matrix(h);
  return h;
}

void nrg_diagonalisations_OpenMP() {
  nrglog('@', "@ nrg_diagonalisations_OpenMP()");
  nrglog('(', "OpenMP diag");
  size_t nr    = NRG::tasks.size();
  size_t itask = 0;
  // cppcheck-suppress unreadVariable symbolName=nth
  int nth      = P::diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
  for (itask = 0; itask < nr; itask++) {
    Invar I  = NRG::tasks[itask];
    Matrix h = nrg_prepare_task_for_diag(I);
    int thid = omp_get_thread_num();
#pragma omp critical
    { nrglog('(', "Diagonalizing " << I << " size=" << h.size1() << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
    Eigen e = diagonalise(h);
#pragma omp critical
    { diag[I] = e; }
  }
}

#ifdef NRG_MPI
const int TAG_HELLO          = 1;
const int TAG_EXIT           = 2;
const int TAG_DIAG           = 3;
const int TAG_SYNC           = 4;
const int TAG_MATRIX         = 5;
const int TAG_INVAR          = 6;
const int TAG_EIGEN          = 7;
const int TAG_MATRIX_SIZE    = 8;
const int TAG_MATRIX_LINE    = 9;
const int TAG_EIGEN_INT      = 10;
const int TAG_EIGEN_VEC      = 11;
const int TAG_EIGEN_RMAXVALS = 12;

void mpi_sync_params();

void check_status(mpi::status &status) {
  if (status.error()) {
    cout << "MPI communication error. rank=" << mpiw->rank() << endl;
    mpienv->abort(1);
  }
}

// NOTE: MPI is limited to message size of 2GB (or 4GB). For big
// problems we thus need to send objects line by line.

//#define MPI_WHOLEMATRIX
#define MPI_LINEBYLINE

#ifdef MPI_WHOLEMATRIX
#define mpi_send_matrix mpi_send_matrix_wholematrix
#define mpi_receive_matrix mpi_receive_matrix_wholematrix
#define mpi_send_eigen mpi_send_eigen_whole
#define mpi_receive_eigen mpi_receive_eigen_whole
#endif

#ifdef MPI_LINEBYLINE
#define mpi_send_matrix mpi_send_matrix_linebyline
#define mpi_receive_matrix mpi_receive_matrix_linebyline
#define mpi_send_eigen mpi_send_eigen_linebyline
#define mpi_receive_eigen mpi_receive_eigen_linebyline
#endif

void mpi_send_matrix_wholematrix(int dest, Matrix &m) { mpiw->send(dest, TAG_MATRIX, m); }

void mpi_receive_matrix_wholematrix(int source, Matrix &m) {
  mpi::status status;
  status = mpiw->recv(source, TAG_MATRIX, m);
  check_status(status);
}

void mpi_send_matrix_linebyline(int dest, const Matrix &m) {
  size_t size1 = m.size1();
  mpiw->send(dest, TAG_MATRIX_SIZE, size1);
  size_t size2 = m.size2();
  mpiw->send(dest, TAG_MATRIX_SIZE, size2);
  nrglog('M', "Sending matrix of size " << size1 << " x " << size2 << " line by line to " << dest);
  for (size_t i = 0; i < size1; i++) {
    ublas::vector<t_matel> vec  = matrix_row<const Matrix>(m, i);
    mpiw->send(dest, TAG_MATRIX_LINE, vec);
  }
}

void mpi_receive_matrix_linebyline(int source, Matrix &m) {
  mpi::status status;
  size_t size1;
  status = mpiw->recv(source, TAG_MATRIX_SIZE, size1);
  check_status(status);
  size_t size2;
  status = mpiw->recv(source, TAG_MATRIX_SIZE, size2);
  check_status(status);
  nrglog('M', "Receiving matrix of size " << size1 << " x " << size2 << " line by line from " << source);
  m = Matrix(size1, size2);
  for (size_t i = 0; i < size1; i++) {
    ublas::vector<t_matel> vec;
    status = mpiw->recv(source, TAG_MATRIX_LINE, vec);
    check_status(status);
    my_assert(vec.size() == size2);
    matrix_row<Matrix>(m, i) = vec;
  }
}

void mpi_send_eigen_whole(int dest, const Eigen &eig) { mpiw->send(dest, TAG_EIGEN, eig); }

void mpi_receive_eigen_whole(int source, Eigen &eig) {
  mpi::status status;
  status = mpiw->recv(source, TAG_EIGEN, eig);
  check_status(status);
}

void mpi_send_eigen_linebyline(int dest, const Eigen &eig) {
  Eigen eigmock; // empty Eigen
  mpiw->send(dest, TAG_EIGEN, eigmock);
  nrglog('M', "Sending eigen from " << mpiw->rank() << " to " << dest);
  mpiw->send(dest, TAG_EIGEN_INT, eig.nr);
  mpiw->send(dest, TAG_EIGEN_INT, eig.rmax);
  mpiw->send(dest, TAG_EIGEN_INT, eig.nrpost);
  mpiw->send(dest, TAG_EIGEN_VEC, eig.value);
  mpiw->send(dest, TAG_EIGEN_VEC, eig.absenergy); // XXX: necessary?
  mpi_send_matrix_linebyline(dest, eig.matrix0);
}

void mpi_receive_eigen_linebyline(int source, Eigen &eig) {
  nrglog('M', "Receiving eigen from " << source << " on " << mpiw->rank());
  mpi::status status;
  Eigen eigmock;
  status = mpiw->recv(source, TAG_EIGEN, eigmock);
  check_status(status);
  status = mpiw->recv(source, TAG_EIGEN_INT, eig.nr);
  check_status(status);
  status = mpiw->recv(source, TAG_EIGEN_INT, eig.rmax);
  check_status(status);
  status = mpiw->recv(source, TAG_EIGEN_INT, eig.nrpost);
  check_status(status);
  status = mpiw->recv(source, TAG_EIGEN_VEC, eig.value);
  check_status(status);
  status = mpiw->recv(source, TAG_EIGEN_VEC, eig.absenergy); // XXX
  check_status(status);
  mpi_receive_matrix_linebyline(source, eig.matrix0);
}

// Read results from a slave process.
Invar read_from(int source) {
  nrglog('M', "Reading results from " << source);
  mpi::status readst;
  Eigen eig;
  mpi_receive_eigen(source, eig);
  Invar Irecv;
  readst = mpiw->recv(source, TAG_INVAR, Irecv);
  check_status(readst);
  nrglog('M', "Received results for subspace " << Irecv << " [nr=" << eig.getnr() << ", dim=" << eig.getrmax() << "]");
  // Some consistency checks
  my_assert(eig.matrix0.size1() == eig.getnr());
  my_assert(eig.matrix0.size2() == eig.getrmax());
  my_assert(eig.matrix0.size1() <= eig.matrix0.size2());
  my_assert(eig.getrmax() == qsrmax[Irecv].total());
  diag[Irecv] = eig;
  return Irecv;
}

void nrg_diagonalisations_MPI() {
  nrglog('@', "@ nrg_diagonalisations_MPI()");
  mpi_sync_params(); // Synchronise parameters
  list<Invar> todo(NRG::tasks); // List of all the remaining tasks
  list<Invar> done; // List of finished tasks.
  // List of the available computation nodes (including the master,
  // which is always at the very beginnig of the deque).
  deque<int> nodes;
  for (size_t i = 0; i < mpiw->size(); i++) nodes.push_back(i);
  size_t nrtasks = NRG::tasks.size();
  nrglog('M', "nrtasks=" << nrtasks << " nrnodes=" << mpiw->size());
  while (!todo.empty()) {
    my_assert(!nodes.empty());
    // i is the node to which the next job will be scheduled
    int i;
    if (todo.size() == 1) {
      // If a single task is left undone, do it on the master node
      // to avoid the unnecessary network copying.
      i = 0;
    } else {
      i = nodes.back();
      nodes.pop_back();
    }
    Invar I;
    if (i == 0) {
      // On master, we take short jobs from the end of the list.
      I = todo.back();
      todo.pop_back();
    } else {
      // On slaves, we take long jobs from the beginning of the
      // list.
      I = todo.front();
      todo.pop_front();
    }
    Matrix h = nrg_prepare_task_for_diag(I);
    nrglog('M',
           "Scheduler: job " << I << " (dim=" << h.size1() << ")"
                             << " on node " << i);
    if (i == 0) {
      // On master, diagonalize immediately.
      diag[I] = diagonalise(h);
      nodes.push_back(0);
      done.push_back(I);
    } else {
      mpiw->send(i, TAG_DIAG, 0);
      mpi_send_matrix(i, h);
      mpiw->send(i, TAG_INVAR, I);
    }
    // Check for terminated jobs
    boost::optional<mpi::status> status;
    while (status = mpiw->iprobe(mpi::any_source, TAG_EIGEN)) {
      nrglog('M', "Receiveing results from " << status->source());
      Invar Irecv = read_from(status->source());
      done.push_back(Irecv);
      // The node is now available for new tasks!
      nodes.push_back(status->source());
    }
  }
  // Keep reading results sent from the slave processes until all
  // tasks have been completed.
  size_t nrdone = done.size();
  while (nrdone != nrtasks) {
    nrglog('M', "So far: " << nrdone << "/" << nrtasks);
    mpi::status status = mpiw->probe(mpi::any_source, TAG_EIGEN);
    Invar Irecv        = read_from(status.source());
    done.push_back(Irecv);
    nrdone++;
  }
}
#endif

// Build matrix H(ri;r'i') in each subspace and diagonalize it
void nrg_diagonalisations() {
  nrglog('@', "@ nrg_diagonalisations()");
  // This needs to be called here, because class Timing is not
  // thread-safe.
  TIME("diag");
  // Call init() again here, because NRG::diagratio might have
  // changed!
  sP.init();
  diag.clear();
#ifdef NRG_MPI
  nrg_diagonalisations_MPI();
#else
  nrg_diagonalisations_OpenMP();
#endif
}

// Determine the list of invariant subspaces in which diagonalisations need
// to be performed.
void nrg_determine_tasks(map<Invar, InvarVec> &ancestors) {
  nrglog('@', "@ nrg_determine_tasks()");
  // Make a list of all subspaces to consider.
  list<Invar> subspaces;
  nrg_make_subspaces_list(subspaces);
  // Auxiliary information: ancestor subspaces and their dimensions.
  qsrmax.clear();
  ancestors.clear();
  // Container holding all the subspaces that appear in the new
  // iteration.
  NRG::tasks.clear();
  for (const auto &I : subspaces) {
    // Determine which subspaces contribute to the Hamiltonian being built
    InvarVec input = input_subspaces();
    for (size_t i = 1; i <= P::combs; i++) input[i].combine(I); // In is the list of differences wrt I
    ancestors[I] = input;
    // Determine the range(s) of index r
    qsrmax[I].determine_ranges(I, input);
    // nr is actually the size of the Hamiltonian submatrix!
    const size_t nr = qsrmax[I].total();
    // Note that NRG::tasks is not the same list as 'subspaces',
    // since some possible subspaces may have dimension 0 and thus
    // do not really exist.
    if (nr) NRG::tasks.push_back(I);
  }
}

// Sort NRG::tasks in the decreasing order of the submatrix size. As
// a side effect, we compute some statistics about matrix sizes.
void sort_task_list() {
  std::vector<pair<size_t, Invar>> tasks_with_sizes;
  for (const auto &i : NRG::tasks) tasks_with_sizes.emplace_back(qsrmax[i].total(), i);
  // Sort in the *decreasing* order!
  sort(rbegin(tasks_with_sizes), rend(tasks_with_sizes));
  auto nr       = tasks_with_sizes.size();
  auto min_size = tasks_with_sizes.back().first;
  auto max_size = tasks_with_sizes.front().first;
  cout << "Stats: nr=" << nr << " min=" << min_size << " max=" << max_size << endl;
  // Report matrix sizes
  if (logletter('S'))
    for (const auto &i : tasks_with_sizes) cout << "size(" << i.second << ")=" << i.first << endl;
  // Update the task list NRG::tasks with the sorted list of subspaces
  transform(begin(tasks_with_sizes), end(tasks_with_sizes), begin(NRG::tasks), [](const auto &p) { return p.second; });
}

// Recalculate irreducible matrix elements for Wilson chains.
// Called from nrg_after_diag().
void nrg_recalc_f(DiagInfo &diag, Opch &opch) {
  nrglog('@', "@ nrg_recalc_f()");
  TIME("recalc f");
  if (!P::substeps) {
    for (size_t i = 0; i < P::channels; i++)
      for (size_t j = 0; j < P::perchannel; j++) opch[i][j].clear(); // Clear all channels
    // In principle, one could also use recalc_doublet() function
    // to simplify the code. OTOH, recalc_irreduc() is probably
    // better because it does not accumulate floating point
    // round-off errors.
    Sym->recalc_irreduc(diag, opch);
  } else {
    size_t Ntrue, M;
    tie(Ntrue, M) = get_Ntrue_M(STAT::N);
    for (size_t i = 0; i < P::channels; i++) {
      if (i == M) {
        for (size_t j = 0; j < P::perchannel; j++) opch[M][j].clear(); // Clear channel M
        Sym->recalc_irreduc_substeps(diag, opch, M);
      } else {
        for (size_t j = 0; j < P::perchannel; j++) {
          MatrixElements &f = opch[i][j];
          MatrixElements opstore;
          opstore.swap(f);
          f.clear();
          Sym->recalc_doublet(diag, opstore, f);
        }
      }
    }
  }
}

void nrg_dump_f(const Opch &opch) {
  cout << endl;
  for (size_t i = 0; i < P::channels; i++)
    for (size_t j = 0; j < P::perchannel; j++) {
      cout << "<f> dump, i=" << i << " j=" << j << endl;
      dump_matrix_elements(opch[i][j]);
    }
  cout << endl;
}

// The absenergy[] values are shifted so that the ground state corresponds
// to zero. This is required in the FDM approach for calculating the
// spectral functions. This is different from subtract_groundstate_energy().
// Called from nrg_do_diag() when diag is loaded from a stored file during
// the second pass of the NRG iteration. Note that the values in absenergy[]
// are not yet shifted when ZnD is calculated in calc_ZnD_one, but it is
// already shifted when calc_generic_FDM(), calc_trace_fdm_discarded/kept()
// are called.
void shift_abs_energies(DiagInfo &diag) {
  LOOP(diag, i)
  for (auto &x : EIGEN(i).absenergy) x -= STAT::GSenergy;
}

// Used in evaluation of vertex functions to speed up the computation.
void calc_boltzmann_factors(DiagInfo &diag) {
  LOOP(diag, i) {
    const size_t len = EIGEN(i).absenergy.size();
    EIGEN(i).boltzmann.resize(len);
    for (size_t j = 0; j < len; j++) EIGEN(i).boltzmann[j] = exp(-EIGEN(i).absenergy[j] / P::T);
  }
}

/* NRG diagonalisation driver: calls nrg_diagionalisations() or
 load_transformations(), as necessary, and performs the truncation. All other
 calculations are done in nrg_after_diag(). Called from nrg_iterate(). */
void nrg_do_diag() {
  nrglog('@', "@ nrg_do_diag()");
  infostring();
  show_coefficients();
  nrg_determine_tasks(iterinfo.ancestors);
  sort_task_list();
  NRG::diagratio = P::diagratio;
  do {
    NRG::notenough = false;
    if (nrgrun) {
      if (!(P::resume && int(STAT::N) <= P::laststored))
        nrg_diagonalisations(); // compute in first run
      else
        load_transformations(STAT::N, diag); // or read from disk
    }
    if (dmnrgrun) {
      load_transformations(STAT::N, diag); // read from disk in second run
      // IMPORTANT: subtract the absolute (!) GS energy in the
      // abs_energy vector. The overall (all shells, all invariant
      // subspaces) lowest abs_energy will thus be equal to zero.
      shift_abs_energies(diag);
      calc_boltzmann_factors(diag);
      if (P::removefiles) remove_transformation_files(STAT::N);
    }
    time_mem::ms.check("diag");
    find_groundstate(diag);
    subtract_groundstate_energy(diag);
    copy_sort_energies(diag, STAT::energies);
    find_clusters(STAT::energies, P::fixeps, STAT::cluster_mapping);
    bool recopy_required = fix_splittings(diag);
    // fix_splittings() returns true if any changes had been made. If so,
    // we determine the clusters again (just to see what the result is).
    if (recopy_required) {
      copy_sort_energies(diag, STAT::energies);
      find_clusters(STAT::energies, P::fixeps, STAT::cluster_mapping);
    }
    nrg_truncate_prepare(diag);
    time_mem::ms.check("trunc");
    if (NRG::notenough) {
      cout << "Insufficient number of states computed." << endl;
      if (P::restart) {
        NRG::diagratio = min(NRG::diagratio * P::restartfactor, 1.0);
        cout << endl
             << "Restarting this iteration step. "
             << "diagratio=" << NRG::diagratio << endl
             << endl;
      }
    }
  } while (nrgrun && P::restart && NRG::notenough);
}

using mapSD = map<string, double>;
std::vector<mapSD> td_data;

// Store all the results from the TD calculation in the form of a map
// (string -> double).
void store_td() {
  mapSD td;
  for (const auto &i : allfields) td[i->name()] = i->rawvalue();
  td_data.push_back(td);
}

// Absolute energies. Must be called in the first NRG run after
// STAT::totalenergy has been updated, but before
// store_transformations(). These energies are initially not
// referrenced to absolute 0. This is done in the second NRG run in
// shift_abs_energies().
void calc_abs_energies(DiagInfo &diag) {
  nrglog('@', "@ calc_abs_energies()");
  LOOP(diag, is) {
    EIGEN(is).absenergy = EIGEN(is).value;
    for (auto &x : EIGEN(is).absenergy) x = STAT::totalenergy + x * STAT::scale;
  }
}

// Perform processing after a successful NRG step.
// At function call:
// - diag contains all information about the eigenstates.
// - STAT::Egs had been computed
// Also called from doZBW() as a final step.
void nrg_after_diag(IterInfo &iterinfo) {
  nrglog('@', "@ nrg_after_diag()");
  // Contribution to the total energy.
  STAT::totalenergy += STAT::Egs * STAT::scale;
  cout << "Total energy=" << setprecision(std::numeric_limits<double>::max_digits10) << STAT::totalenergy << "  Egs=" << STAT::Egs << endl;
  if (nrgrun) calc_abs_energies(diag);
  if (P::dm && nrgrun) {
    // Store eigenenergies and eigenvectors in all subspaces.
    // NOTE: we store before the truncation!
    if (!(P::resume && int(STAT::N) <= P::laststored)) store_transformations(STAT::N, diag);
  }
  // Logging of ALL states (not only those that remain after truncation)
  if (P::dumpenergies) dumptofile(diag, Fenergies);
  // Measurements are performed before the truncation!
  nrg_perform_measurements(diag);
  // Consistency checks
  LOOP_const(diag, i) {
    my_assert(EIGEN(i).matrix0.size1() <= EIGEN(i).matrix0.size2());
    if (!P::ZBW) my_assert(EIGEN(i).matrix0.size2() == qsrmax[INVAR(i)].total());
  }
  if (!P::ZBW) {
    split_in_blocks(diag);
    time_mem::ms.check("split");
  }
  if (do_recalc_all()) { // Either ...
    nrg_recalculate_operators(diag, iterinfo);
    nrg_calculate_spectral_and_expv(diag, iterinfo);
  }
  // Actual truncation occurs at this point
  if (!P::ZBW) nrg_truncate_perform(diag);
  // Store information about subspaces and states
  size_t nrall  = 0;
  size_t nrkept = 0;
  LOOP_const(diag, i) {
    const Invar I  = INVAR(i);
    dm[STAT::N][I] = DimSub(NRSTATES(i), RMAX(i), qsrmax[I], EIGEN(i).value, EIGEN(i).absenergy);
    nrall += RMAX(i);
    nrkept += NRSTATES(i);
  }
  double ratio = double(nrkept) / nrall;
  cout << "Kept: " << nrkept << " out of " << nrall << ", ratio=" << setprecision(3) << ratio << endl;
  if (!LAST_ITERATION()) {
    nrg_recalc_f(diag, iterinfo.opch);
    if (P::dump_f) nrg_dump_f(iterinfo.opch);
  }
  if (do_recalc_kept()) { // ... or ...
    nrg_recalculate_operators(diag, iterinfo);
    nrg_calculate_spectral_and_expv(diag, iterinfo);
  }
  if (do_no_recalc()) { // ... or this
    nrg_calculate_spectral_and_expv(diag, iterinfo);
  }
  if (P::checksumrules) check_operator_sumrules(diag, iterinfo);
  time_mem::ms.check("recalc");
  if (!P::ZBW) {
    // Free up memory that contains information we no longer need
    nrg_trim_matrices(diag, iterinfo);
    nrg_clear_eigenvectors(diag);
    time_mem::ms.check("trim");
    diagprev.swap(diag); // IMPORTANT: we need to retain the eigenenergies!
  }
  // Store TD data (all outfields)
  store_td();
}

// Perform one iteration step
void nrg_iterate(IterInfo &iterinfo) {
  nrg_do_diag();
  nrg_after_diag(iterinfo);
  time_mem::memory_time_brief_report();
}

void docalc0ht(unsigned int extra_steps) {
  for (int i = -int(extra_steps); i <= -1; i++) {
    STAT::set_N(int(P::Ninit) - 1 + i);
    double E_rescale_factor = pow(P::Lambda, i / 2.0); // NOLINT
    nrg_calculate_TD(diagprev, E_rescale_factor);
  }
}

// Perform calculations with quantities from 'data' file
void docalc0(const IterInfo &iterinfo) {
  nrglog('@', "@ docalc0()");
  STAT::set_N(int(P::Ninit) - 1); // in the usual case with Ninit=0, this will result in N=-1
  cout << endl << "Before NRG iteration";
  cout << " (N=" << STAT::N << ")" << endl;
  nrg_perform_measurements(diagprev);
  nrg_calculate_spectral_and_expv(diagprev, iterinfo);
  // Logging of ALL states (prior to truncation!)
  if (P::dumpenergies) dumptofile(diagprev, Fenergies);
  if (P::checksumrules) check_operator_sumrules(diagprev, iterinfo);
}

// doZBW() takes the place of nrg_iterate() called from
// nrg_main_loop() in the case of zero-bandwidth calculation.
// Thus it replaces nrg_do_diag() + nrg_after_diag().
void doZBW(IterInfo &iterinfo) {
  cout << endl << "Zero bandwidth calculation" << endl;
  // TRICK: scale will be that for N=Ninit-1, but STAT::N=Ninit.
  STAT::set_N(int(P::Ninit) - 1);
  STAT::N = P::Ninit; // this is a hack!
  // begin nrg_do_diag() equivalent
  if (nrgrun) diag = diagprev;
  if (dmnrgrun) {
    load_transformations(STAT::N, diag);
    shift_abs_energies(diag);
    calc_boltzmann_factors(diag); // !!
    remove_transformation_files(STAT::N);
  }
  find_groundstate(diag);
  subtract_groundstate_energy(diag);
  copy_sort_energies(diag, STAT::energies); // required in nrg_truncate_prepare()
  nrg_truncate_prepare(diag);               // determine # of kept and discarded states
  // end nrg_do_diag() equivalent
  nrg_after_diag(iterinfo);
}

// ****************************  Main NRG loop ****************************

void nrg_main_loop(IterInfo &iterinfo) {
  nrglog('@', "@ nrg_main_loop()");

  // N denotes the order of the Hamiltonian. N=0 corresponds to H_0, i.e.
  // the initial Hamiltonian (cf. Krishna-Murty I)
  for (size_t N = P::Ninit; N < P::Nmax; N++) { // here size_t, because Ninit>=0. Note that N is int.
    if (nrgrun && P::forcestop == int(N)) exit1("*** Stop request at iteration " << P::forcestop);
    STAT::set_N(N);
    nrg_iterate(iterinfo);
  }
}

// Estimate the Kondo temperature (Wilson's definition). This
// calculation only makes sense for the SIAM and the Kondo model and
// even then only if Tmin << TKW, B=0, and Sz^2 is computed.
void calc_TKW() {
  const size_t len = td_data.size();
  if (!(outfield_exists("<Sz^2>") && len > 2)) return;
  const double T0chiT0 = td_data[len - 1]["<Sz^2>"];
  const double y       = T0chiT0 + 0.07;
  // Linear interpolation
  double T = 0.0;
  for (size_t i = len - 1; i >= 1; i--) {
    const double x0 = td_data[i]["T"];          // lower T
    const double x1 = td_data[i - 1]["T"];      // higher T
    const double y0 = td_data[i]["<Sz^2>"];     // lower Tchi(T)
    const double y1 = td_data[i - 1]["<Sz^2>"]; // higher Tchi(T)
    if (y0 <= y && y <= y1) {
      T = x0 + (x1 - x0) / (y1 - y0) * (y - y0);
      break;
    }
  }
  cout << "TKW=" << T << endl;
}

// Total number of states (symmetry taken into account)
size_t count_states(const DiagInfo &dg) {
  size_t states = 0;
  LOOP_const(dg, i) states += Sym->mult(INVAR(i)) * NRSTATES(i);
  return states;
}

// Count non-empty subspaces
size_t count_subspaces(const DiagInfo &dg) {
  size_t subspaces = 0;
  LOOP_const(dg, i) if (NRSTATES(i)) subspaces++;
  return subspaces;
}

// Dump information about states.
void states_report(const DiagInfo &dg, ostream &fout = cout) {
  fout << "Number of invariant subspaces: " << count_subspaces(dg) << endl;
  LOOP_const(dg, is) if (NRSTATES(is)) fout << "(" << INVAR(is) << ") " << NRSTATES(is) << " states: " << is.second.value << endl;
  fout << "Number of states (multiplicity taken into account): " << count_states(dg) << endl << endl;
}

void start_run(IterInfo &iterinfo) {
  nrglog('@', "@ start_run()");
  states_report(diagprev);
  open_output_files(iterinfo);
  // If setting calc0 is set to true, a calculation of TD quantities
  // is performed before starting the NRG iteration.
  if (nrgrun && P::calc0 && !P::ZBW) {
    docalc0ht(P::tdht);
    docalc0(iterinfo);
  }
  if (P::ZBW) doZBW(iterinfo); // in both NRG and DMNRG runs
  nrg_main_loop(iterinfo);
  close_output_files();
  cout << endl << "** Iteration completed." << endl << endl;
}

// === AFTER THE NRG ITERATION HAD COMPLETED =================================

// Processing performed both after NRG and after DM runs.
void finalize_common() {
  nrglog('@', "@ finalize_common()");
  TIME("broaden");
  for (auto &i : allspectra) i->clear(); // processing happens in the destructor
}

// Save a dump of all subspaces, with dimension info, etc.
void dump_subspace_information() {
  ofstream O(FN_SUBSPACES);
  for (size_t N = P::Ninit; N < P::Nmax; N++) {
    O << "Iteration " << N << endl;
    O << "len_dm=" << dm[N].size() << endl;
    for (const auto &i : dm[N])
      O << "I=" << i.first << " len=" << i.second.eigenvalue.size() << " kept=" << i.second.kept << " total=" << i.second.total << endl;
    O << endl;
  }
}

// Called after the first NRG run.
void finalize_nrg() {
  nrglog('@', "@ finalize_nrg()");
  finalize_common();
  cout << endl << "Total energy: " << HIGHPREC(STAT::totalenergy) << endl;
  // True ground state energy is just the value of totalenergy at the end
  // of the iteration. This is the energy of the lowest energy state in the
  // last iteration. All other states (incl. from previous shells)
  // obviously have higher energies.
  STAT::GSenergy = STAT::totalenergy;
  calc_TKW();
  if (P::fdm) calc_ZnD();
  if (P::dumpsubspaces) dump_subspace_information();
}

// Called after the second NRG run.
void finalize_dmnrg() {
  nrglog('@', "@ finalize_dmnrg()");
  finalize_common();
  // These two should match if value_raw and value vectors are
  // handled correctly. GSenergy was computed in the first NRG run,
  // while totalenergy is recomputed in the second DMNRG run.
  my_assert(num_equal(STAT::GSenergy, STAT::totalenergy));
}

/************************ MAIN ****************************************/

void print_about_message(ostream &s) {
  s << "NRG Ljubljana - (c) rok.zitko@ijs.si" << endl;
  s << "Timestamp: " << __TIMESTAMP__ << endl;
  s << "Compiled on " << __DATE__ << " at " << __TIME__ << endl << endl;
}

// Called immediately after the symmetry type is determined from the data file. This ensures that Invar can be parsed correctly.
void set_symmetry(const string &sym_string) {
  if (all_syms.count(sym_string) != 1) exit1("Unknown symmetry " << sym_string);
  cout << "SYMMETRY TYPE: " << sym_string << endl;
  Sym = all_syms[sym_string];
  TD::init();
  my_assert(P::channels > 0); // must be set at this point
  Sym->init();
  Sym->set_combs(P::combs);
  Sym->set_channels(P::channels);
  Sym->set_substeps(P::substeps);
  Sym->load(); // This call actually initializes the data structures.
  Sym->report();
}

void prep_run(RUNTYPE t)
{
  STAT::runtype = t;
  read_data(iterinfo);
}

void start_calculation(IterInfo &iterinfo) {
  nrglog('@', "@ start_calculation()");
  prep_run(RUNTYPE::NRG);
  dm = AllSteps(P::Nlen);
  start_run(iterinfo);
  finalize_nrg();
  if (string(P::stopafter) == "nrg") exit1("*** Stopped after the first sweep.");
  if (!P::dm) return;
  init_rho(diagprev, rho);
  if (P::fdm) init_rho_FDM(rhoFDM, STAT::N);
  if (!P::ZBW) {
    calc_densitymatrix(rho);
    if (P::fdm) calc_fulldensitymatrix(rhoFDM);
  }
  if (string(P::stopafter) == "rho") exit1("*** Stopped after the DM calculation.");
  prep_run(RUNTYPE::DMNRG);
  start_run(iterinfo);
  finalize_dmnrg();
}

// Copy parameters to an object that is synchronised accross all
// processes. Warning: init() has to be called at the beginning of
// the program (after parsing the parameters in P), but also before
// each series of diagonalizations, because NRG::diagratio might have
// changed!
void sharedparam::init() {
  diagroutine = P::diagroutine;
  diagratio   = NRG::diagratio;
  dsyevrlimit = P::dsyevrlimit;
  logall      = P::logall;
  log         = P::log;
}

#ifdef NRG_MPI
int mpidebuglevel = 0;

void mpidebug(const string &str) {
  if (mpidebuglevel > 0) cout << "MPI process " << mpiw->rank() << " " << str << endl;
}

void mpi_sync_params() {
  // Synchronize global parameters
  if (mpiw->rank() == 0) {
    sP.init();
    for (size_t i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_SYNC, 0);
  }
  mpi::broadcast(*mpiw, sP, 0);
}
#endif

// What is the last iteration completed in the previous NRG runs?
void init_laststored() {
  if (P::resume) {
    P::laststored = -1;
    for (size_t N = P::Ninit; N < P::Nmax; N++) {
      const string fn = unitaryfn(N);
      ifstream F(fn.c_str()); // open file
      if (F.good())           // successful?
        P::laststored = N;
    }
    cout << "Last unitary file found: " << P::laststored << endl;
  }
}

// Master process does most of the i/o and passes calculations to the slaves.
void run_nrg_master() {
  cout << setprecision(std::numeric_limits<double>::max_digits10); // ensure no precision is lost
  read_parameters();
  validate_parameters();
  calculate_invariants();
  init_laststored();
  dump_parameters();
  sP.init(); // copy parameters for MPI
#ifdef NRG_MPI
  for (int i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_HELLO, 0);
#endif
  // IterInfo iterinfo; // TO DO: remove global variable
  start_calculation(iterinfo);
#ifdef NRG_MPI
  cout << "Master done. Terminating slave processes." << endl;
  for (int i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_EXIT, 0);
  cout << "Master exiting." << endl;
#endif
  if (P::done) { ofstream D("DONE"); } // Indicate completion by creating a flag file
  remove_workdir(); // Ok to make this call multiple times
}

#ifdef NRG_MPI
// Handle a diagonalisation request:
void slave_diag() {
  mpi::status status;
  const int MASTER = 0;
  // 1. receive the matrix and the subspace identification
  mpidebug("recv");
  Matrix m;
  mpi_receive_matrix(MASTER, m);
  Invar I;
  status = mpiw->recv(MASTER, TAG_INVAR, I);
  check_status(status);
  // 2. preform the diagonalisation
  mpidebug("diagonalise");
  Eigen eig = diagonalise(m);
  // 3. send back the results
  mpidebug("send");
  mpi_send_eigen(MASTER, eig);
  mpiw->send(MASTER, TAG_INVAR, I);
}

void run_nrg_slave() {
  set_new_handler(outOfMemory);
  cout << "MPI slave rank " << mpiw->rank() << endl;
  const int MASTER = 0;
  bool done        = false;
  while (!done) {
    boost::optional<mpi::status> probestatus;
    if (probestatus = mpiw->iprobe(MASTER, mpi::any_tag)) {
      // A message can be received.
      int task;
      mpi::status status;
      status = mpiw->recv(MASTER, mpi::any_tag, task);
      check_status(status);
      nrglog('M', "Slave " << mpiw->rank() << " received message with tag " << status.tag());
      switch (status.tag()) {
        case TAG_HELLO: mpidebug("ready"); break;
        case TAG_EXIT:
          mpidebug("exiting");
          done = true;
          break;
        case TAG_DIAG:
          mpidebug("diag");
          slave_diag();
          break;
        case TAG_SYNC:
          mpidebug("sync");
          mpi_sync_params();
          break;
        default: cout << "MPI error: unknown tag on " << mpiw->rank() << endl; break;
      } // switch
    } else {
      // No message received. We sleep for a while to reduce the
      // load on the computer. (OpenMPI "feature" workaround)
      usleep(100);
    }
  } // while(!done)
}
#else
void run_nrg_slave() {}
#endif
