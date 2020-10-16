/*
 "NRG Ljubljana" - Numerical renormalization group for multiple
 impurities and an arbitrary number of channels

 Copyright (C) 2005-2020 Rok Zitko

   This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
   License as published by the Free Software Foundation; either version 2 of the License, or (at your option) any
   later version.

   This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied
   warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
   details.

   You should have received a copy of the GNU General Public License along with this program; if not, write to the
   Free Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

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
#include "mp.h"

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
using t_expv = cmpl; // we allow the calculation of expectation values of non-Hermitian operators!
inline cmpl CONJ_ME(cmpl z) { return conj(z); }
#endif

using t_weight = cmpl; // spectral weight accumulators (complex in general)

enum class RUNTYPE { NRG, DMNRG };

#include "params.cc"
#include "outfield.cc"

// This is included in the library only. Should not be used if cblas library is available.
#ifdef CBLAS_WORKAROUND
 #define ADD_
 #include "cblas_globals.c"
 #include "cblas_dgemm.c"
 #include "cblas_zgemm.c"
 #include "cblas_xerbla.c"
#endif

inline const size_t MAX_NDX = 1000; // max index number

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
}

#ifdef NRG_MPI
mpi::environment *mpienv;
mpi::communicator *mpiw;
int myrank() { return mpiw->rank(); } // used in diag.h, time_mem.h
#else
int myrank() { return 0; }
#endif

// Quantum number types defined to enforce type checking
using Number = int;
using Ispin = int;
using Sspin = int;
using Tangmom = int;
using SZspin = int;

#include "invar.cc"

// NOTE: Row major is the C array format: A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], etc. The default in UBLAS is
// row major, while LAPACK routines expect column major matrices. Of course, this is of no concern for symmetric
// matrices. Default storage type is unbounded_array<T>.
//
// Thus, as always:
//  first index - row
//  second index - column
//
// when accessing columns - stride=m.size2()
// when accessing rows - stide=1
//
// Optimization rule: use stride 1 sequential access where possible. ublas default matrix storage is row major (i.e.
// C-like). The rule is "right index the same as inner loop variable".

using Matrix = ublas::matrix<t_matel>;

template<typename T> auto range0(const T b) { return boost::irange(T{0}, b); }
template<typename T> auto range1(const T b) { return boost::irange(T{1}, b+1); }

template <typename T>
  void save(boost::archive::binary_oarchive &oa, const ublas::matrix<T> &m) {
    oa << m.size1() << m.size2();
    for (const auto i : range0(m.size1()))
      oa << ublas::vector<T>(ublas::matrix_row<const ublas::matrix<T>>(m, i));
  }

template <typename T>
  void load(boost::archive::binary_iarchive &ia, ublas::matrix<T> &m) {
    size_t size1, size2;
    ia >> size1 >> size2;
    m = ublas::matrix<T>(size1, size2);
    for (const auto i : range0(size1)) {
      ublas::vector<T> vec;
      ia >> vec;
      ublas::matrix_row<ublas::matrix<T>>(m, i) = vec;
    }
  }

#include "numerics.h"

// Result of a diagonalisation: eigenvalues and eigenvectors
struct RawEigen {
  using EVEC = ublas::vector<t_eigen>;
  EVEC value_orig;     // eigenvalues as computed
  Matrix matrix;       // eigenvectors
  RawEigen() {}
  RawEigen(size_t nr, size_t dim) {
    my_assert(nr <= dim);
    value_orig.resize(nr);
    matrix.resize(nr, dim);
  }
  auto getnrcomputed() const { return value_orig.size(); } // number of computed eigenpairs
  auto getdim() const { return matrix.size2(); }           // valid also after the split_in_blocks_Eigen() call
};
  
// Augments RawEigen with the information about truncation and block structure of the eigenvectors.
class Eigen : public RawEigen {
 private:
  long nrpost = -1;  // number of eigenpairs after truncation (-1: keep all)
 public:
  EVEC value_zero;   // eigenvalues with Egs subtracted
  size_t getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }
  auto getnrstored() const  { return value_zero.size(); }                   // number of stored states
  auto getnrall() const { return getnrcomputed(); }                         // all = all computed
  auto getnrkept() const { return getnrpost(); }
  auto getnrdiscarded() const { return getnrcomputed()-getnrpost(); }
  auto all() const { return range0(getnrcomputed()); }                           // iterator over all states
  auto kept() const { return range0(getnrpost()); }                              // iterator over kept states
  auto discarded() const { return boost::irange(getnrpost(), getnrcomputed()); } // iterator over discarded states
  auto stored() const { return range0(getnrstored()); }                          // iterator over all stored states
  // NOTE: "absolute" energy means that it is expressed in the absolute energy scale rather than SCALE(N).
  EVEC absenergy;      // absolute energies
  EVEC absenergyG;     // absolute energies (0 is the absolute ground state of the system) [SAVED TO FILE]
  EVEC absenergyN;     // absolute energies (referenced to the lowest energy in the N-th step)
  // 'blocks' contains eigenvectors separated according to the invariant
  // subspace from which they originate. This separation is required for
  // using the efficient BLAS routines when performing recalculations of
  // the matrix elements.
  std::vector<Matrix> blocks;
  Eigen() : RawEigen() {}
  Eigen(size_t nr, size_t rmax) : RawEigen(nr, rmax) {}
  // Truncate to nrpost states.
  void truncate_prepare_subspace(size_t nrpost_) {
    nrpost = nrpost_;
    my_assert(nrpost <= getnrstored());
  }
  void truncate_perform() {
    for (auto &i : blocks) {
      my_assert(nrpost <= i.size1());
      i.resize(nrpost, i.size2());
    }
    value_zero.resize(nrpost);
  }
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors form an identity matrix. This is used to
  // represent the spectral decomposition in the eigenbasis itself.
  void diagonal(const EVEC &v) {
    value_orig = value_zero = v;
    matrix   = ublas::identity_matrix<t_eigen>(v.size());
  }
  void subtract_Egs(double Egs) {
    value_zero = value_orig;
    for (auto &x : value_zero) x -= Egs;
    my_assert(value_zero[0] >= 0);
  }
  void subtract_GS_energy(double GS_energy) {
    for (auto &x : absenergyG) x -= GS_energy;
    my_assert(absenergyG[0] >= 0);
  }
  void save(boost::archive::binary_oarchive &oa) const {
    // RawEigen
    oa << value_orig;
    ::save(oa, matrix);
    // Eigen
    oa << value_zero << nrpost << absenergy << absenergyG << absenergyN;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    // RawEigen
    ia >> value_orig;
    ::load(ia, matrix);
    // Eigen
    ia >> value_zero >> nrpost >> absenergy >> absenergyG >> absenergyN;
  } 
};

// Full information after diagonalizations (eigenspectra in all subspaces)
class DiagInfo : public std::map<Invar, Eigen> {
 public:
   explicit DiagInfo() {}
   DiagInfo(ifstream &fdata, const size_t nsubs, const Params &P) {
     for (const auto i : range1(nsubs)) {
       Invar I;
       fdata >> I;
       auto energies = read_vector<double>(fdata);
       if (!P.data_has_rescaled_energies && !P.absolute)
         energies /= P.SCALE(P.Ninit); // rescale to the suitable energy scale
       (*this)[I].diagonal(energies);
     }
     my_assert(size() == nsubs);
   }
   auto subspaces() const { return *this | boost::adaptors::map_keys; }
   auto eigs() const { return *this | boost::adaptors::map_values; }
   auto eigs() { return *this | boost::adaptors::map_values; }
   t_eigen find_groundstate() const {
     auto [Iground, eig] = *ranges::min_element(*this, [](const auto a, const auto b) { return a.second.value_orig(0) < b.second.value_orig(0); });
     auto Egs = eig.value_orig(0);
     return Egs;
   }
   void subtract_Egs(const t_eigen Egs) {
     ranges::for_each(this->eigs(), [Egs](auto &eig)       { eig.subtract_Egs(Egs); });
   }
   void subtract_GS_energy(const t_eigen GS_energy) {
     ranges::for_each(this->eigs(), [GS_energy](auto &eig) { eig.subtract_GS_energy(GS_energy); });
   }
   std::vector<t_eigen> sorted_energies() const {
     std::vector<t_eigen> energies;
     for (const auto &eig: this->eigs())
       energies.insert(energies.end(), eig.value_zero.begin(), eig.value_zero.end());
     return energies | ranges::move | ranges::actions::sort;
   }
   void dump_value_zero(ostream &F) const {
     for (const auto &[I, eig]: *this)
       F << "Subspace: " << I << std::endl << eig.value_zero << std::endl;
   }
   void truncate_perform() {
     for (auto &[I, eig] : *this) eig.truncate_perform(); // Truncate subspace to appropriate size
   }
   size_t size_subspace(const Invar &I) const {
     const auto f = this->find(I);
     return f != this->cend() ? f->second.getnrstored() : 0;
   }
   void clear_eigenvectors() {
     for (auto &eig : this->eigs())
       for (auto &m : eig.blocks) 
         m = Matrix(0, 0);
   }
   // Total number of states (symmetry taken into account)
   template <typename MF> auto count_states(MF && mult) const {
     return ranges::accumulate(*this, 0, [mult](auto n, const auto &x) { const auto &[I, eig] = x; return n + mult(I)*eig.getnrstored(); });
   }
   auto count_subspaces() const {    // Count non-empty subspaces
     return ranges::count_if(this->eigs(), [](const auto &eig) { return eig.getnrstored()>0; });
   }
   template <typename MF>
     void states_report(MF && mult) const {
       fmt::print("Number of invariant subspaces: {}\n", count_subspaces());
       for (const auto &[I, eig]: *this) 
         if (eig.getnrstored()) 
           fmt::print("({}) {} states: {}\n", I.str(), eig.getnrstored(), eig.value_orig);
       fmt::print("Number of states (multiplicity taken into account): {}\n\n", count_states(mult));
     }
   void save(const size_t N) const {
     const string fn = workdir.unitaryfn(N);
     ofstream MATRIXF(fn, ios::binary | ios::out);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
     boost::archive::binary_oarchive oa(MATRIXF);
     oa << this->size();
     for(const auto &[I, eig]: *this) {
       oa << I;
       eig.save(oa);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn)); // Check after each write.
     }
   }
   void load(const size_t N, const bool remove_files = false) {
     const string fn = workdir.unitaryfn(N);
     std::ifstream MATRIXF(fn, ios::binary | ios::in);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
     boost::archive::binary_iarchive ia(MATRIXF);
     size_t nr; // Number of subspaces
     ia >> nr;
     for (const auto cnt : range0(nr)) {
       Invar inv;
       ia >> inv;
       (*this)[inv].load(ia);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
     }
     if (remove_files) remove(fn);
   }
   explicit DiagInfo(const size_t N, const bool remove_files = false) { load(N, remove_files); }
};

class MatrixElements : public std::map<Twoinvar, Matrix> {
 public:
   MatrixElements() {}
   MatrixElements(ifstream &fdata, const DiagInfo &diag) {
     size_t nf; // Number of I1 x I2 combinations
     fdata >> nf;
     for (const auto i : range0(nf)) {
       Invar I1, I2;
       fdata >> I1 >> I2;
       if (const auto it1 = diag.find(I1), it2 = diag.find(I2); it1 != diag.end() && it2 != diag.end())
         read_matrix(fdata, (*this)[{I1, I2}], it1->second.getnrstored(), it2->second.getnrstored());
       else
         throw std::runtime_error("Corrupted input file.");
     }
     my_assert(this->size() == nf);
   }
   std::ostream &insertor(std::ostream &os) const { 
     for (const auto &[II, mat] : *this)
       os << "----" << II << "----" << endl << mat << endl;
     return os;
   }
};
std::ostream &operator<<(std::ostream &os, const MatrixElements &m) { return m.insertor(os); }

class DensMatElements : public std::map<Invar, Matrix> {
 public:
   template <typename MF>
     double trace(MF mult) const {
       return ranges::accumulate(*this, 0.0, [mult](double acc, const auto z) { const auto &[I, mat] = z; 
                                return acc + mult(I) * trace_real_nochecks(mat); });
     }
   void save(size_t N, const string &prefix) const {
     const string fn = workdir.rhofn(prefix, N);
     std::ofstream MATRIXF(fn, ios::binary | ios::out);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
     boost::archive::binary_oarchive oa(MATRIXF);
     oa << this->size();
     for (const auto &[I, mat] : *this) {
       oa << I;
       ::save(oa, mat);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn));  // Check each time
     }
     MATRIXF.close();
   }
   void load(size_t N, const string &prefix, const bool remove_files) {
     const string fn = workdir.rhofn(prefix, N);
     std::ifstream MATRIXF(fn, ios::binary | ios::in);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
     boost::archive::binary_iarchive ia(MATRIXF);
     size_t nr;
     ia >> nr;
     for (const auto cnt : range0(nr)) {
       Invar inv;
       ia >> inv;
       ::load(ia, (*this)[inv]);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));  // Check each time
     }
     MATRIXF.close();
     if (remove_files)
       if (remove(fn)) throw std::runtime_error(fmt::format("Error removing {}", fn));
   }
};

// Map of operators matrices
using CustomOp = std::map<std::string, MatrixElements>;

// Vector containing irreducible matrix elements of f operators.
using OpchChannel = std::vector<MatrixElements>;
// Each channel contains P.perchannel OpchChannel matrices.
class Opch : public std::vector<OpchChannel> {
 public:
   Opch() {}
   explicit Opch(const size_t nrch) { this->resize(nrch); }
   Opch(ifstream &fdata, const DiagInfo &diag, const Params &P) {
     this->resize(P.channels);
     for (const auto i : range0(size_t(P.channels))) {
       (*this)[i] = OpchChannel(P.perchannel);
       for (const auto j : range0(size_t(P.perchannel))) {
         char ch;
         size_t iread, jread;
         fdata >> ch >> iread >> jread;
         my_assert(ch == 'f' && i == iread && j == jread);
         (*this)[i][j] = MatrixElements(fdata, diag);
       }
     }
   }
   void dump() {
     std::cout << std::endl;
     for (const auto &&[i, ch] : *this | ranges::views::enumerate)
       for (const auto &&[j, mat] : ch | ranges::views::enumerate)
         std::cout << fmt::format("<f> dump, i={} j={}\n", i, j) << mat << endl;
     std::cout << std::endl;
   }
};

class Symmetry;

// Dimensions of the invariant subspaces |r,1>, |r,2>, |r,3>, etc. The name "rmax" comes from the maximal value of
// the index "r" which ranges from 1 through rmax.

class Rmaxvals {
 private:
   std::vector<size_t> values;
   shared_ptr<Symmetry> Sym;
 public:
   Rmaxvals() = default;
   Rmaxvals(const Invar &I, const InvarVec &In, const DiagInfo &diagprev, shared_ptr<Symmetry> Sym);
   auto combs() const { return values.size(); }
   auto rmax(const size_t i) const {
     my_assert(i < combs());
     return values[i];
   }
   auto exists(const size_t i) const {
     my_assert(i < combs());
     return values[i] > 0; 
   }
   auto offset(const size_t i) const {
     my_assert(i < combs());
     return ranges::accumulate(begin(values), begin(values) + i, size_t{0});
   }
   auto operator[](const size_t i) const { return rmax(i); }
   auto total() const { return ranges::accumulate(values, 0); } // total number of states
 private:
   friend ostream &operator<<(ostream &os, const Rmaxvals &rmax) {
     for (const auto &x : rmax.values) os << x << ' ';
     return os;
   }
   template <class Archive> void serialize(Archive &ar, const unsigned int version) { ar &values; }
   friend class boost::serialization::access;
};

class QSrmax : public std::map<Invar, Rmaxvals> {
 public:
   QSrmax() {}
   QSrmax(const DiagInfo &, shared_ptr<Symmetry>);
   // List of invariant subspaces in which diagonalisations need to be performed
   std::vector<Invar> task_list() const {
     std::vector<pair<size_t, Invar>> tasks_with_sizes;
     for (const auto &[I, rm] : *this)
       if (rm.total())
         tasks_with_sizes.emplace_back(rm.total(), I);
     ranges::sort(tasks_with_sizes, std::greater<>()); // sort in the *decreasing* order!
     auto nr       = tasks_with_sizes.size();
     auto min_size = tasks_with_sizes.back().first;
     auto max_size = tasks_with_sizes.front().first;
     cout << "Stats: nr=" << nr << " min=" << min_size << " max=" << max_size << endl;
     return tasks_with_sizes | ranges::views::transform( [](const auto &p) { return p.second; } ) | ranges::to<std::vector>();
   }
   void dump() const {
     for(const auto &[I, rm]: *this)
       std::cout << "rmaxvals(" << I << ")=" << rm << " total=" << rm.total() << std::endl;
   }
   auto at_or_null(const Invar &I) const {
     const auto i = this->find(I);
     return i == this->cend() ? Rmaxvals() : i->second;
   }
};

// Information about the number of states, kept and discarded, rmax, and eigenenergies. Required for the
// density-matrix construction.
template<typename M> struct DimSubGen {
  using EVEC = ublas::vector<M>;
  size_t kept  = 0;
  size_t total = 0;
  Rmaxvals rmax;
  Eigen eig;
  bool is_last = false;
  size_t min() const { return is_last ? 0 : kept; } // min(), max() return the range of D states to be summed over in FDM
  size_t max() const { return total; }
  auto all() const { return boost::irange(min(), max()); }
};

using DimSub = DimSubGen<t_eigen>;

// Full information about the number of states and matrix dimensions
// Example: dm[N].rmax[I] etc.
using Subs = map<Invar, DimSub>;

class AllSteps : public std::vector<Subs> {
 public:
   size_t Nbegin, Nend; // range of valid indexes
   AllSteps(size_t Nbegin, size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend ? Nend : 1); } // at least 1 for ZBW
   auto Nall() const { return boost::irange(Nbegin, Nend); }
   void dump_absenergyG(ostream &F) const {
     for (const auto N : Nall()) {
       F << std::endl << "===== Iteration number: " << N << std::endl;
       for (const auto &[I, ds]: this->at(N))
         F << "Subspace: " << I << std::endl << ds.eig.absenergyG << std::endl;
     }
   }
   void dump_all_absolute_energies(std::string filename = "absolute_energies.dat"s) {
     std::ofstream F(filename);
     this->dump_absenergyG(F);
   }
   // Save a dump of all subspaces, with dimension info, etc.
   void dump_subspaces(const std::string filename = "subspaces.dat"s) const {
     ofstream O(filename);
     for (const auto N : Nall()) {
       O << "Iteration " << N << std::endl;
       O << "len_dm=" << this->at(N).size() << std::endl;
       for (const auto &[I, DS] : this->at(N))
         O << "I=" << I << " kept=" << DS.kept << " total=" << DS.total << std::endl;
       O << std::endl;
     }
   }
   void shift_abs_energies(const double GS_energy) {
     for (const auto N : Nall())
       for (auto &ds : this->at(N) | boost::adaptors::map_values)
         ds.eig.subtract_GS_energy(GS_energy);
   }
   void store(const size_t ndx, const DiagInfo &diag, const QSrmax &qsrmax, const bool last) {
     my_assert(Nbegin <= ndx && ndx < Nend);
     for (const auto &[I, eig]: diag)
       (*this)[ndx][I] = { eig.getnrkept(), eig.getdim(), qsrmax.at_or_null(I), eig, last };
   }
};

class Step {
 private:
   // N denotes the order of the Hamiltonian. N=0 corresponds to H_0, i.e. the initial Hamiltonian
   int trueN; // "true N", sets the energy scale, it may be negative, trueN <= ndxN
   size_t ndxN; // "index N", iteration step, used as an array index, ndxN >= 0
   const Params &P; // reference to parameters (beta, T)
   
 public:
   RUNTYPE runtype; // NRG vs. DM-NRG run
   void set(int newN) {
     trueN = newN;
     ndxN = std::max(newN, 0);
   }
   void init() { set(P.Ninit); }
   Step(const Params &P_, RUNTYPE runtype_) : P(P_), runtype(runtype_) { init(); }
   void next() { trueN++; ndxN++; }
   size_t N() const { return ndxN; }
   size_t ndx() const { return ndxN; }
   double energyscale() const { return P.SCALE(trueN+1); } // current energy scale in units of bandwidth D
   double scale() const { // scale factor as used in the calculation
     return P.absolute ? 1.0 : energyscale();
   }
   double unscale() const { // 'unscale' parameter for dimensionless quantities
     return P.absolute ? energyscale() : 1.0;
   }
   double Teff() const { return energyscale()/P.betabar; }  // effective temperature for thermodynamic calculations
   double TD_factor() const { return P.betabar / unscale(); }
   double scT() const { return scale()/P.T; } // scT = scale*P.T, scaled physical temperature that appears in the exponents in spectral function calculations (Boltzmann weights)
   pair<size_t, size_t> NM() const {
     const size_t N = ndxN / P.channels;
     const size_t M = ndxN - N*P.channels; // M ranges 0..channels-1
     return {N, M};
   }
   void infostring() const {
     auto info = fmt::format(" ***** [{}] Iteration {}/{} (scale {}) ***** ", runtype == RUNTYPE::NRG ? "NRG"s : "DM"s, 
                             ndxN+1, int(P.Nmax), energyscale());
     info += P.substeps ? fmt::format(" step {} substep {}", NM().first+1, NM().second+1) : "";
     fmt::print(fmt::emphasis::bold, "\n{}\n", info);
   }
   void set_ZBW() {
     trueN = P.Ninit - 1; // if Ninit=0, trueN will be -1 (this is the only exceptional case)
     ndxN = P.Ninit;
   }
   // Return true if the spectral-function merging is to be performed at the current step
   bool N_for_merging() const {
     if (P.NN1) return true;
     if (P.NN2avg) return true;
     return P.NN2even ? IS_EVEN(ndxN) : IS_ODD(ndxN);
   }
   size_t firstndx() const { return P.Ninit; }
   size_t lastndx() const { return P.ZBW ? P.Ninit : P.Nmax-1; }
   // Return true if this is the first step of the NRG iteration
   bool first() const { return ndxN == firstndx(); }
   // Return true if N is the last step of the NRG iteration
   bool last(int N) const {
     return N == lastndx() || (P.ZBW && N == firstndx()); // special case!
   }
   bool last() const { return last(ndxN); }
   bool end() const { return ndxN >= P.Nmax; } // ndxN is outside the allowed range
   // NOTE: for ZBWcalculations, Ninit=0 and Nmax=0, so that first() == true and last() == true for ndxN=0.
   bool nrg() const { return runtype == RUNTYPE::NRG; }
   bool dmnrg() const { return runtype == RUNTYPE::DMNRG; }
   // Index 'n' of the last site in the existing chain, f_n (at iteration 'N'). The site being added is f_{n+1}. This
   // is the value that we use in building the matrix, cf. nrg-make_matrix-ISO.cc
   int getnn() const { return ndxN; }
};

// Namespace for storing various statistical quantities calculated during iteration.
class Stats {
 public:
   t_eigen Egs{};
   
   // ** Thermodynamic quantities
   double Z{};
   double Zft{};   // grand-canonical partition function (at shell n)
   double Zgt{};   // grand-canonical partition function for computing G(T)
   double Zchit{}; // grand-canonical partition function for computing chi(T)

   TD td;
   
   //  ** Expectation values
   map<string, t_expv> expv;    // expectation values of custom operators
   map<string, t_expv> fdmexpv; // Expectation values computed using the FDM algorithm
   
   // ** Energies
   // "total_energy" is the total energy of the ground state at the current iteration. This is the sum of all the 
   // zero state energies (eigenvalue shifts converted to absolute energies) for all the iteration steps so far.
   t_eigen total_energy{};
   // GS_energy is the energy of the ground states in absolute units. It is equal to the value of the variable
   // "total_energy" at the end of the iteration.
   t_eigen GS_energy{};
   std::vector<double> rel_Egs;        // Values of 'Egs' for all NRG steps.
   std::vector<double> abs_Egs;        // Values of 'Egs' (multiplied by the scale, i.e. in absolute scale) for all NRG steps.
   std::vector<double> energy_offsets; // Values of "total_energy" for all NRG steps.
   
   // Containers related to the FDM-NRG approach
   // ==========================================
   // Consult A. Weichselbaum, J. von Delft, PRL 99, 076402 (2007).
   vmpf ZnDG;                    // Z_n^D=\sum_s^D exp(-beta E^n_s), sum over **discarded** states at shell n
   vmpf ZnDN;                    // Z'_n^D=Z_n^D exp(beta E^n_0)=\sum_s^D exp[-beta(E^n_s-E^n_0)]
   std::vector<double> ZnDNd;    // 
   std::vector<double> wn;       // Weights w_n. They sum to 1.
   std::vector<double> wnfactor; // wn/ZnDG

   double ZZG{};                 // grand-canonical partition function with energies referred to the ground state energy

   double Z_fdm{};               // grand-canonical partition function (full-shell) at temperature T
   double F_fdm{};               // free-energy at temperature T
   double E_fdm{};               // energy at temperature T
   double C_fdm{};               // heat capacity at temperature T
   double S_fdm{};               // entropy at temperature T
   
   TD_FDM td_fdm;

   explicit Stats(const Params &P, const std::string filename_td = "td"s, const std::string filename_tdfdm = "tdfdm"s) : 
     td(P, filename_td), rel_Egs(MAX_NDX), abs_Egs(MAX_NDX), energy_offsets(MAX_NDX), 
     ZnDG(MAX_NDX), ZnDN(MAX_NDX), ZnDNd(MAX_NDX), wn(MAX_NDX), wnfactor(MAX_NDX), td_fdm(P, filename_tdfdm) {}
};

class ChainSpectrum;
class BaseSpectrum;
using spCS_t = shared_ptr<ChainSpectrum>;

// Wrapper class for NRG spectral-function algorithms
class Algo {
 private:
 public:
   const Params &P;
   Algo() = delete;
   Algo(const Algo&) = delete;
   explicit Algo(const Params &P) : P(P) {}
   virtual ~Algo() = default;
   virtual std::shared_ptr<ChainSpectrum> make_cs(const BaseSpectrum &) = 0;
   virtual void calc(const Step &step, const Eigen &, const Eigen &, const Matrix &, const Matrix &, 
                     const BaseSpectrum &, t_factor, std::shared_ptr<ChainSpectrum>, const Invar &,
                     const Invar &, const DensMatElements &, const Stats &stats) const {};
   virtual string name() = 0;
   virtual string merge() { return ""; }    // what merging rule to use
   virtual string rho_type() { return ""; } // what rho type is required
};

void dump_diagonal_matrix(const Matrix &m, const size_t max_nr, ostream &F)
{
  for (const auto r : range0(std::min(m.size1(), max_nr)))
    F << m(r,r) << ' ';
  F << std::endl;
}

void dump_diagonal_op(const std::string &name, const MatrixElements &n, const size_t max_nr, ostream &F) {
  F << "Diagonal matrix elements of operator " << name << std::endl;
  for (const auto &[II, mat] : n) {
    const auto & [I1, I2] = II;
    if (I1 == I2) {
      F << I1 << ": ";
      dump_diagonal_matrix(mat, max_nr, F);
    }
  }
}

// We trim the matrices containing the irreducible matrix elements of the operators to the sizes that are actually
// required in the next iterations. This saves memory and leads to better cache usage in recalc_general()
// recalculations. Note: this is only needed for strategy=all; copying is avoided for strategy=kept.
void trim_matel(const DiagInfo &diag, MatrixElements &op) {
  for (auto &[II, mat] : op) {
    const auto &[I1, I2] = II;
    // Current matrix dimensions
    const auto size1 = mat.size1();
    const auto size2 = mat.size2();
    if (size1 == 0 || size2 == 0) continue;
    // Target matrix dimensions
    const auto nr1 = diag.at(I1).getnrstored();
    const auto nr2 = diag.at(I2).getnrstored();
    my_assert(nr1 <= size1 && nr2 <= size2);
    if (nr1 == size1 && nr2 == size2) // Trimming not necessary!!
      continue;
    ublas::matrix_range<Matrix> m2(mat, ublas::range(0, nr1), ublas::range(0, nr2));
    Matrix m2new = m2;
    mat.swap(m2new);
  }
}

void trim_op(const DiagInfo &diag, CustomOp &allops) {
  for (auto &[name, op] : allops) 
    trim_matel(diag, op);
}

// Object of class IterInfo cotains full information about matrix representations when entering stage N of the NRG
// iteration.
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

   void dump_diagonal(const size_t max_nr, ostream &F = std::cout) const {
     if (max_nr) {
       for (const auto &[name, m] : ops)  dump_diagonal_op(name, m, max_nr, F);
       for (const auto &[name, m] : opsg) dump_diagonal_op(name, m, max_nr, F);
     }
   }
   void trim_matrices(const DiagInfo &diag) {
     trim_op(diag, ops);
     trim_op(diag, opsp);
     trim_op(diag, opsg);
     trim_op(diag, opd);
     trim_op(diag, opt);
     trim_op(diag, opot);
     trim_op(diag, opq);
   }
};

#include "spectral.h"

#include "coef.cc"
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
#include "sym-SPU1LR.cc"
#include "sym-SU2.cc"
#include "sym-QSLR.cc"
#include "sym-QST.cc"
#include "sym-QSTZ.cc"
#include "sym-QSZTZ.cc"
#include "sym-QSZLR.cc"
#include "sym-QJ.cc"
#include "sym-U1.cc"
 #ifdef NRG_COMPLEX
 #include "sym-SPSU2C3.cc"
 #include "sym-QSC3.cc"
 #endif
#endif

// Operator sumrules.
template<typename F> double norm(const MatrixElements &m, shared_ptr<Symmetry> Sym, F factor_fnc, int SPIN) {
  weight_bucket sum;
  for (const auto &[II, mat] : m) {
    const auto & [I1, Ip] = II;
    if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
    sum += factor_fnc(Ip, I1) * frobenius_norm(mat);
  }
  return 2.0 * cmpl(sum).real(); // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
}

void operator_sumrules(const IterInfo &a, shared_ptr<Symmetry> Sym) {
  // We check sum rules wrt some given spin (+1/2, by convention). For non-spin-polarized calculations, this is
  // irrelevant (0).
  const int SPIN = Sym->isfield() ? 1 : 0;
  for (const auto &[name, m] : a.opd)
    cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensFactorFnc(), SPIN) << std::endl;
  for (const auto &[name, m] : a.opq)
    cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensquadFactorFnc(), 0) << std::endl;
}

#include "read-input.cc"

template <typename T>
  std::string formatted_output(T x, const Params &P) {
    return fmt::format("{x:>{width}}", "x"_a=x, "width"_a=P.width_custom);
  }

std::string formatted_output(double x, const Params &P) {
  return fmt::format("{x:>{width}.{prec}}", "x"_a=x, "prec"_a=P.prec_custom, "width"_a=P.width_custom);
}

template <typename T>
  bool negligible_imag_part(std::complex<T> z, const double output_imag_eps = 1e-13) {
    return abs(z.imag()) < abs(z.real()) * output_imag_eps; 
  }

// The output format for complex values is X+IY or X-IY, where X and Y are real and imaginary part, respectively. The
// imaginary part is only shown where its value relative to the real part is sufficiently large. No space is used in
// the outputted string in order to simplify parsing.

std::string formatted_output(cmpl z, const Params &P) {
  const auto [r, i] = reim(z);
  const auto str = P.noimag || negligible_imag_part(z) ?
    fmt::format("{r:.{prec}f}", "r"_a=r, "prec"_a=P.prec_custom) :
    fmt::format("{r:.{prec}f}{s}I{absi:.{prec}f}", "r"_a=r, "s"_a=(i>0 ? "+" : "-"), "absi"_a=abs(i), "prec"_a=P.prec_custom);
  return fmt::format("{str:>{width}}", "str"_a=str, "width"_a=P.width_custom); // the width for the whole X+iY string
}

#include "bins.h"

// Object of class 'ChainSpectrum' will contain information about the spectral density calculated at a given stage of
// the NRG run, i.e. for a finite Wilson chain.  We then merge it in an object of class 'Spectrum' which holds the
// spectral information for the entire run (i.e. the physical spectral density).

class ChainSpectrum {
 private:
   const Params &P;
public:
   explicit ChainSpectrum(const Params &P) : P(P) {}
   virtual void add(double energy, t_weight weight) = 0;
   virtual ~ChainSpectrum() = default; // required, because there are virtual members
};

class ChainSpectrumBinning : public ChainSpectrum {
 private:
   Bins spos, sneg;
 public:
   explicit ChainSpectrumBinning(const Params &P) : ChainSpectrum(P), spos(P), sneg(P) {}
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
   explicit ChainSpectrumTemp(const Params &P) : ChainSpectrum(P), v(P) {}
   void add(double T, t_weight value) override { v.add_value(T, value); }
   friend class SpectrumTemp;
};

#include "matsubara.h"

class ChainSpectrumMatsubara : public ChainSpectrum {
 private:
   Matsubara m;
 public:
   ChainSpectrumMatsubara(const Params &P, matstype mt) : ChainSpectrum(P), m(P.mats, mt, P.T){};
   void add(size_t n, t_weight w) { m.add(n, w); }
   void add(double energy, t_weight w) override { my_assert_not_reached(); }
   t_weight total_weight() const { return m.total_weight(); }
   friend class SpectrumMatsubara;
};

// Object of class spectrum will contain everything that we know about a spectral density.
class Spectrum {
 public:
   string opname, filename;
   shared_ptr<Algo> algotype;
   const Params &P;
   Spectrum(const string &opname, const string &filename, shared_ptr<Algo> algotype, const Params &P) :
     opname(opname), filename(filename), algotype(algotype), P(P) {}; // NOLINT
   virtual ~Spectrum()= default; // required (the destructor saves the results to a file)
   virtual void merge(std::shared_ptr<ChainSpectrum>, const Step &) = 0; // called from spec.cc as the very last step
   string name() { return opname; }
};

#include "spectrumrealfreq.cc"

// G(T) type of results, i.e. not a real spectrum
class SpectrumTemp : public Spectrum {
 private:
   Spikes results;
 public:
   SpectrumTemp(const string &opname, const string &filename, shared_ptr<Algo> algotype, const Params &P) : 
     Spectrum(opname, filename, algotype, P) {}
   void merge(std::shared_ptr<ChainSpectrum> cs, const Step &) override {
     auto t = dynamic_pointer_cast<ChainSpectrumTemp>(cs);
     std::copy(begin(t->v), end(t->v), std::back_inserter(results));
   }
   ~SpectrumTemp() override {
     cout << "Spectrum: " << opname << " " << algotype->name() << endl;
     ranges::sort(results, sortfirst());
     results.save(safe_open(filename + ".dat"), P.prec_xy, P.reim);
   }
};

// This container actually holds the GF on the Matsubara axis, not a spectral function.
class SpectrumMatsubara : public Spectrum {
 private:
   Matsubara results;
 public:
   SpectrumMatsubara(const string &opname, const string &filename, shared_ptr<Algo> algotype, matstype mt, const Params &P)
     : Spectrum(opname, filename, algotype, P), results(P.mats, mt, P.T) {}
   void merge(std::shared_ptr<ChainSpectrum> cs, const Step &) override {
     auto t = dynamic_pointer_cast<ChainSpectrumMatsubara>(cs);
     for (const auto n : range0(results.v.size())) results.v[n].second += t->m.v[n].second;
   }     
   ~SpectrumMatsubara() override { 
     fmt::print(fmt::emphasis::bold, "Spectrum: {} -> {}\n", opname, algotype->name());
     results.save(safe_open(filename + ".dat"), P.prec_xy);
   }
};

// Check if the trace of the density matrix equals 'ref_value'.
void check_trace_rho(const DensMatElements &m, shared_ptr<Symmetry> Sym, double ref_value = 1.0) {
  if (!num_equal(m.trace(Sym->multfnc()), ref_value))
    throw std::runtime_error("check_trace_rho() failed");
}

enum class axis { RealFreq, Temp, Matsubara };

std::ostream & operator<<(std::ostream &os, const axis &a) {
  if (a == axis::RealFreq)  os << "RealFreq";
  if (a == axis::Temp)      os << "Temp";
  if (a == axis::Matsubara) os << "Matsubara";
  return os;
}

inline std::string to_string(std::complex<double> &z) {
  ostringstream s;
  s << z;
  return s.str();
}

// All information about calculating a spectral function: pointers to the operator data, raw spectral data acccumulators,
// algorithm, etc.
class BaseSpectrum {
 public:
   string name;
   string prefix; // "dens", "corr", etc.
   const MatrixElements &op1, &op2;
   shared_ptr<Spectrum> spec;  
   shared_ptr<Algo>     algotype; // Algo_FDM, Algo_DMNRG,...
   axis a;              // axis::RealFreq, axis::Temp, axis::Matsubara, etc.
   matstype mt;         // matstype::bosonic, matstype::fermionic, etc.
   int spin{};          // -1 or +1, or 0 where irrelevant
   std::string fullname() const {
     std::string s = fmt::format("{} {} {} {}", name, prefix, algotype ? algotype->name() : "[Algorithm not set]", a);
     if (a == axis::Matsubara) s += " " + matstypestring(mt);
     return s;
   }
   BaseSpectrum(const MatrixElements &op1, const MatrixElements &op2, const string name, const string prefix, const matstype mt, const int spin) :
     name(name), prefix(prefix), op1(op1), op2(op2), a(axis::RealFreq), mt(mt), spin(spin) {}
};
using speclist = std::list<BaseSpectrum>;

#include "spec.cc"
#include "dmnrg.h"
#include "splitting.cc"

// Determine the ranges of index r
Rmaxvals::Rmaxvals(const Invar &I, const InvarVec &InVec, const DiagInfo &diagprev, shared_ptr<Symmetry> Sym) {
  for (const auto &[i, In] : InVec | ranges::views::enumerate)
    values.push_back(Sym->triangle_inequality(I, In, Sym->QN_subspace(i)) ? diagprev.size_subspace(In) : 0);
}

// Formatted output of the computed expectation values
class ExpvOutput {
 private:
   ofstream F;                // output stream
   map<string, t_expv> &m;    // reference to the name->value mapping
   const list<string> fields; // list of fields to be output (may be a subset of the fields actually present in m)
   const Params &P;
   void field_numbers() {     // Consecutive numbers for the columns
     F << '#' << formatted_output(1, P) << ' ';
     for (const auto ctr : range1(fields.size())) F << formatted_output(1 + ctr, P) << ' ';
     F << endl;
   }
   // Label and field names. Label is the first column (typically the temperature).
   void field_names(string labelname = "T") {
     F << '#' << formatted_output(labelname, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(op, P); });
     F << endl;
   }
 public:
   // Output the current values for the label and for all the fields
   void field_values(double labelvalue, bool cout_dump = true) {
     F << ' ' << formatted_output(labelvalue, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(m[op], P); });
     F << endl;
     if (cout_dump)
       for (const auto &op: fields)
         fmt::print(fmt::emphasis::bold | fg(fmt::color::red), "<{}>={}\n", op, to_string(m[op]));
   }
   ExpvOutput(const string &fn, map<string, t_expv> &m_, const list<string> &fields_, const Params &P_) : m(m_), fields(fields_), P(P_) {
     F.open(fn);
     field_numbers();
     field_names();
   }
};

// open_files() and open_files_spec() open the output files and establishe the data structures for storing spectral
// information.
void open_files(speclist &sl, BaseSpectrum &spec, shared_ptr<Algo> algotype, const axis a, const Params &P) {
  const string fn = spec.prefix + "_" + (algotype ? algotype->name() : "OOPS") + "_dens_" + spec.name; // no suffix (.dat vs. .bin)
  switch (a) {
    case axis::RealFreq:  spec.spec = make_shared<SpectrumRealFreq>(spec.name, fn, algotype, P); break;
    case axis::Temp:      spec.spec = make_shared<SpectrumTemp>(spec.name, fn, algotype, P); break;
    case axis::Matsubara: spec.spec = make_shared<SpectrumMatsubara>(spec.name, fn, algotype, spec.mt, P); break;
    default: my_assert_not_reached();
  }
  spec.algotype = algotype;
  spec.a        = a;
  sl.push_back(spec);
}

void open_files_spec(const RUNTYPE &runtype, speclist &sl, BaseSpectrum &spec, const Params &P) {
  if (spec.prefix == "gt") {
    if (runtype == RUNTYPE::NRG) open_files(sl, spec, make_shared<Algo_GT>(P), axis::Temp, P);
    return;
  }
  if (spec.prefix == "i1t") {
    if (runtype == RUNTYPE::NRG) open_files(sl, spec, make_shared<Algo_I1T>(P), axis::Temp, P);
    return;
  }
  if (spec.prefix == "i2t") {
    if (runtype == RUNTYPE::NRG) open_files(sl, spec, make_shared<Algo_I2T>(P), axis::Temp, P);
    return;
  }
  if (spec.prefix == "chit") {
    if (runtype == RUNTYPE::NRG) open_files(sl, spec, make_shared<Algo_CHIT>(P), axis::Temp, P);
    return;
  }
  // If we did not return from this funciton by this point, what we are computing is the spectral function. There are
  // several possibilities in this case, all of which may be enabled at the same time.
  if (runtype == RUNTYPE::NRG) {
    if (P.finite) open_files(sl, spec, make_shared<Algo_FT>(P), axis::RealFreq, P);
    if (P.finitemats) open_files(sl, spec, make_shared<Algo_FTmats>(P), axis::Matsubara, P);
  }
  if (runtype == RUNTYPE::DMNRG) {
    if (P.dmnrg) open_files(sl, spec, make_shared<Algo_DMNRG>(P), axis::RealFreq, P);
    if (P.dmnrgmats) open_files(sl, spec, make_shared<Algo_DMNRGmats>(P), axis::Matsubara, P);
    if (P.cfs) open_files(sl, spec, make_shared<Algo_CFS>(P), axis::RealFreq, P);
    if (P.cfsgt) open_files(sl, spec, make_shared<Algo_CFSgt>(P), axis::RealFreq, P);
    if (P.cfsls) open_files(sl, spec, make_shared<Algo_CFSls>(P), axis::RealFreq, P);
    if (P.fdm) open_files(sl, spec, make_shared<Algo_FDM>(P), axis::RealFreq, P);
    if (P.fdmgt) open_files(sl, spec, make_shared<Algo_FDMgt>(P), axis::RealFreq, P);
    if (P.fdmls) open_files(sl, spec, make_shared<Algo_FDMls>(P), axis::RealFreq, P);
    if (P.fdmmats) open_files(sl, spec, make_shared<Algo_FDMmats>(P), axis::Matsubara, P);
  }
}

template <typename T> ostream & operator<<(ostream &os, const std::set<T> &x) {
  std::copy(x.cbegin(), x.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}

class Oprecalc {
 public:
   // The following lists hold the names of operators which need to be recomputed. The default behavior is to
   // recompute all the operators that are required to calculate the requested spectral densities, see function
   // open_files(). In addition, singlet operators are always recomputed in the first NRG run, so that we can
   // calculate the expectation values.
   set<string> s, p, g, d, v, t, q, ot;

   speclist spectraD, spectraS, spectraT, spectraQ, spectraGT, spectraI1T, spectraI2T, spectraK, spectraCHIT, spectraC, spectraOT;

   // Calculate spectral densities
   void spectral_densities(const Step &step, const DiagInfo &diag, DensMatElements &rho, DensMatElements &rhoFDM, const Stats &stats, shared_ptr<Symmetry> Sym) {
     TIME("spec");
     for (auto &i : spectraS)    calc_generic(i, step, diag, Sym->CorrelatorFactorFnc(),   Sym->TrivialCheckSpinFnc(),  rho, rhoFDM, stats);
     for (auto &i : spectraCHIT) calc_generic(i, step, diag, Sym->CorrelatorFactorFnc(),   Sym->TrivialCheckSpinFnc(),  rho, rhoFDM, stats);
     for (auto &i : spectraD)    calc_generic(i, step, diag, Sym->SpecdensFactorFnc(),     Sym->SpecdensCheckSpinFnc(), rho, rhoFDM, stats);
     for (auto &i : spectraT)    calc_generic(i, step, diag, Sym->SpinSuscFactorFnc(),     Sym->TrivialCheckSpinFnc(),  rho, rhoFDM, stats);
     for (auto &i : spectraOT)   calc_generic(i, step, diag, Sym->OrbSuscFactorFnc(),      Sym->TrivialCheckSpinFnc(),  rho, rhoFDM, stats);
     for (auto &i : spectraQ)    calc_generic(i, step, diag, Sym->SpecdensquadFactorFnc(), Sym->TrivialCheckSpinFnc(),  rho, rhoFDM, stats);
     for (auto &i : spectraGT)   calc_generic(i, step, diag, Sym->SpecdensFactorFnc(),     Sym->SpecdensCheckSpinFnc(), rho, rhoFDM, stats);
     for (auto &i : spectraI1T)  calc_generic(i, step, diag, Sym->SpecdensFactorFnc(),     Sym->SpecdensCheckSpinFnc(), rho, rhoFDM, stats);
     for (auto &i : spectraI2T)  calc_generic(i, step, diag, Sym->SpecdensFactorFnc(),     Sym->SpecdensCheckSpinFnc(), rho, rhoFDM, stats);
   }

   void report(ostream &F, const string &name, const set<string> &x) {
     F << name << "=[" << x << "]" << endl;
   }

   void report(ostream &F = cout) {
     F << endl << "Computing the following operators:" << endl;
     report(F, "s", s);
     report(F, "p", p);
     report(F, "g", g);
     report(F, "d", d);
     report(F, "v", v);
     report(F, "t", t);
     report(F, "q", q);
     report(F, "ot", ot);
   }
   
   bool do_s(const string &name, const Params &P, const Step &step) {
     if (step.nrg()) return true;                                          // for computing <O> 
     if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
     return s.count(name);
   }
   
   bool do_g(const string &name, const Params &P, const Step &step) {
     if (step.nrg()) return true;                                          // for computing <O>
     if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
     return g.count(name);
   }
   
   // Wrapper routine for recalculations
   template <typename RecalcFnc>
     MatrixElements recalc_common(const MatrixElements &mold, RecalcFnc recalc_fnc, const Step &step, const DiagInfo &diag,
                                  const QSrmax &qsrmax, const std::string name, const string &tip, shared_ptr<Symmetry> Sym, const Params &P) {
       nrglog('0', "Recalculate " << tip << " " << name);
       auto mnew = recalc_fnc(diag, qsrmax, mold);
       if (tip == "g") Sym->recalc_global(step, diag, qsrmax, name, mnew);
       return mnew;
     }
   
   template <typename ... Args>
     MatrixElements recalc_or_clear(bool recalc, Args&& ... args) {
       return recalc ? recalc_common(std::forward<Args>(args)...) : MatrixElements();
     }

   // Recalculate operator matrix representations
   ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
     void recalculate_operators(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, IterInfo &a, shared_ptr<Symmetry> Sym, const Params &P) {
       for (auto &[name, m] : a.ops)
         m = recalc_or_clear(do_s(name, P, step), m, [Sym](const auto &... pr) { return Sym->recalc_singlet(pr..., 1);  }, step, diag, qsrmax, name, "s", Sym, P);
       for (auto &[name, m] : a.opsp)
         m = recalc_or_clear(p.count(name),       m, [Sym](const auto &... pr) { return Sym->recalc_singlet(pr..., -1); }, step, diag, qsrmax, name, "p", Sym, P);
       for (auto &[name, m] : a.opsg) 
         m = recalc_or_clear(do_g(name, P, step), m, [Sym](const auto &... pr) { return Sym->recalc_singlet(pr...,  1); }, step, diag, qsrmax, name, "g", Sym, P);
       for (auto &[name, m] : a.opd)
         m = recalc_or_clear(d.count(name),       m, [Sym](const auto &... pr) { return Sym->recalc_doublet(pr...);     }, step, diag, qsrmax, name, "d", Sym, P);
       for (auto &[name, m] : a.opt)
         m = recalc_or_clear(t.count(name),       m, [Sym](const auto &... pr) { return Sym->recalc_triplet(pr...);     }, step, diag, qsrmax, name, "t", Sym, P);
       for (auto &[name, m] : a.opot)
         m = recalc_or_clear(ot.count(name),      m, [Sym](const auto &... pr) { return Sym->recalc_orb_triplet(pr...); }, step, diag, qsrmax, name, "ot", Sym, P);
       for (auto &[name, m] : a.opq)
         m = recalc_or_clear(q.count(name),       m, [Sym](const auto &... pr) { return Sym->recalc_quadruplet(pr...);  }, step, diag, qsrmax, name, "q", Sym, P);
     }

   // Construct the suffix of the filename for spectral density files: 'A_?-A_?'.
   // If SPIN == 1 or SPIN == -1, '-u' or '-d' is appended to the string.
   string sdname(const string &a, const string &b, int spin = 0) {
     return a + "-" + b + (spin == 0 ? "" : (spin == 1 ? "-u" : "-d"));
   }

   void loopover(const RUNTYPE &runtype, const Params &P,
                 const CustomOp &set1, const CustomOp &set2,
                 const string_token &stringtoken, speclist &spectra, const string &prefix,
                 set<string> &rec1, set<string> &rec2, matstype mt, int spin = 0) {
    for (const auto &[name1, op1] : set1) {
      for (const auto &[name2, op2] : set2) {
        if (const auto name = sdname(name1, name2, spin); stringtoken.find(name)) {
          BaseSpectrum spec(op1, op2, name, prefix, mt, spin);
          open_files_spec(runtype, spectra, spec, P);
          std::cout << "Spectrum: " << spec.fullname() << std::endl;
          rec1.insert(name1);
          rec2.insert(name2);
        }
      }
    }
  }

  // Reset lists of operators which need to be iterated
  Oprecalc(const RUNTYPE &runtype, const IterInfo &a, shared_ptr<Symmetry> Sym, const Params &P) {
    cout << endl << "Computing the following spectra:" << endl;
    // Correlators (singlet operators of all kinds)
    string_token sts(P.specs);
    loopover(runtype, P, a.ops,  a.ops,  sts, spectraS, "corr", s, s, matstype::bosonic);
    loopover(runtype, P, a.opsp, a.opsp, sts, spectraS, "corr", p, p, matstype::bosonic);
    loopover(runtype, P, a.opsg, a.opsg, sts, spectraS, "corr", g, g, matstype::bosonic);
    loopover(runtype, P, a.ops,  a.opsg, sts, spectraS, "corr", s, g, matstype::bosonic);
    loopover(runtype, P, a.opsg, a.ops,  sts, spectraS, "corr", g, s, matstype::bosonic);
    // Global susceptibilities (global singlet operators)
    string_token stchit(P.specchit);
    loopover(runtype, P, a.ops,  a.ops,  stchit, spectraCHIT, "chit", s, s, matstype::bosonic);
    loopover(runtype, P, a.ops,  a.opsg, stchit, spectraCHIT, "chit", s, g, matstype::bosonic);
    loopover(runtype, P, a.opsg, a.ops,  stchit, spectraCHIT, "chit", g, s, matstype::bosonic);
    loopover(runtype, P, a.opsg, a.opsg, stchit, spectraCHIT, "chit", g, g, matstype::bosonic);
    // Dynamic spin susceptibilities (triplet operators)
    string_token stt(P.spect);
    loopover(runtype, P, a.opt, a.opt, stt, spectraT, "spin", t, t, matstype::bosonic);
    string_token stot(P.specot);
    loopover(runtype, P, a.opot, a.opot, stot, spectraOT, "orbspin", ot, ot, matstype::bosonic);
    const int varmin = (Sym->isfield() ? -1 : 0);
    const int varmax = (Sym->isfield() ? +1 : 0);
    // Spectral functions (doublet operators)
    string_token std(P.specd);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, std, spectraD, "spec", d, d, matstype::fermionic, SPIN);
    string_token stgt(P.specgt);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, stgt, spectraGT, "gt", d, d, matstype::fermionic, SPIN);
    string_token sti1t(P.speci1t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, sti1t, spectraI1T, "i1t", d, d, matstype::fermionic, SPIN);
    string_token sti2t(P.speci2t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, sti2t, spectraI2T, "i2t", d, d, matstype::fermionic, SPIN);
    // Spectral functions (quadruplet operators)
    string_token stq(P.specq);
    loopover(runtype, P, a.opq, a.opq, stq, spectraQ, "specq", q, q, matstype::fermionic);
    report();
  }
};

// Store eigenvalue & quantum numbers information (flow diagrams)
class Annotated {
 private:
   std::ofstream F;
   
   // scaled = true -> output scaled energies (i.e. do not multiply by the rescale factor)
   inline t_eigen scaled_energy(t_eigen e, const Step &step, const Stats &stats, bool scaled = true, bool absolute = false) {
     return e * (scaled ? 1.0 : step.scale()) + (absolute ? stats.total_energy : 0.0);
   }
   
   const Params &P;

 public:
   explicit Annotated(const Params &P) : P(P) {}
   
   void dump(const Step &step, const DiagInfo &diag, const Stats &stats, shared_ptr<Symmetry> Sym, const std::string filename = "annotated.dat") {
     if (!P.dumpannotated) return;
     if (!F.is_open()) { // open output file
       F.open(filename);
       F << setprecision(P.dumpprecision);
     }
     std::vector<pair<t_eigen, Invar>> seznam;
     for (const auto &[I, eig] : diag)
       for (const auto e : eig.value_zero)
         seznam.emplace_back(e, I);
     ranges::sort(seznam);
     size_t len = std::min<size_t>(seznam.size(), P.dumpannotated); // non-const
     // If states are clustered, we dump the full cluster
     while (len < seznam.size()-1 && my_fcmp(seznam[len].first, seznam[len-1].first, P.grouptol) == 0) len++;
     auto scale = [&step, &stats, this](auto x) { return scaled_energy(x, step, stats, P.dumpscaled, P.dumpabs); };
     if (P.dumpgroups) {
       // Group by degeneracies
       for (size_t i = 0; i < len;) { // i increased in the while loop below
         const auto [e0, I0] = seznam[i];
         F << scale(e0);
         std::vector<string> QNstrings;
         size_t total_degeneracy = 0; // Total number of levels (incl multiplicity)
         while (i < len && my_fcmp(seznam[i].first, e0, P.grouptol) == 0) {
           const auto [e, I] = seznam[i];
           QNstrings.push_back(to_string(I));
           total_degeneracy += Sym->mult(I);
           i++;
         }
         ranges::sort(QNstrings);
         for (const auto &j : QNstrings) F << " (" << j << ")";
         F << " [" << total_degeneracy << "]" << endl;
       }
     } else {
       seznam.resize(len); // truncate!
       for (const auto &[e, I] : seznam) 
         F << scale(e) << " " << I << endl;
     }
     F << std::endl; // Consecutive iterations are separated by an empty line
   }
};

// Handle all output
struct Output {
  RUNTYPE runtype;
  const Params &P;
  Annotated annotated;
  ofstream Fenergies;  // all energies (different file for NRG and for DMNRG)
  unique_ptr<ExpvOutput> custom;
  unique_ptr<ExpvOutput> customfdm;
  
  Output(const RUNTYPE &runtype, const IterInfo &iterinfo, Stats &stats, const Params &P,
         const std::string filename_energies= "energies.nrg"s,
         const std::string filename_custom = "custom", 
         const std::string filename_customfdm = "customfdm")
    : runtype(runtype), P(P), annotated(P) {
      // We dump all energies to separate files for NRG and DM-NRG runs. This is a very convenient way to check if both
      // runs produce the same results.
      if (P.dumpenergies && runtype == RUNTYPE::NRG) Fenergies.open(filename_energies);
      list<string> ops;
      for (const auto &name : iterinfo.ops  | boost::adaptors::map_keys) ops.push_back(name);
      for (const auto &name : iterinfo.opsg | boost::adaptors::map_keys) ops.push_back(name);
      if (runtype == RUNTYPE::NRG)
        custom = make_unique<ExpvOutput>(filename_custom, stats.expv, ops, P);
      else if (runtype == RUNTYPE::DMNRG && P.fdmexpv) 
        customfdm = make_unique<ExpvOutput>(filename_customfdm, stats.fdmexpv, ops, P);
    }

  // Dump all energies in diag to a file
  void dump_all_energies(const DiagInfo &diag, int N) {
    if (!Fenergies) return;
    Fenergies << endl << "===== Iteration number: " << N << endl;
    diag.dump_value_zero(Fenergies);
  }
};

CONSTFNC t_expv calc_trace_singlet(const Step &step, const DiagInfo &diag, const MatrixElements &n, shared_ptr<Symmetry> Sym) {
  matel_bucket tr; // note: t_matel = t_expv
  for (const auto &[I, eig] : diag) {
    const auto & nI = n.at({I,I});
    const auto dim = eig.getnrstored();
    my_assert(dim == nI.size2());
    matel_bucket sum;
    for (const auto r : range0(dim)) sum += exp(-step.TD_factor() * eig.value_zero(r)) * nI(r, r);
    tr += Sym->mult(I) * t_matel(sum);
  }
  return tr;
}

// Measure thermodynamic expectation values of singlet operators
void measure_singlet(const Step &step, Stats &stats, const DiagInfo &diag, const IterInfo &a, Output &output, shared_ptr<Symmetry> Sym, const Params &P) {
  const auto Z = ranges::accumulate(diag, t_expv{}, [&Sym, &step](auto total, const auto &d) { const auto &[I, eig] = d;
    return total + Sym->mult(I) * ranges::accumulate(eig.value_zero, t_expv{},
                                                     [f=step.TD_factor()](auto sum, const auto &x) { return sum + exp(-f*x); }); });
  for (const auto &[name, m] : a.ops)  stats.expv[name] = calc_trace_singlet(step, diag, m, Sym) / Z;
  for (const auto &[name, m] : a.opsg) stats.expv[name] = calc_trace_singlet(step, diag, m, Sym) / Z;
  output.custom->field_values(step.Teff());
}

template<typename T>
  T trace_contract(const ublas::matrix<T> &A, const ublas::matrix<T> &B, const size_t range)
{
  T sum{};
  for (const auto i : range0(range))
       for (const auto j : range0(range)) 
      sum += A(i, j) * B(j, i);
  return sum;
}

CONSTFNC t_expv calc_trace_fdm_kept(size_t ndx, const MatrixElements &n, const DensMatElements &rhoFDM, const AllSteps &dm, shared_ptr<Symmetry> Sym) {
  matel_bucket tr;
  for (const auto &[I, rhoI] : rhoFDM)
    tr += Sym->mult(I) * trace_contract(rhoI, n.at({I,I}), dm[ndx].at(I).kept); // over kept states ONLY
  return tr;
}

// Expectation values using the FDM algorithm
void measure_singlet_fdm(const Step &step, Stats &stats, const DiagInfo &diag, const IterInfo &a, Output &output, 
                         const DensMatElements &rhoFDM, const AllSteps &dm, shared_ptr<Symmetry> Sym, const Params &P) {
  for (const auto &[name, m] : a.ops)  stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, dm, Sym);
  for (const auto &[name, m] : a.opsg) stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, dm, Sym);
  output.customfdm->field_values(P.T);
}

// Calculate grand canonical partition function at current NRG energy shell. This is not the same as the true
// partition function of the full problem! Instead this is the Z_N that is used to initialize the density matrix,
// i.e. rho = 1/Z_N \sum_{l} exp{-beta E_l} |l;N> <l;N|.  calc_grand_canonical_Z() is also used to calculate
// stats.Zft, that is used to compute the spectral function with the conventional approach, as well as stats.Zgt for
// G(T) calculations, stats.Zchit for chi(T) calculations.
double calc_grand_canonical_Z(const Step &step, const DiagInfo &diag, shared_ptr<Symmetry> Sym, const double factor = 1.0) {
  bucket ZN;
  for (const auto &[I, eig]: diag) 
    for (const auto &i : eig.kept()) // sum over all kept states
      ZN += Sym->mult(I) * exp(-eig.value_zero(i) * step.scT() * factor);
  my_assert(ZN >= 1.0);
  return ZN;
}

Matrix diagonal_exp(const Eigen &eig, const double factor)
{
  const auto dim = eig.getnrstored();
  Matrix m(dim, dim, 0);
  for (const auto i: range0(dim)) 
      m(i, i) = exp(-eig.value_zero(i) * factor);
  return m;
}

// Calculate rho_N, the density matrix at the last NRG iteration. It is
// normalized to 1. Note: in CFS approach, we consider all states in the
// last iteration to be "discarded".
// For the details on the full Fock space approach see:
// F. B. Anders, A. Schiller, Phys. Rev. Lett. 95, 196801 (2005).
// F. B. Anders, A. Schiller, Phys. Rev. B 74, 245113 (2006).
// R. Peters, Th. Pruschke, F. B. Anders, Phys. Rev. B 74, 245114 (2006).
DensMatElements init_rho(const Step &step, const DiagInfo &diag, shared_ptr<Symmetry> Sym) {
  DensMatElements rho;
  for (const auto &[I, eig]: diag)
    rho[I] = diagonal_exp(eig, step.scT()) / calc_grand_canonical_Z(step, diag, Sym);
  check_trace_rho(rho, Sym);
  return rho;
}

// Determine the number of states to be retained. Returns Emax - the highest energy to still be retained.
t_eigen highest_retained_energy(const Step &step, const DiagInfo &diag, const Params &P) {
  auto energies = diag.sorted_energies();
  my_assert(energies.front() == 0.0); // check for the subtraction of Egs
  const size_t totalnumber = energies.size();
  size_t nrkeep;
  if (P.keepenergy <= 0.0) {
    nrkeep = P.keep;
  } else {
    double keepenergy = P.keepenergy * step.unscale();
    // We add 1 for historical reasons. We thus keep states with E<=Emax, and one additional state which has E>Emax.
    nrkeep = 1 + ranges::count_if(energies, [=](double e) { return e <= keepenergy; });
    nrkeep = std::clamp<size_t>(nrkeep, P.keepmin, P.keep);
  }
  // Check for near degeneracy and ensure that the truncation occurs in a "gap" between nearly-degenerate clusters of
  // eigenvalues.
  if (P.safeguard > 0.0) {
    size_t cnt_extra = 0;
    while (nrkeep < totalnumber && (energies[nrkeep] - energies[nrkeep - 1]) <= P.safeguard && cnt_extra < P.safeguardmax) {
      nrkeep++;
      cnt_extra++;
    }
    if (cnt_extra) debug("Safeguard: keep additional " << cnt_extra << " states");
  }
  nrkeep = std::clamp<size_t>(nrkeep, 1, totalnumber);
  return energies[nrkeep - 1];
}

struct truncate_stats {
  size_t nrall, nrallmult, nrkept, nrkeptmult;
  truncate_stats(const DiagInfo &diag, shared_ptr<Symmetry> Sym) {
    nrall = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getdim(); });
    nrallmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getdim(); });
    nrkept = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getnrkept(); });
    nrkeptmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getnrkept(); });
  }
  void report() {
    nrgdump4(nrkept, nrkeptmult, nrall, nrallmult) << std::endl;
  }
};

struct NotEnough : public std::exception {};

// Compute the number of states to keep in each subspace. Returns true if an insufficient number of states has been
// obtained in the diagonalization and we need to compute more states.
void truncate_prepare(const Step &step, DiagInfo &diag, shared_ptr<Symmetry> Sym, const Params &P) {
  const auto Emax = highest_retained_energy(step, diag, P);
  for (auto &[I, eig] : diag)
    diag[I].truncate_prepare_subspace(step.last() && P.keep_all_states_in_last_step() ? eig.getnrcomputed() :
                                      ranges::count_if(eig.value_zero, [Emax](const double e) { return e <= Emax; }));
  std::cout << "Emax=" << Emax/step.unscale() << " ";
  truncate_stats ts(diag, Sym);
  ts.report();
  if (ranges::any_of(diag, [Emax](const auto &d) { const auto &[I, eig] = d; 
    return eig.getnrkept() == eig.getnrcomputed() && eig.value_zero(eig.getnrcomputed()-1) != Emax && eig.getnrcomputed() < eig.getdim(); }))
      throw NotEnough();
  const double ratio = double(ts.nrkept) / ts.nrall;
  fmt::print(FMT_STRING("Kept: {} out of {}, ratio={:.3}\n"), ts.nrkept, ts.nrall, ratio);
}

// Calculate partial statistical sums, ZnD*, and the grand canonical Z (stats.ZZG), computed with respect to absolute
// energies. calc_ZnD() must be called before the second NRG run.
void calc_ZnD(const AllSteps &dm, Stats &stats, shared_ptr<Symmetry> Sym, const double T) {
  mpf_set_default_prec(400); // this is the number of bits, not decimal digits!
  for (const auto N : dm.Nall()) {
    my_mpf ZnDG, ZnDN; // arbitrary-precision accumulators to avoid precision loss
    mpf_set_d(ZnDG, 0.0);
    mpf_set_d(ZnDN, 0.0);
    for (const auto &[I, ds] : dm[N])
      for (const auto i : ds.all()) {
        my_mpf g, n;
        mpf_set_d(g, Sym->mult(I) * exp(-ds.eig.absenergyG[i]/T)); // absenergyG >= 0.0
        mpf_set_d(n, Sym->mult(I) * exp(-ds.eig.absenergyN[i]/T)); // absenergyN >= 0.0
        mpf_add(ZnDG, ZnDG, g);
        mpf_add(ZnDN, ZnDN, n);
      }
    mpf_set(stats.ZnDG[N], ZnDG);
    mpf_set(stats.ZnDN[N], ZnDN);
    stats.ZnDNd[N] = mpf_get_d(stats.ZnDN[N]);
  }
  // Note: for ZBW, Nlen=Nmax+1. For Ninit=Nmax=0, index 0 will thus be included here.
  my_mpf ZZG;
  mpf_set_d(ZZG, 0.0);
  for (const auto N : dm.Nall()) {
    my_mpf a;
    mpf_set(a, stats.ZnDG[N]);
    my_mpf b;
    mpf_set_d(b, Sym->nr_combs());
    mpf_pow_ui(b, b, dm.Nend - N - 1);
    my_mpf c;
    mpf_mul(c, a, b);
    mpf_add(ZZG, ZZG, c);
  }
  stats.ZZG = mpf_get_d(ZZG);
  cout << "ZZG=" << HIGHPREC(stats.ZZG) << endl;
  for (const auto N : dm.Nall()) {
    const double w  = pow(Sym->nr_combs(), int(dm.Nend - N - 1)) / stats.ZZG;
    stats.wnfactor[N] = w; // These ratios enter the terms for the spectral function.
    stats.wn[N] = w * mpf_get_d(stats.ZnDG[N]); // This is w_n defined after Eq. (8) in the WvD paper.
  }
  const auto sumwn = ranges::accumulate(stats.wn, 0.0);
  cout << "sumwn=" << sumwn << " sumwn-1=" << sumwn - 1.0 << endl;
  my_assert(num_equal(sumwn, 1.0));  // Check the sum-rule.
}

void report_ZnD(Stats &stats, const Params &P) {
  for (const auto N : P.Nall())
    cout << "ZG[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDG[N])) << endl;
  for (const auto N : P.Nall())
    cout << "ZN[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDN[N])) << endl;
  for (const auto N : P.Nall())
    cout << "w[" << N << "]=" << HIGHPREC(stats.wn[N]) << endl;
  for (const auto N : P.Nall())
    cout << "wfactor[" << N << "]=" << HIGHPREC(stats.wnfactor[N]) << endl;
}

// TO DO: use Boost.Multiprecision instead of low-level GMP calls
// https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html
void fdm_thermodynamics(const AllSteps &dm, Stats &stats, shared_ptr<Symmetry> Sym, const double T)
{
  stats.td_fdm.T = T;
  stats.Z_fdm = stats.ZZG*exp(-stats.GS_energy/T); // this is the true partition function
  stats.td_fdm.F = stats.F_fdm = -log(stats.ZZG)*T+stats.GS_energy; // F = -k_B*T*log(Z)
  // We use multiple precision arithmetics to ensure sufficient accuracy in the calculation of
  // the variance of energy and thus the heat capacity.
  my_mpf E, E2;
  mpf_set_d(E, 0.0);
  mpf_set_d(E2, 0.0);
  for (const auto N : dm.Nall())
    if (stats.wn[N] > 1e-16) 
      for (const auto &[I, ds] : dm[N]) 
        for (const auto i : ds.all()) {
          my_mpf weight;
          mpf_set_d(weight, stats.wn[N] * Sym->mult(I) * exp(-ds.eig.absenergyN[i]/T));
          mpf_div(weight, weight, stats.ZnDN[N]);
          my_mpf e;
          mpf_set_d(e, ds.eig.absenergy[i]);
          my_mpf e2;
          mpf_mul(e2, e, e);
          mpf_mul(e, e, weight);
          mpf_mul(e2, e2, weight);
          mpf_add(E, E, e);
          mpf_add(E2, E2, e2);
        }
  stats.td_fdm.E = stats.E_fdm = mpf_get_d(E);
  my_mpf sqrE;
  mpf_mul(sqrE, E, E);
  my_mpf varE;
  mpf_sub(varE, E2, sqrE);
  stats.td_fdm.C = stats.C_fdm = mpf_get_d(varE)/pow(T,2);
  stats.td_fdm.S = stats.S_fdm = (stats.E_fdm-stats.F_fdm)/T;
  cout << endl;
  cout << "Z_fdm=" << HIGHPREC(stats.Z_fdm) << endl;
  cout << "F_fdm=" << HIGHPREC(stats.F_fdm) << endl;
  cout << "E_fdm=" << HIGHPREC(stats.E_fdm) << endl;
  cout << "C_fdm=" << HIGHPREC(stats.C_fdm) << endl;
  cout << "S_fdm=" << HIGHPREC(stats.S_fdm) << endl;
  cout << endl;
  stats.td_fdm.save_values();
}

template<typename F>
  double trace(F fnc, const double rescale_factor, const DiagInfo &diag, shared_ptr<Symmetry> Sym) {
    bucket b;
    for (const auto &[I, eig] : diag)
      b += Sym->mult(I) * ranges::accumulate(eig.value_zero, 0.0, [fnc, rescale_factor](double acc, const auto x) { 
        const double betaE = rescale_factor * x; return acc + fnc(betaE) * exp(-betaE); });
    return b;
  }

// We calculate thermodynamic quantities before truncation to make better use of the available states. Here we
// compute quantities which are defined for all symmetry types. Other calculations are performed by calculate_TD
// member functions defined in symmetry.cc.
void calculate_TD(const Step &step, const DiagInfo &diag, Stats &stats, Output &output, shared_ptr<Symmetry> Sym, double additional_factor = 1.0) {
  // Rescale factor for energies. The energies are expressed in units of omega_N, thus we need to appropriately
  // rescale them to calculate the Boltzmann weights at the temperature scale Teff (Teff=scale/betabar).
  double rescale_factor = step.TD_factor() * additional_factor;
  const auto Z = trace([](double x) { return 1; }, rescale_factor, diag, Sym); // partition function
  const auto E = trace([](double x) { return x; }, rescale_factor, diag, Sym); // Tr[beta H]
  const auto E2 = trace([](double x) { return sqr(x); }, rescale_factor, diag, Sym); // Tr[(beta H)^2]
  stats.Z = Z;
  stats.td.T  = step.Teff();
  stats.td.E  = E/Z;             // beta <H>
  stats.td.E2 = E2/Z;            // beta^2 <H^2>
  stats.td.C  = E2/Z - sqr(E/Z); // C/k_B=beta^2(<H^2>-<H>^2)
  stats.td.F  = -log(Z);         // F/(k_B T)=-ln(Z)
  stats.td.S  = E/Z+log(Z);      // S/k_B=beta<H>+ln(Z)
  Sym->calculate_TD(step, diag, stats, rescale_factor);  // symmetry-specific calculation routine
  stats.td.save_values();
}

void calculate_spectral_and_expv(const Step &step, Stats &stats, Output &output, Oprecalc &oprecalc, const DiagInfo &diag, 
                                 const IterInfo &iterinfo, const AllSteps &dm, shared_ptr<Symmetry> Sym, const Params &P) {
  // Zft is used in the spectral function calculations using the conventional approach. We calculate it here, in
  // order to avoid recalculations later on.
  stats.Zft = calc_grand_canonical_Z(step, diag, Sym);
  if (string(P.specgt) != "" || string(P.speci1t) != "" || string(P.speci2t) != "")
    stats.Zgt = calc_grand_canonical_Z(step, diag, Sym, 1.0/(P.gtp*step.scT()) ); // exp(-x*gtp)
  if (string(P.specchit) != "") 
    stats.Zchit = calc_grand_canonical_Z(step, diag, Sym, 1.0/(P.chitp*step.scT()) ); // exp(-x*chitp)
  DensMatElements rho, rhoFDM;
  if (step.dmnrg()) {
    if (P.need_rho()) {
      rho.load(step.ndx(), FN_RHO, P.removefiles);
      check_trace_rho(rho, Sym); // Check if Tr[rho]=1, i.e. the normalization
    }
    if (P.need_rhoFDM()) 
      rhoFDM.load(step.ndx(), FN_RHOFDM, P.removefiles);
  }
  oprecalc.spectral_densities(step, diag, rho, rhoFDM, stats, Sym);
  if (step.nrg()) {
    measure_singlet(step, stats, diag, iterinfo, output, Sym, P);
    iterinfo.dump_diagonal(P.dumpdiagonal);
  }
  if (step.dmnrg() && P.fdmexpv && step.N() == P.fdmexpvn) measure_singlet_fdm(step, stats, diag, iterinfo, output, rhoFDM, dm, Sym, P);
}

// Perform calculations of physical quantities. Called prior to NRG iteration (if calc0=true) and after each NRG
// step.
void perform_basic_measurements(const Step &step, const DiagInfo &diag, shared_ptr<Symmetry> Sym, Stats &stats, Output &output) {
  output.dump_all_energies(diag, step.ndx());
  calculate_TD(step, diag, stats, output, Sym);
  output.annotated.dump(step, diag, stats, Sym);
}

// Subspaces for the new iteration
auto new_subspaces(const DiagInfo &diagprev, shared_ptr<Symmetry> Sym) {
  std::set<Invar> subspaces;
  for (const auto &I : diagprev.subspaces()) {
    const auto all = Sym->new_subspaces(I);
    const auto non_empty = all | ranges::views::filter([&Sym](const auto &In) { return Sym->Invar_allowed(In); }) | ranges::to<std::vector>();
    std::copy(non_empty.begin(), non_empty.end(), std::inserter(subspaces, subspaces.end()));
  }
  return subspaces;
}

Matrix prepare_task_for_diag(const Step &step, const Invar &I, const Opch &opch, const Coef &coef, 
                             const DiagInfo &diagprev, shared_ptr<Symmetry> Sym, const Params &P) {
  const auto anc = Sym->ancestors(I);
  const Rmaxvals rm{I, anc, diagprev, Sym};
  Matrix h(rm.total(), rm.total(), 0);   // H_{N+1}=\lambda^{1/2} H_N+\xi_N (hopping terms)
  for (const auto i : Sym->combs())
    for (const auto r : range0(rm.rmax(i)))
      h(rm.offset(i) + r, rm.offset(i) + r) = P.nrg_step_scale_factor() * diagprev.at(anc[i]).value_zero(r);
  Sym->make_matrix(h, step, rm, I, anc, opch, coef);  // Symmetry-type-specific matrix initialization steps
  if (P.logletter('m')) dump_matrix(h);
  return h;
}

DiagInfo diagonalisations_OpenMP(const Step &step, const Opch &opch, const Coef &coef, const DiagInfo &diagprev, 
                                 const std::vector<Invar> &tasks, const DiagParams &DP, shared_ptr<Symmetry> Sym, const Params &P) {
  DiagInfo diagnew;
  size_t nr = tasks.size();
  size_t itask = 0;
  // cppcheck-suppress unreadVariable symbolName=nth
  int nth = P.diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
  for (itask = 0; itask < nr; itask++) {
    const Invar I  = tasks[itask];
    auto h = prepare_task_for_diag(step, I, opch, coef, diagprev, Sym, P);
    int thid = omp_get_thread_num();
#pragma omp critical
    { nrglog('(', "Diagonalizing " << I << " size=" << h.size1() << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
    Eigen e = diagonalise(h, DP);
#pragma omp critical
    { diagnew[I] = e; }
  }
  return diagnew;
}

#ifdef NRG_MPI
enum TAG : int { TAG_EXIT = 1, TAG_DIAG, TAG_SYNC, TAG_MATRIX, TAG_INVAR, TAG_MATRIX_SIZE, TAG_MATRIX_LINE, TAG_EIGEN_INT, TAG_EIGEN_VEC };

void mpi_send_params(const DiagParams &DP) {
  mpilog("Sending diag parameters " << DP.diag << " " << DP.diagratio);
  for (size_t i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_SYNC, 0);
  auto DPcopy = DP;
  mpi::broadcast(*mpiw, DPcopy, 0);
}

DiagParams mpi_receive_params() {
  DiagParams DP;
  mpi::broadcast(*mpiw, DP, 0);
  mpilog("Received diag parameters " << DP.diag << " " << DP.diagratio);
  return DP;
}

void check_status(mpi::status status) {
  if (status.error()) {
    cout << "MPI communication error. rank=" << mpiw->rank() << endl;
    mpienv->abort(1);
  }
}

// NOTE: MPI is limited to message size of 2GB (or 4GB). For big problems we thus need to send objects line by line.

void mpi_send_matrix(int dest, const Matrix &m) {
  auto size1 = m.size1();
  mpiw->send(dest, TAG_MATRIX_SIZE, size1);
  auto size2 = m.size2();
  mpiw->send(dest, TAG_MATRIX_SIZE, size2);
  mpilog("Sending matrix of size " << size1 << " x " << size2 << " line by line to " << dest);
  for (const auto i: range0(size1)) {
    ublas::vector<t_matel> vec = ublas::matrix_row<const Matrix>(m, i);
    mpiw->send(dest, TAG_MATRIX_LINE, vec);
  }
}

auto mpi_receive_matrix(int source) {
  size_t size1;
  check_status(mpiw->recv(source, TAG_MATRIX_SIZE, size1));
  size_t size2;
  check_status(mpiw->recv(source, TAG_MATRIX_SIZE, size2));
  Matrix m(size1, size2);
  mpilog("Receiving matrix of size " << size1 << " x " << size2 << " line by line from " << source);
  for (const auto i: range0(size1)) {
    ublas::vector<t_matel> vec;
    check_status(mpiw->recv(source, TAG_MATRIX_LINE, vec));
    my_assert(vec.size() == size2);
    ublas::matrix_row<Matrix>(m, i) = vec;
  }
  return m;
}

void mpi_send_eigen(int dest, const Eigen &eig) {
  mpilog("Sending eigen from " << mpiw->rank() << " to " << dest);
  mpiw->send(dest, TAG_EIGEN_VEC, eig.value_orig);
  mpi_send_matrix(dest, eig.matrix);
}

auto mpi_receive_eigen(int source) {
  mpilog("Receiving eigen from " << source << " on " << mpiw->rank());
  Eigen eig;
  check_status(mpiw->recv(source, TAG_EIGEN_VEC, eig.value_orig));
  eig.matrix = mpi_receive_matrix(source);
  return eig;
}

// Read results from a slave process.
std::pair<Invar, Eigen> read_from(int source) {
  mpilog("Reading results from " << source);
  auto eig = mpi_receive_eigen(source);
  Invar Irecv;
  check_status(mpiw->recv(source, TAG_INVAR, Irecv));
  mpilog("Received results for subspace " << Irecv << " [nr=" << eig.getnrstored() << ", dim=" << eig.getdim() << "]");
  my_assert(eig.value_orig.size() == eig.matrix.size1());
  my_assert(eig.matrix.size1() <= eig.matrix.size2());
  return {Irecv, eig};
}

DiagInfo diagonalisations_MPI(const Step &step, const Opch &opch, const Coef &coef, const DiagInfo &diagprev, 
                              const std::vector<Invar> &tasks, const DiagParams &DP, shared_ptr<Symmetry> Sym, const Params &P) {
  DiagInfo diagnew;
  mpi_send_params(DP); // Synchronise parameters
  list<Invar> todo; // List of all the tasks to handle
  copy(begin(tasks), end(tasks), back_inserter(todo));
  list<Invar> done; // List of finished tasks.
  // List of the available computation nodes (including the master,
  // which is always at the very beginnig of the deque).
  deque<int> nodes;
  for (const auto i: range0(mpiw->size())) nodes.push_back(i);
  nrglog('M', "nrtasks=" << tasks.size() << " nrnodes=" << mpiw->size());
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
    auto h = prepare_task_for_diag(step, I, opch, coef, diagprev, Sym, P);
    nrglog('M', "Scheduler: job " << I << " (dim=" << h.size1() << ")" << " on node " << i);
    if (i == 0) {
      // On master, diagonalize immediately.
      diagnew[I] = diagonalise(h, DP);
      nodes.push_back(0);
      done.push_back(I);
    } else {
      mpiw->send(i, TAG_DIAG, 0);
      mpi_send_matrix(i, h);
      mpiw->send(i, TAG_INVAR, I);
    }
    // Check for terminated jobs
    while (auto status = mpiw->iprobe(mpi::any_source, TAG_EIGEN_VEC)) {
      nrglog('M', "Receiveing results from " << status->source());
      auto [Irecv, eig] = read_from(status->source());
      diagnew[Irecv] = eig;
      done.push_back(Irecv);
      // The node is now available for new tasks!
      nodes.push_back(status->source());
    }
  }
  // Keep reading results sent from the slave processes until all tasks have been completed.
  while (done.size() != tasks.size()) {
    auto status = mpiw->probe(mpi::any_source, TAG_EIGEN_VEC);
    auto [Irecv, eig]  = read_from(status.source());
    diagnew[Irecv] = eig;
    done.push_back(Irecv);
  }
  return diagnew;
}
#endif

// Build matrix H(ri;r'i') in each subspace and diagonalize it
DiagInfo diagonalisations(const Step &step, const Opch &opch, const Coef &coef, const DiagInfo &diagprev, 
                          const std::vector<Invar> &tasks, double diagratio, shared_ptr<Symmetry> Sym, const Params &P) {
  TIME("diag");
#ifdef NRG_MPI
  return diagonalisations_MPI(step, opch, coef, diagprev, tasks, DiagParams(P, diagratio), Sym, P);
#else
  return diagonalisations_OpenMP(step, opch, coef, diagprev, tasks, DiagParams(P, diagratio), Sym, P);
#endif
}

// Determine the structure of matrices in the new NRG shell
QSrmax::QSrmax(const DiagInfo &diagprev, shared_ptr<Symmetry> Sym) {
  for (const auto &I : new_subspaces(diagprev, Sym))
    (*this)[I] = Rmaxvals{I, Sym->ancestors(I), diagprev, Sym};
}

// Recalculate irreducible matrix elements for Wilson chains.
void recalc_irreducible(const Step &step, const DiagInfo &diag, const QSrmax &qsrmax, Opch &opch, shared_ptr<Symmetry> Sym, const Params &P) {
  TIME("recalc f");
  if (!P.substeps) {
    opch = Sym->recalc_irreduc(step, diag, qsrmax);
  } else {
    const auto [N, M] = step.NM();
    for (const auto i: range0(size_t(P.channels)))
      if (i == M) {
        opch[i] = Sym->recalc_irreduc_substeps(step, diag, qsrmax, i);
      } else {
        for (const auto j: range0(size_t(P.perchannel)))
          opch[i][j] = Sym->recalc_doublet(diag, qsrmax, opch[i][j]);
      }
  }
}

DiagInfo do_diag(const Step &step, IterInfo &iterinfo, const Coef &coef, Stats &stats, const DiagInfo &diagprev, 
                 QSrmax &qsrmax, shared_ptr<Symmetry> Sym, const Params &P) {
  step.infostring();
  Sym->show_coefficients(step, coef);
  auto tasks = qsrmax.task_list();
  double diagratio = P.diagratio;
  DiagInfo diag;
  while (true) {
    try {
      if (step.nrg()) {
        if (!(P.resume && int(step.ndx()) <= P.laststored))
          diag = diagonalisations(step, iterinfo.opch, coef, diagprev, tasks, diagratio, Sym, P); // compute in first run
        else
          diag = DiagInfo(step.ndx(), false); // or read from disk
      }
      if (step.dmnrg()) {
        diag = DiagInfo(step.ndx(), P.removefiles); // read from disk in second run
        diag.subtract_GS_energy(stats.GS_energy);
      }
      stats.Egs = diag.find_groundstate();
      if (step.nrg()) // should be done only once!
        diag.subtract_Egs(stats.Egs);
      auto cluster_mapping = find_clusters(diag.sorted_energies(), P.fixeps);
      fix_splittings(diag, cluster_mapping);
      truncate_prepare(step, diag, Sym, P);
      break;
    }
    catch (NotEnough &e) {
      fmt::print(fmt::emphasis::bold | fg(fmt::color::yellow), "Insufficient number of states computed.\n");
      if (!(step.nrg() && P.restart)) break;
      diagratio = min(diagratio * P.restartfactor, 1.0);
      fmt::print(fmt::emphasis::bold | fg(fmt::color::yellow), "\nRestarting this iteration step. diagratio={}\n\n", diagratio);
    }
  }
  return diag;
}

// Absolute energies. Must be called in the first NRG run after stats.total_energy has been updated, but before
// store_transformations(). absenergyG is updated to its correct values (referrenced to absolute 0) in
// shift_abs_energies().
void calc_abs_energies(const Step &step, DiagInfo &diag, const Stats &stats) {
  for (auto &eig : diag.eigs()) {
    eig.absenergyN = eig.value_zero * step.scale();        // referenced to the lowest energy in current NRG step (not modified later on)
    eig.absenergy = eig.absenergyN;
    for (auto &x : eig.absenergy) x += stats.total_energy; // absolute energies (not modified later on)
    eig.absenergyG = eig.absenergy;                        // referenced to the absolute 0 (updated by shft_abs_energies())
  }
}

// Perform processing after a successful NRG step. Also called from doZBW() as a final step.
void after_diag(const Step &step, IterInfo &iterinfo, Stats &stats, DiagInfo &diag, Output &output,
                QSrmax &qsrmax, AllSteps &dm, Oprecalc &oprecalc, shared_ptr<Symmetry> Sym, const Params &P) {
  stats.total_energy += stats.Egs * step.scale(); // stats.Egs has already been initialized
  cout << "Total energy=" << HIGHPREC(stats.total_energy) << "  Egs=" << HIGHPREC(stats.Egs) << endl;
  stats.rel_Egs[step.ndx()] = stats.Egs;
  stats.abs_Egs[step.ndx()] = stats.Egs * step.scale();
  stats.energy_offsets[step.ndx()] = stats.total_energy;
  if (step.nrg()) {
    calc_abs_energies(step, diag, stats);  // only in the first run, in the second one the data is loaded from file!
    if (P.dm && !(P.resume && int(step.ndx()) <= P.laststored))
      diag.save(step.ndx());
    perform_basic_measurements(step, diag, Sym, stats, output); // Measurements are performed before the truncation!
  }
  if (!P.ZBW)
    split_in_blocks(diag, qsrmax);
  if (P.do_recalc_all(step.runtype)) { // Either ...
    oprecalc.recalculate_operators(step, diag, qsrmax, iterinfo, Sym, P);
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, dm, Sym, P);
  }
  if (!P.ZBW)
    diag.truncate_perform();                        // Actual truncation occurs at this point
  dm.store(step.ndx(), diag, qsrmax, step.last());  // Store information about subspaces and states for DM algorithms
  if (!step.last()) {
    recalc_irreducible(step, diag, qsrmax, iterinfo.opch, Sym, P);
    if (P.dump_f) iterinfo.opch.dump();
  }
  if (P.do_recalc_kept(step.runtype)) { // ... or ...
    oprecalc.recalculate_operators(step, diag, qsrmax, iterinfo, Sym, P);
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, dm, Sym, P);
  }
  if (P.do_recalc_none())  // ... or this
    calculate_spectral_and_expv(step, stats, output, oprecalc, diag, iterinfo, dm, Sym, P);
  if (P.checksumrules) operator_sumrules(iterinfo, Sym);
}

// Perform one iteration step
DiagInfo iterate(const Step &step, IterInfo &iterinfo, const Coef &coef, Stats &stats, const DiagInfo &diagprev,
                 Output &output, AllSteps &dm, Oprecalc &oprecalc, shared_ptr<Symmetry> Sym, const Params &P) {
  QSrmax qsrmax{diagprev, Sym};
  auto diag = do_diag(step, iterinfo, coef, stats, diagprev, qsrmax, Sym, P);
  after_diag(step, iterinfo, stats, diag, output, qsrmax, dm, oprecalc, Sym, P);
  iterinfo.trim_matrices(diag);
  diag.clear_eigenvectors();
  time_mem::memory_time_brief_report();
  return diag;
}

// Perform calculations with quantities from 'data' file
void docalc0(Step &step, const IterInfo &iterinfo, const DiagInfo &diag0, Stats &stats, Output &output, 
             Oprecalc &oprecalc, shared_ptr<Symmetry> Sym, const Params &P) {
  step.set(P.Ninit - 1); // in the usual case with Ninit=0, this will result in N=-1
  cout << endl << "Before NRG iteration";
  cout << " (N=" << step.N() << ")" << endl;
  perform_basic_measurements(step, diag0, Sym, stats, output);
  AllSteps empty_dm(0, 0);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag0, iterinfo, empty_dm, Sym, P);
  if (P.checksumrules) operator_sumrules(iterinfo, Sym);
}

// doZBW() takes the place of iterate() called from main_loop() in the case of zero-bandwidth calculation.
// It replaces do_diag() and calls after_diag() as the last step.
DiagInfo nrg_ZBW(Step &step, IterInfo &iterinfo, Stats &stats, const DiagInfo &diag0, Output &output, 
                 AllSteps &dm, Oprecalc &oprecalc, shared_ptr<Symmetry> Sym, const Params &P) {
  cout << endl << "Zero bandwidth calculation" << endl;
  step.set_ZBW();
  // --- begin do_diag() equivalent
  DiagInfo diag;
  if (step.nrg()) 
    diag = diag0;
  if (step.dmnrg()) {
    diag = DiagInfo(step.ndx(), P.removefiles);
    diag.subtract_GS_energy(stats.GS_energy);
  }
  stats.Egs = diag.find_groundstate();
  if (step.nrg())      
    diag.subtract_Egs(stats.Egs);
  truncate_prepare(step, diag, Sym, P); // determine # of kept and discarded states
  // --- end do_diag() equivalent
  QSrmax qsrmax{};
  after_diag(step, iterinfo, stats, diag, output, qsrmax, dm, oprecalc, Sym, P);
  return diag;
}

// ****************************  Main NRG loop ****************************

DiagInfo nrg_loop(Step &step, IterInfo &iterinfo, const Coef &coef, Stats &stats, const DiagInfo &diag0, 
                  Output &output, AllSteps &dm, Oprecalc &oprecalc, shared_ptr<Symmetry> Sym, const Params &P) {
  DiagInfo diag = diag0;
  for (step.init(); !step.end(); step.next())
    diag = iterate(step, iterinfo, coef, stats, diag, output, dm, oprecalc, Sym, P);
  step.set(step.lastndx());
  return diag;
}

DiagInfo run_nrg(Step &step, IterInfo &iterinfo, const Coef &coef, Stats &stats, const DiagInfo &diag0, 
                 AllSteps &dm, shared_ptr<Symmetry> Sym, const Params &P) {
  diag0.states_report(Sym->multfnc());
  auto oprecalc = Oprecalc(step.runtype, iterinfo, Sym, P);
  auto output = Output(step.runtype, iterinfo, stats, P);
  // If calc0=true, a calculation of TD quantities is performed before starting the NRG iteration.
  if (step.nrg() && P.calc0 && !P.ZBW)
    docalc0(step, iterinfo, diag0, stats, output, oprecalc, Sym, P);
  DiagInfo diag = P.ZBW ? nrg_ZBW(step, iterinfo, stats, diag0, output, dm, oprecalc, Sym, P) 
    : nrg_loop(step, iterinfo, coef, stats, diag0, output, dm, oprecalc, Sym, P);
  fmt::print(fmt::emphasis::bold | fg(fmt::color::red), FMT_STRING("\nTotal energy: {:.18}\n"), stats.total_energy);
  stats.GS_energy = stats.total_energy;
  if (step.nrg() && P.dumpsubspaces) dm.dump_subspaces();
  cout << endl << "** Iteration completed." << endl << endl;
  return diag;
}

void print_about_message() {
  fmt::print(fmt::emphasis::bold, "NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print(fmt::emphasis::bold, "Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print(fmt::emphasis::bold, "Compiled on {} at {}\n\n", __DATE__, __TIME__);
}

std::unique_ptr<Symmetry> get(const std::string &sym_string, const Params &P, Allfields &allfields)
{
  if (sym_string == "QS")     return std::make_unique<SymmetryQS>(P, allfields);
  if (sym_string == "QSZ")    return std::make_unique<SymmetryQSZ>(P, allfields);  
#ifdef NRG_SYM_MORE
  if (sym_string == "ISO")    return std::make_unique<SymmetryISO>(P, allfields);
  if (sym_string == "ISO2")   return std::make_unique<SymmetryISO2>(P, allfields);
  if (sym_string == "ISOSZ")  return std::make_unique<SymmetryISOSZ>(P, allfields);
  if (sym_string == "SPSU2")  return std::make_unique<SymmetrySPSU2>(P, allfields);
  if (sym_string == "SPU1")   return std::make_unique<SymmetrySPU1>(P, allfields);
#endif
#ifdef NRG_SYM_ALL
  if (sym_string == "DBLISOSZ")  return std::make_unique<SymmetryDBLISOSZ>(P, allfields);
  if (sym_string == "DBLSU2")    return std::make_unique<SymmetryDBLSU2>(P, allfields);
  if (sym_string == "ISOLR")     return std::make_unique<SymmetryISOLR>(P, allfields);
  if (sym_string == "ISO2LR")    return std::make_unique<SymmetryISO2LR>(P, allfields);
  if (sym_string == "ISOSZLR")   return std::make_unique<SymmetryISOSZLR>(P, allfields);
  if (sym_string == "NONE")      return std::make_unique<SymmetryNONE>(P, allfields);
  if (sym_string == "P")         return std::make_unique<SymmetryP>(P, allfields);
  if (sym_string == "PP")        return std::make_unique<SymmetryPP>(P, allfields);
  if (sym_string == "QJ")        return std::make_unique<SymmetryQJ>(P, allfields);
  if (sym_string == "QSLR")      return std::make_unique<SymmetryQSLR>(P, allfields); 
  if (sym_string == "QST")       return std::make_unique<SymmetryQST>(P, allfields);
  if (sym_string == "QSTZ")      return std::make_unique<SymmetryQSTZ>(P, allfields);
  if (sym_string == "QSZLR")     return std::make_unique<SymmetryQSZLR>(P, allfields);
  if (sym_string == "QSZTZ")     return std::make_unique<SymmetryQSZTZ>(P, allfields);
  if (sym_string == "SL")        return std::make_unique<SymmetrySL>(P, allfields);
  if (sym_string == "SL3")       return std::make_unique<SymmetrySL3>(P, allfields);
  if (sym_string == "SPSU2LR")   return std::make_unique<SymmetrySPSU2LR>(P, allfields);
  if (sym_string == "SPSU2T")    return std::make_unique<SymmetrySPSU2T>(P, allfields);
  if (sym_string == "SPU1LR")    return std::make_unique<SymmetrySPU1LR>(P, allfields);
  if (sym_string == "SU2")       return std::make_unique<SymmetrySU2>(P, allfields);
  if (sym_string == "U1")        return std::make_unique<SymmetryU1>(P, allfields);
 #ifdef NRG_COMPLEX
  if (sym_string == "QSC3")      return std::make_unique<SymmetryQSC3>(P, allfields);
  if (sym_string == "SPSU2C3")   return std::make_unique<SymmetrySPSU2C3>(P, allfields);
 #endif
#endif 
  throw std::runtime_error("Unknown symmetry " + sym_string);
}

// Called immediately after parsing the information about the number of channels from the data file. This ensures
// that Invar can be parsed correctly.
shared_ptr<Symmetry> set_symmetry(const Params &P, Stats &stats) {
  my_assert(P.channels > 0 && P.combs > 0); // must be set at this point
  cout << "SYMMETRY TYPE: " << P.symtype.value() << endl;
  auto Sym = get(P.symtype.value(), P, stats.td.allfields);
  Sym->load();
  Sym->erase_first();
  return Sym;
}

void calculation() {
  // Workdir workdir;
  Params P("param", "param", workdir);
  Stats stats(P);
  auto [diag0, iterinfo, coef, Sym] = read_data(P, stats);
  Step step{P, RUNTYPE::NRG};
  AllSteps dm(P.Ninit, P.Nlen);
  auto diag = run_nrg(step, iterinfo, coef, stats, diag0, dm, Sym, P);
  if (string(P.stopafter) == "nrg") exit1("*** Stopped after the first sweep.");
  dm.shift_abs_energies(stats.GS_energy); // we call this here, to enable a file dump
  if (P.dumpabsenergies)
    dm.dump_all_absolute_energies();
  if (P.dm) {
    if (P.need_rho()) {
      auto rho = init_rho(step, diag, Sym);
      rho.save(step.lastndx(), FN_RHO);
      if (!P.ZBW) calc_densitymatrix(rho, dm, Sym, P);
    }
    if (P.need_rhoFDM()) {
      calc_ZnD(dm, stats, Sym, P.T);
      if (P.logletter('w')) 
        report_ZnD(stats, P);
      fdm_thermodynamics(dm, stats, Sym, P.T);
      auto rhoFDM = init_rho_FDM(step.lastndx(), dm, stats, Sym, P.T);
      rhoFDM.save(step.lastndx(), FN_RHOFDM);
      if (!P.ZBW) calc_fulldensitymatrix(step, rhoFDM, dm, stats, Sym, P);
    }
    if (string(P.stopafter) == "rho") exit1("*** Stopped after the DM calculation.");
    auto [diag0_dm, iterinfo_dm, coef_dm, Sym_dm] = read_data(P, stats);
    Step step_dmnrg{P, RUNTYPE::DMNRG};
    run_nrg(step_dmnrg, iterinfo_dm, coef_dm, stats, diag0_dm, dm, Sym_dm, P);
    my_assert(num_equal(stats.GS_energy, stats.total_energy));
  }
  if (P.done) { ofstream D("DONE"); } // Indicate completion by creating a flag file
}

// Master process does most of the i/o and passes calculations to the slaves.
void run_nrg_master() {
  calculation();
#ifdef NRG_MPI
  for (int i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_EXIT, 0);
#endif
}

#ifdef NRG_MPI
// Handle a diagonalisation request:
void slave_diag(const int master, const DiagParams &DP) {
  // 1. receive the matrix and the subspace identification
  auto m = mpi_receive_matrix(master);
  Invar I;
  check_status(mpiw->recv(master, TAG_INVAR, I));
  // 2. preform the diagonalisation
  Eigen eig = diagonalise(m, DP);
  // 3. send back the results
  mpi_send_eigen(master, eig);
  mpiw->send(master, TAG_INVAR, I);
}

void run_nrg_slave() {
  constexpr auto master = 0;
  DiagParams DP;
  for (;;) {
    if (mpiw->iprobe(master, mpi::any_tag)) { // message can be received.
      int task;
      auto status = mpiw->recv(master, mpi::any_tag, task);
      mpilog("Slave " << mpiw->rank() << " received message with tag " << status.tag());
      check_status(status);
      switch (status.tag()) {
        case TAG_SYNC:
          DP = mpi_receive_params();
          break;
        case TAG_DIAG:
          slave_diag(master, DP);
          break;
        case TAG_EXIT:
          return; // exit from run_slave()
        default: 
          cout << "MPI error: unknown tag on " << mpiw->rank() << endl; 
          break;
      }
    } else usleep(100); // sleep to reduce the load on the computer. (OpenMPI "feature" workaround)
  }
}
#else
void run_nrg_slave() {}
#endif
