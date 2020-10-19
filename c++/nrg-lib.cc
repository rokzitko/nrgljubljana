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
using scalar = double;
using t_matel = double;                       // type for the matrix elements
using t_coef = double;                        // type for the Wilson chain coefficients
using t_coef = double;                      // type for various prefactors in recalculations
using t_expv = double;                        // type for expectation values of operators
#endif

#ifdef NRG_COMPLEX
using scalar = cmpl;
using t_matel = cmpl;
using t_coef = cmpl;
using t_coef = cmpl;
using t_expv = cmpl; // we allow the calculation of expectation values of non-Hermitian operators!
#endif

using t_eigen = double;  // type for the eigenvalues (always double)
using t_weight = cmpl;   // spectral weight accumulators (always complex)

template <typename S> struct traits {};

template <> struct traits<double> {
  using t_matel = double;  // type for the matrix elements
  using t_coef = double;   // type for the Wilson chain coefficients & various prefactors
  using t_expv = double;   // type for expectation values of operators
  using t_eigen = double;  // type for the eigenvalues (always double)
  using t_weight = cmpl;   // spectral weight accumulators (always complex)
  using Matrix = ublas::matrix<t_matel>;
};

template <> struct traits<cmpl> {
  using t_matel = cmpl;
  using t_coef = cmpl;
  using t_expv = cmpl;     // we allow the calculation of expectation values of non-Hermitian operators!
  using t_eigen = double;  // type for the eigenvalues (always double)
  using t_weight = cmpl;   // spectral weight accumulators (always complex)
  using Matrix = ublas::matrix<t_matel>;
};

inline cmpl conj_me(const cmpl &z) { return conj(z); } // conjugation
inline double conj_me(const double x) { return x; }    // no op

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
   std::cout << "Memory used: " << long(ms.used() / 1024) << " MB "; // NOLINT
#endif
   std::cout << "Time elapsed: " << prec3(tm.total_in_seconds()) << " s" << std::endl;
 }
}

#ifdef NRG_MPI
mpi::environment *mpienv;
mpi::communicator *mpiw;
int myrank() { return mpiw->rank(); } // used in diag.h, time_mem.h
#else
int myrank() { return 0; }
#endif

// Quantum number types
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
template <typename S> class Eigen_tmpl {
public:
  using t_eigen = typename traits<S>::t_eigen;
  using Matrix2 = typename traits<S>::Matrix;
  using EVEC = ublas::vector<t_eigen>;
  EVEC value_orig;               // eigenvalues as computed
  Matrix2 matrix; // eigenvectors
  Eigen_tmpl() {}
  Eigen_tmpl(const size_t nr, const size_t dim) {
    my_assert(nr <= dim);
    value_orig.resize(nr);
    matrix.resize(nr, dim);
  }
  auto getnrcomputed() const { return value_orig.size(); } // number of computed eigenpairs
  auto getdim() const { return matrix.size2(); }           // valid also after the split_in_blocks_Eigen() call
  // Now add information about truncation and block structure of the eigenvectors
 private:
  long nrpost = -1;  // number of eigenpairs after truncation (-1: keep all)
 public:
  EVEC value_zero;   // eigenvalues with Egs subtracted
  auto getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }
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
  std::vector<Matrix2> blocks;
  // Truncate to nrpost states.
  void truncate_prepare_subspace(const size_t nrpost_) {
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
  void subtract_Egs(const double Egs) {
    value_zero = value_orig;
    for (auto &x : value_zero) x -= Egs;
    my_assert(value_zero[0] >= 0);
  }
  void subtract_GS_energy(const double GS_energy) {
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
using Eigen = Eigen_tmpl<scalar>;

// Full information after diagonalizations (eigenspectra in all subspaces)
template<typename S>
class DiagInfo_tmpl : public std::map<Invar, Eigen_tmpl<S>> {
 public:
   explicit DiagInfo_tmpl() {}
   DiagInfo_tmpl(std::ifstream &fdata, const size_t nsubs, const Params &P) {
     for (const auto i : range1(nsubs)) {
       Invar I;
       fdata >> I;
       auto energies = read_vector<double>(fdata);
       if (!P.data_has_rescaled_energies && !P.absolute)
         energies /= P.SCALE(P.Ninit); // rescale to the suitable energy scale
       (*this)[I].diagonal(energies);
     }
     my_assert(this->size() == nsubs);
   }
   auto subspaces() const { return *this | boost::adaptors::map_keys; }
   auto eigs() const { return *this | boost::adaptors::map_values; }
   auto eigs() { return *this | boost::adaptors::map_values; }
   auto find_groundstate() const {
     const auto [Iground, eig] = *ranges::min_element(*this, [](const auto a, const auto b) { return a.second.value_orig(0) < b.second.value_orig(0); });
     const auto Egs = eig.value_orig(0);
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
   void dump_value_zero(std::ostream &F) const {
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
     const std::string fn = workdir.unitaryfn(N);
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
     const std::string fn = workdir.unitaryfn(N);
     std::ifstream MATRIXF(fn, std::ios::binary | std::ios::in);
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
   explicit DiagInfo_tmpl(const size_t N, const bool remove_files = false) { load(N, remove_files); }
};
using DiagInfo = DiagInfo_tmpl<scalar>;

template<typename S>
class MatrixElements_tmpl : public std::map<Twoinvar, typename traits<S>::Matrix> {
 public:
   MatrixElements_tmpl() {}
   MatrixElements_tmpl(std::ifstream &fdata, const DiagInfo_tmpl<S> &diag) {
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
       os << "----" << II << "----" << std::endl << mat << std::endl;
     return os;
   }
};
template<typename S>
std::ostream &operator<<(std::ostream &os, const MatrixElements_tmpl<S> &m) { return m.insertor(os); }
using MatrixElements = MatrixElements_tmpl<scalar>;

template<typename S>
class DensMatElements_tmpl : public std::map<Invar, typename traits<S>::Matrix> {
 public:
   template <typename MF>
     auto trace(MF mult) const {
       return ranges::accumulate(*this, 0.0, [mult](double acc, const auto z) { const auto &[I, mat] = z; 
         return acc + mult(I) * trace_real(mat); });
     }
   void save(const size_t N, const std::string &prefix) const {
     const auto fn = workdir.rhofn(prefix, N);
     std::ofstream MATRIXF(fn, std::ios::binary | std::ios::out);
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
   void load(const size_t N, const string &prefix, const bool remove_files) {
     const auto fn = workdir.rhofn(prefix, N);
     std::ifstream MATRIXF(fn, std::ios::binary | std::ios::in);
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
using DensMatElements = DensMatElements_tmpl<scalar>;

// Map of operators matrices
template<typename S>
using CustomOp_tmpl = std::map<std::string, MatrixElements_tmpl<S>>;
using CustomOp = CustomOp_tmpl<scalar>;

// Vector containing irreducible matrix elements of f operators.
template<typename S>
using OpchChannel_tmpl = std::vector<MatrixElements_tmpl<S>>;
using OpchChannel = OpchChannel_tmpl<scalar>;

// Each channel contains P.perchannel OpchChannel matrices.
template<typename S>
class Opch_tmpl : public std::vector<OpchChannel_tmpl<S>> {
 public:
   Opch_tmpl() {}
   explicit Opch_tmpl(const size_t nrch) { this->resize(nrch); }
   Opch_tmpl(std::ifstream &fdata, const DiagInfo_tmpl<S> &diag, const Params &P) {
     this->resize(P.channels);
     for (const auto i : range0(size_t(P.channels))) {
       (*this)[i] = OpchChannel_tmpl<S>(P.perchannel);
       for (const auto j : range0(size_t(P.perchannel))) {
         char ch;
         size_t iread, jread;
         fdata >> ch >> iread >> jread;
         my_assert(ch == 'f' && i == iread && j == jread);
         (*this)[i][j] = MatrixElements_tmpl<S>(fdata, diag);
       }
     }
   }
   void dump() {
     std::cout << std::endl;
     for (const auto &&[i, ch] : *this | ranges::views::enumerate)
       for (const auto &&[j, mat] : ch | ranges::views::enumerate)
         std::cout << fmt::format("<f> dump, i={} j={}\n", i, j) << mat << std::endl;
     std::cout << std::endl;
   }
};
using Opch = Opch_tmpl<scalar>;

template<typename S> class Symmetry_tmpl;
using Symmetry = Symmetry_tmpl<scalar>;

// Dimensions of the invariant subspaces |r,1>, |r,2>, |r,3>, etc. The name "rmax" comes from the maximal value of
// the index "r" which ranges from 1 through rmax.

class Rmaxvals {
 private:
   std::vector<size_t> values;
   std::shared_ptr<Symmetry> Sym;
 public:
   Rmaxvals() = default;
   template<typename S>
     Rmaxvals(const Invar &I, const InvarVec &In, const DiagInfo_tmpl<S> &diagprev, shared_ptr<Symmetry> Sym);
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
   // *** Mathematica interfacing: i1,j1 are 1-based
   bool offdiag_contributes(const size_t i1, const size_t j1) const { // i,j are 1-based (Mathematica interface)
     my_assert(1 <= i1 && i1 <= combs() && 1 <= j1 && j1 <= combs());
     my_assert(i1 != j1);
     return exists(i1-1) && exists(j1-1); // shift by 1
   }
   auto chunk(const size_t i1) const {
     return std::make_pair(offset(i1-1), rmax(i1-1));
   }
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
   template<typename S> QSrmax(const DiagInfo_tmpl<S> &, shared_ptr<Symmetry>);
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
     std::cout << "Stats: nr=" << nr << " min=" << min_size << " max=" << max_size << std::endl;
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
template<typename S> struct DimSub_tmpl {
  size_t kept  = 0;
  size_t total = 0;
  Rmaxvals rmax;
  Eigen_tmpl<S> eig;
  bool is_last = false;
  auto min() const { return is_last ? 0 : kept; } // min(), max() return the range of D states to be summed over in FDM
  auto max() const { return total; }
  auto all() const { return boost::irange(min(), max()); }
};
using DimSub = DimSub_tmpl<scalar>;

// Full information about the number of states and matrix dimensions
// Example: dm[N].rmax[I] etc.
template<typename S>
using Subs = std::map<Invar, DimSub_tmpl<S>>;

template<typename S>
class AllSteps_tmpl : public std::vector<Subs<S>> {
 public:
   const size_t Nbegin, Nend; // range of valid indexes
   AllSteps_tmpl(const size_t Nbegin, const size_t Nend) : Nbegin(Nbegin), Nend(Nend) { this->resize(Nend ? Nend : 1); } // at least 1 for ZBW
   auto Nall() const { return boost::irange(Nbegin, Nend); }
   void dump_absenergyG(std::ostream &F) const {
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
     std::ofstream O(filename);
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
using AllSteps = AllSteps_tmpl<scalar>;

class Step {
 private:
   // N denotes the order of the Hamiltonian. N=0 corresponds to H_0, i.e. the initial Hamiltonian
   int trueN; // "true N", sets the energy scale, it may be negative, trueN <= ndxN
   size_t ndxN; // "index N", iteration step, used as an array index, ndxN >= 0
   const Params &P; // reference to parameters (beta, T)
   
 public:
   const RUNTYPE runtype; // NRG vs. DM-NRG run
   void set(const int newN) {
     trueN = newN;
     ndxN = std::max(newN, 0);
   }
   void init() { set(P.Ninit); }
   Step(const Params &P_, const RUNTYPE runtype_) : P(P_), runtype(runtype_) { init(); }
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
template<typename S>
class Stats_tmpl {
 public:
   using t_eigen = typename traits<S>::t_eigen;
   using t_expv  = typename traits<S>::t_expv;
   t_eigen Egs{};
   
   // ** Thermodynamic quantities
   double Z{};
   double Zft{};   // grand-canonical partition function (at shell n)
   double Zgt{};   // grand-canonical partition function for computing G(T)
   double Zchit{}; // grand-canonical partition function for computing chi(T)

   TD td;
   
   //  ** Expectation values
   std::map<std::string, t_expv> expv;    // expectation values of custom operators
   std::map<std::string, t_expv> fdmexpv; // Expectation values computed using the FDM algorithm
   
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

   explicit Stats_tmpl(const Params &P, const std::string filename_td = "td"s, const std::string filename_tdfdm = "tdfdm"s) : 
     td(P, filename_td), rel_Egs(MAX_NDX), abs_Egs(MAX_NDX), energy_offsets(MAX_NDX), 
     ZnDG(MAX_NDX), ZnDN(MAX_NDX), ZnDNd(MAX_NDX), wn(MAX_NDX), wnfactor(MAX_NDX), td_fdm(P, filename_tdfdm) {}
};
using Stats = Stats_tmpl<scalar>;

// Wrapper class for NRG spectral-function algorithms
class Algo {
 private:
 public:
   const Params &P;
   Algo() = delete;
   Algo(const Algo&) = delete;
   explicit Algo(const Params &P) : P(P) {}
   virtual ~Algo() = default;
   virtual void begin(const Step &) {}
   virtual void calc(const Step &, const Eigen &, const Eigen &, const Matrix &, const Matrix &, 
                     const t_coef, const Invar &, const Invar &, const DensMatElements &, const Stats &stats) {} // XXX: =0 ?
   virtual void end(const Step &) {}
   virtual std::string rho_type() { return ""; } // what rho type is required
};

template<typename M> void dump_diagonal_matrix(const ublas::matrix<M> &m, const size_t max_nr, std::ostream &F) {
  for (const auto r : range0(std::min(m.size1(), max_nr)))
    F << m(r,r) << ' ';
  F << std::endl;
}

template<typename S>
void dump_diagonal_op(const std::string &name, const MatrixElements_tmpl<S> &n, const size_t max_nr, std::ostream &F) {
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
template<typename S> void trim_matel(const DiagInfo_tmpl<S> &diag, MatrixElements_tmpl<S> &op) {
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
    ublas::matrix_range<typename traits<S>::Matrix> m2(mat, ublas::range(0, nr1), ublas::range(0, nr2));
    typename traits<S>::Matrix m2new = m2;
    mat.swap(m2new);
  }
}

template<typename S> void trim_op(const DiagInfo_tmpl<S> &diag, CustomOp_tmpl<S> &allops) {
  for (auto &[name, op] : allops) 
    trim_matel(diag, op);
}

// Object of class IterInfo cotains full information about matrix representations when entering stage N of the NRG
// iteration.
template<typename S> class IterInfo_tmpl {
 public:
   Opch_tmpl<S> opch;     // f operators (channels)
   CustomOp_tmpl<S> ops;  // singlet operators (even parity)
   CustomOp_tmpl<S> opsp; // singlet operators (odd parity)
   CustomOp_tmpl<S> opsg; // singlet operators [global op]
   CustomOp_tmpl<S> opd;  // doublet operators (spectral functions)
   CustomOp_tmpl<S> opt;  // triplet operators (dynamical spin susceptibility)
   CustomOp_tmpl<S> opq;  // quadruplet operators (spectral functions for J=3/2)
   CustomOp_tmpl<S> opot; // orbital triplet operators

   void dump_diagonal(const size_t max_nr, std::ostream &F = std::cout) const {
     if (max_nr) {
       for (const auto &[name, m] : ops)  dump_diagonal_op(name, m, max_nr, F);
       for (const auto &[name, m] : opsg) dump_diagonal_op(name, m, max_nr, F);
     }
   }
   void trim_matrices(const DiagInfo_tmpl<S> &diag) {
     trim_op(diag, ops);
     trim_op(diag, opsp);
     trim_op(diag, opsg);
     trim_op(diag, opd);
     trim_op(diag, opt);
     trim_op(diag, opot);
     trim_op(diag, opq);
   }
};
using IterInfo = IterInfo_tmpl<scalar>;

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
template<typename S, typename F> 
double norm(const MatrixElements_tmpl<S> &m, shared_ptr<Symmetry> Sym, F factor_fnc, const int SPIN) {
  weight_bucket sum;
  for (const auto &[II, mat] : m) {
    const auto & [I1, Ip] = II;
    if (!Sym->check_SPIN(I1, Ip, SPIN)) continue;
    sum += factor_fnc(Ip, I1) * frobenius_norm(mat);
  }
  return 2.0 * cmpl(sum).real(); // Factor 2: Tr[d d^\dag + d^\dag d] = 2 \sum_{i,j} A_{i,j}^2 !!
}

template<typename S>
void operator_sumrules(const IterInfo_tmpl<S> &a, shared_ptr<Symmetry> Sym) {
  // We check sum rules wrt some given spin (+1/2, by convention). For non-spin-polarized calculations, this is
  // irrelevant (0).
  const int SPIN = Sym->isfield() ? 1 : 0;
  for (const auto &[name, m] : a.opd)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensFactorFnc(), SPIN) << std::endl;
  for (const auto &[name, m] : a.opq)
    std::cout << "norm[" << name << "]=" << norm(m, Sym, Sym->SpecdensquadFactorFnc(), 0) << std::endl;
}

#include "read-input.cc"

template <typename T>
  std::string formatted_output(const T x, const Params &P) {
    return fmt::format("{x:>{width}}", "x"_a=x, "width"_a=P.width_custom);
  }

std::string formatted_output(const double x, const Params &P) {
  return fmt::format("{x:>{width}.{prec}}", "x"_a=x, "prec"_a=P.prec_custom, "width"_a=P.width_custom);
}

template <typename T>
  bool negligible_imag_part(const std::complex<T> &z, const double output_imag_eps = 1e-13) {
    return abs(z.imag()) < abs(z.real()) * output_imag_eps; 
  }

// The output format for complex values is X+IY or X-IY, where X and Y are real and imaginary part, respectively. The
// imaginary part is only shown where its value relative to the real part is sufficiently large. No space is used in
// the outputted string in order to simplify parsing.

std::string formatted_output(const cmpl z, const Params &P) {
  const auto [r, i] = reim(z);
  const auto str = P.noimag || negligible_imag_part(z) ?
    fmt::format("{r:.{prec}f}", "r"_a=r, "prec"_a=P.prec_custom) :
    fmt::format("{r:.{prec}f}{s}I{absi:.{prec}f}", "r"_a=r, "s"_a=(i>0 ? "+" : "-"), "absi"_a=abs(i), "prec"_a=P.prec_custom);
  return fmt::format("{str:>{width}}", "str"_a=str, "width"_a=P.width_custom); // the width for the whole X+iY string
}

// #### GF's

enum class gf_type { bosonic, fermionic };

std::string gf_typestring(gf_type gt) { // XXX: used anywhere?
  switch (gt) {
  case gf_type::bosonic: return "bosonic";
  case gf_type::fermionic: return "fermionic";
  default: my_assert_not_reached();
  }
}

// Sign factor in GFs for bosonic/fermionic operators.
inline constexpr auto S_FERMIONIC = -1;
inline constexpr auto S_BOSONIC   = 1;
int gf_sign(const gf_type gt)
{
  return gt == gf_type::bosonic ? S_BOSONIC : S_FERMIONIC;
}

#include "bins.h"
#include "matsubara.h"

class ChainBinning {
 private:
   const Params &P;
   Bins spos, sneg;
 public:
   explicit ChainBinning(const Params &P) : P(P), spos(P), sneg(P) {}
   void add(const double energy, const t_weight weight) {
     if (energy >= 0.0)
       spos.add(energy, weight);
     else
       sneg.add(-energy, weight);
   }
   auto total_weight() const { return spos.total_weight() + sneg.total_weight(); }
   friend class SpectrumRealFreq;
};

class ChainMatsubara {
 private:
   const Params &P;
   Matsubara m;
 public:
   explicit ChainMatsubara(const Params &P, const gf_type gt) : P(P), m(P.mats, gt, P.T){};
   void add(const size_t n, const t_weight w) { m.add(n, w); }
   auto total_weight() const { return m.total_weight(); }
   friend class GFMatsubara;
};

class ChainTempDependence {
 private:
   const Params &P;
   Temp v;
 public:
   explicit ChainTempDependence(const Params &P) : P(P), v(P) {}
   void add(const double T, const t_weight value) { v.add_value(T, value); }
   friend class TempDependence;
};

#include "spectrumrealfreq.cc"

class GFMatsubara {
 private:
   const std::string name, algoname, filename;
   const Params &P;
   Matsubara results;
 public:
   GFMatsubara(const string &name, const std::string &algoname, const std::string &filename, gf_type gt, const Params &P) : 
     name(name), algoname(algoname), filename(filename), P(P), results(P.mats, gt, P.T) {}
   void merge(const ChainMatsubara &cm) {
     results.merge(cm.m);
   }     
   ~GFMatsubara() {
     fmt::print(fmt::emphasis::bold, "GF Matsubara: {} {} -> {}\n", name, algoname, filename);
     results.save(safe_open(filename + ".dat"), P.prec_xy);
   }
};

class TempDependence {
 private:
   const std::string name, algoname, filename;
   const Params &P;
   Spikes results;
 public:
   TempDependence(const std::string &name, const std::string &algoname, const std::string &filename, const Params &P) : 
     name(name), algoname(algoname), filename(filename),  P(P) {}
   void merge(const ChainTempDependence &ctd) {
     std::copy(ctd.v.begin(), ctd.v.end(), std::back_inserter(results));
   }
   ~TempDependence() {
     fmt::print(fmt::emphasis::bold, "Temperature dependence: {} {} -> {}\n", name, algoname, filename);
     ranges::sort(results, sortfirst());
     results.save(safe_open(filename + ".dat"), P.prec_xy, P.reim);
   }
};


// Check if the trace of the density matrix equals 'ref_value'.
template<typename S>
void check_trace_rho(const DensMatElements_tmpl<S> &m, std::shared_ptr<Symmetry> Sym, const double ref_value = 1.0) {
  if (!num_equal(m.trace(Sym->multfnc()), ref_value))
    throw std::runtime_error("check_trace_rho() failed");
}

enum class axis { RealFreq, Temp, Matsubara };

std::ostream & operator<<(std::ostream &os, const axis &a) { // XXX: keep this?
  if (a == axis::RealFreq)  os << "RealFreq";
  if (a == axis::Temp)      os << "Temp";
  if (a == axis::Matsubara) os << "Matsubara";
  return os;
}

inline std::string to_string(std::complex<double> &z) { // XXX: i/o?
  std::ostringstream s;
  s << z;
  return s.str();
}

// All information about calculating a spectral function: pointers to the operator data, raw spectral data
// acccumulators, algorithm, etc.
class BaseSpectrum {
 public:
   const MatrixElements &op1, &op2;
   int spin{};                      // -1 or +1, or 0 where irrelevant
   std::shared_ptr<Algo> algo;      // Algo_FDM, Algo_DMNRG,...
   BaseSpectrum(const MatrixElements &op1, const MatrixElements &op2, const int spin) :
     op1(op1), op2(op2), spin(spin) {}
};
using speclist = std::list<BaseSpectrum>;

#include "spec.cc"
#include "dmnrg.h"
#include "splitting.cc"

// Determine the ranges of index r
template<typename S>
Rmaxvals::Rmaxvals(const Invar &I, const InvarVec &InVec, const DiagInfo_tmpl<S> &diagprev, std::shared_ptr<Symmetry> Sym) {
  for (const auto &[i, In] : InVec | ranges::views::enumerate)
    values.push_back(Sym->triangle_inequality(I, In, Sym->QN_subspace(i)) ? diagprev.size_subspace(In) : 0);
}

// Formatted output of the computed expectation values
template<typename S>
class ExpvOutput_tmpl {
 private:
   using t_expv = typename traits<S>::t_expv;
   std::ofstream F;                     // output stream
   std::map<std::string, t_expv> &m;    // reference to the name->value mapping
   const std::list<std::string> fields; // list of fields to be output (may be a subset of the fields actually present in m)
   const Params &P;
   void field_numbers() {     // Consecutive numbers for the columns
     F << '#' << formatted_output(1, P) << ' ';
     for (const auto ctr : range1(fields.size())) F << formatted_output(1 + ctr, P) << ' ';
     F << std::endl;
   }
   // Label and field names. Label is the first column (typically the temperature).
   void field_names(const std::string labelname = "T") {
     F << '#' << formatted_output(labelname, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(op, P); });
     F << std::endl;
   }
 public:
   // Output the current values for the label and for all the fields
   void field_values(const double labelvalue, const bool cout_dump = true) {
     F << ' ' << formatted_output(labelvalue, P) << ' ';
     std::transform(fields.cbegin(), fields.cend(), std::ostream_iterator<std::string>(F, " "), [this](const auto op) { return formatted_output(m[op], P); });
     F << std::endl;
     if (cout_dump)
       for (const auto &op: fields)
         fmt::print(fmt::emphasis::bold | fg(fmt::color::red), "<{}>={}\n", op, to_string(m[op]));
   }
   ExpvOutput_tmpl(const string &fn, map<string, t_expv> &m_, const list<string> &fields_, const Params &P_) : m(m_), fields(fields_), P(P_) {
     F.open(fn);
     field_numbers();
     field_names();
   }
};
using ExpvOutput = ExpvOutput_tmpl<scalar>;

// Establish the data structures for storing spectral information [and prepare output files].
template<typename M>
void prepare_spec_algo(speclist &sl, M && op1, M && op2, int spin, std::string name, std::string prefix, std::string algoname, gf_type gt, const Params &P) {
  fmt::print("Spectrum: {} {} {}\n", name, prefix, algoname);
  const auto filename = prefix + "_" + algoname + "_dens_" + name; // no suffix (.dat vs. .bin)
  const auto sign = gf_sign(gt);
  BaseSpectrum spec(std::forward<M>(op1), std::forward<M>(op2), spin);
  if (algoname == "FT")
    spec.algo = std::make_shared<Algo_FT>(SpectrumRealFreq(name,algoname,filename,P), sign, P);
  sl.push_back(spec);
}

template<typename M>
void prepare_spec(const RUNTYPE &runtype, speclist &sl, M && op1, M && op2, std::string name, std::string prefix, gf_type gt, int spin, const Params &P) { 
  // If we did not return from this funciton by this point, what we are computing is the spectral function. There are
  // several possibilities in this case, all of which may be enabled at the same time.
  if (runtype == RUNTYPE::NRG) {
    if (P.finite) prepare_spec_algo(sl, std::forward<M>(op1), std::forward<M>(op2), spin, name, prefix, "FT", gt, P); // XXX make_shared drugje!
  }
}

/*
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
*/

template <typename T> ostream & operator<<(ostream &os, const std::set<T> &x) {
  std::copy(x.cbegin(), x.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}

template<typename S>
class Oprecalc_tmpl {
 public:
   // The following lists hold the names of operators which need to be recomputed. The default behavior is to
   // recompute all the operators that are required to calculate the requested spectral densities, see function
   // open_files(). In addition, singlet operators are always recomputed in the first NRG run, so that we can
   // calculate the expectation values.
   std::set<std::string> s, p, g, d, v, t, q, ot;

   speclist spectraD, spectraS, spectraT, spectraQ, spectraGT, spectraI1T, spectraI2T, spectraK, spectraCHIT, spectraC, spectraOT;

   // Calculate spectral densities
   void spectral_densities(const Step &step, const DiagInfo_tmpl<S> &diag, DensMatElements_tmpl<S> &rho, DensMatElements_tmpl<S> &rhoFDM, 
                           const Stats &stats, shared_ptr<Symmetry> Sym) {
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

   void report(std::ostream &F, const std::string &name, const std::set<std::string> &x) {
     F << name << "=[" << x << "]" << std::endl;
   }

   void report(std::ostream &F = std::cout) {
     F << std::endl << "Computing the following operators:" << std::endl;
     report(F, "s", s);
     report(F, "p", p);
     report(F, "g", g);
     report(F, "d", d);
     report(F, "v", v);
     report(F, "t", t);
     report(F, "q", q);
     report(F, "ot", ot);
   }
   
   bool do_s(const std::string &name, const Params &P, const Step &step) {
     if (step.nrg()) return true;                                          // for computing <O> 
     if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
     return s.count(name);
   }
   
   bool do_g(const std::string &name, const Params &P, const Step &step) {
     if (step.nrg()) return true;                                          // for computing <O>
     if (step.dmnrg() && P.fdmexpv && step.N() <= P.fdmexpvn) return true; // for computing <O> using FDM algorithm
     return g.count(name);
   }
   
   // Wrapper routine for recalculations
   template <typename RecalcFnc>
     MatrixElements_tmpl<S> recalc_common(const MatrixElements_tmpl<S> &mold, RecalcFnc recalc_fnc, const Step &step, const DiagInfo_tmpl<S> &diag,
                                  const QSrmax &qsrmax, const std::string name, const std::string &tip, shared_ptr<Symmetry> Sym, const Params &P) {
       nrglog('0', "Recalculate " << tip << " " << name);
       auto mnew = recalc_fnc(diag, qsrmax, mold);
       if (tip == "g") Sym->recalc_global(step, diag, qsrmax, name, mnew);
       return mnew;
     }
   
   template <typename ... Args>
     MatrixElements_tmpl<S> recalc_or_clear(bool recalc, Args&& ... args) {
       return recalc ? recalc_common(std::forward<Args>(args)...) : MatrixElements_tmpl<S>();
     }

   // Recalculate operator matrix representations
   ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO // avoid false positives
     void recalculate_operators(const Step &step, const DiagInfo_tmpl<S> &diag, const QSrmax &qsrmax, 
                                IterInfo_tmpl<S> &a, shared_ptr<Symmetry> Sym, const Params &P) {
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
   auto sdname(const std::string &a, const std::string &b, const int spin = 0) {
     return a + "-" + b + (spin == 0 ? "" : (spin == 1 ? "-u" : "-d"));
   }

   void loopover(const RUNTYPE &runtype, const Params &P,
                 const CustomOp_tmpl<S> &set1, const CustomOp_tmpl<S> &set2,
                 const string_token &stringtoken, speclist &spectra, const std::string &prefix,
                 std::set<std::string> &rec1, std::set<std::string> &rec2, gf_type mt, const int spin = 0) { // mt -> gt
    for (const auto &[name1, op1] : set1) {
      for (const auto &[name2, op2] : set2) {
        if (const auto name = sdname(name1, name2, spin); stringtoken.find(name)) {
          prepare_spec(runtype, spectra, op1, op2, name, prefix, mt, spin, P);
          rec1.insert(name1);
          rec2.insert(name2);
        }
      }
    }
  }

  // Reset lists of operators which need to be iterated
  Oprecalc_tmpl(const RUNTYPE &runtype, const IterInfo_tmpl<S> &a, shared_ptr<Symmetry> Sym, const Params &P) {
    std::cout << std::endl << "Computing the following spectra:" << std::endl;
    // Correlators (singlet operators of all kinds)
    string_token sts(P.specs);
    loopover(runtype, P, a.ops,  a.ops,  sts, spectraS, "corr", s, s, gf_type::bosonic);
    loopover(runtype, P, a.opsp, a.opsp, sts, spectraS, "corr", p, p, gf_type::bosonic);
    loopover(runtype, P, a.opsg, a.opsg, sts, spectraS, "corr", g, g, gf_type::bosonic);
    loopover(runtype, P, a.ops,  a.opsg, sts, spectraS, "corr", s, g, gf_type::bosonic);
    loopover(runtype, P, a.opsg, a.ops,  sts, spectraS, "corr", g, s, gf_type::bosonic);
    // Global susceptibilities (global singlet operators)
    string_token stchit(P.specchit);
    loopover(runtype, P, a.ops,  a.ops,  stchit, spectraCHIT, "chit", s, s, gf_type::bosonic);
    loopover(runtype, P, a.ops,  a.opsg, stchit, spectraCHIT, "chit", s, g, gf_type::bosonic);
    loopover(runtype, P, a.opsg, a.ops,  stchit, spectraCHIT, "chit", g, s, gf_type::bosonic);
    loopover(runtype, P, a.opsg, a.opsg, stchit, spectraCHIT, "chit", g, g, gf_type::bosonic);
    // Dynamic spin susceptibilities (triplet operators)
    string_token stt(P.spect);
    loopover(runtype, P, a.opt, a.opt, stt, spectraT, "spin", t, t, gf_type::bosonic);
    string_token stot(P.specot);
    loopover(runtype, P, a.opot, a.opot, stot, spectraOT, "orbspin", ot, ot, gf_type::bosonic);
    const int varmin = (Sym->isfield() ? -1 : 0);
    const int varmax = (Sym->isfield() ? +1 : 0);
    // Spectral functions (doublet operators)
    string_token std(P.specd);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, std, spectraD, "spec", d, d, gf_type::fermionic, SPIN);
    string_token stgt(P.specgt);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, stgt, spectraGT, "gt", d, d, gf_type::fermionic, SPIN);
    string_token sti1t(P.speci1t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, sti1t, spectraI1T, "i1t", d, d, gf_type::fermionic, SPIN);
    string_token sti2t(P.speci2t);
    for (int SPIN = varmin; SPIN <= varmax; SPIN += 2)
      loopover(runtype, P, a.opd, a.opd, sti2t, spectraI2T, "i2t", d, d, gf_type::fermionic, SPIN);
    // Spectral functions (quadruplet operators)
    string_token stq(P.specq);
    loopover(runtype, P, a.opq, a.opq, stq, spectraQ, "specq", q, q, gf_type::fermionic);
    report();
  }
};
using Oprecalc = Oprecalc_tmpl<scalar>;

// Store eigenvalue & quantum numbers information (RG flow diagrams)
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
   template<typename S> void dump(const Step &step, const DiagInfo_tmpl<S> &diag, const Stats &stats, 
                                  shared_ptr<Symmetry> Sym, const std::string filename = "annotated.dat") {
     if (!P.dumpannotated) return;
     if (!F.is_open()) { // open output file
       F.open(filename);
       F << std::setprecision(P.dumpprecision);
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
         F << " [" << total_degeneracy << "]" << std::endl;
       }
     } else {
       seznam.resize(len); // truncate!
       for (const auto &[e, I] : seznam) 
         F << scale(e) << " " << I << std::endl;
     }
     F << std::endl; // Consecutive iterations are separated by an empty line
   }
};

// Handle all output
template<typename S>
struct Output_tmpl {
  const RUNTYPE runtype;
  const Params &P;
  Annotated annotated;
  std::ofstream Fenergies;  // all energies (different file for NRG and for DMNRG)
  std::unique_ptr<ExpvOutput_tmpl<S>> custom;
  std::unique_ptr<ExpvOutput_tmpl<S>> customfdm;
  Output_tmpl(const RUNTYPE &runtype, const IterInfo_tmpl<S> &iterinfo, Stats_tmpl<S> &stats, const Params &P,
              const std::string filename_energies= "energies.nrg"s,
              const std::string filename_custom = "custom", 
              const std::string filename_customfdm = "customfdm")
    : runtype(runtype), P(P), annotated(P) {
      // We dump all energies to separate files for NRG and DM-NRG runs. This is a very convenient way to check if both
      // runs produce the same results.
      if (P.dumpenergies && runtype == RUNTYPE::NRG) Fenergies.open(filename_energies);
      std::list<std::string> ops;
      for (const auto &name : iterinfo.ops  | boost::adaptors::map_keys) ops.push_back(name);
      for (const auto &name : iterinfo.opsg | boost::adaptors::map_keys) ops.push_back(name);
      if (runtype == RUNTYPE::NRG)
        custom = std::make_unique<ExpvOutput_tmpl<S>>(filename_custom, stats.expv, ops, P);
      else if (runtype == RUNTYPE::DMNRG && P.fdmexpv) 
        customfdm = std::make_unique<ExpvOutput_tmpl<S>>(filename_customfdm, stats.fdmexpv, ops, P);
    }
  // Dump all energies in diag to a file
  void dump_all_energies(const DiagInfo_tmpl<S> &diag, const int N) {
    if (!Fenergies) return;
    Fenergies << std::endl << "===== Iteration number: " << N << std::endl;
    diag.dump_value_zero(Fenergies);
  }
};
using Output = Output_tmpl<scalar>;

template<typename S>
CONSTFNC auto calc_trace_singlet(const Step &step, const DiagInfo_tmpl<S> &diag, const MatrixElements_tmpl<S> &n, std::shared_ptr<Symmetry> Sym) {
  typename traits<S>::t_matel tr{};
  for (const auto &[I, eig] : diag) {
    const auto & nI = n.at({I,I});
    const auto dim = eig.getnrstored();
    my_assert(dim == nI.size2());
    typename traits<S>::t_matel sum{};
    for (const auto r : range0(dim)) sum += exp(-step.TD_factor() * eig.value_zero(r)) * nI(r, r);
    tr += Sym->mult(I) * sum;
  }
  return tr; // note: t_expv = t_matel
}

// Measure thermodynamic expectation values of singlet operators
template<typename S>
void measure_singlet(const Step &step, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diag, const IterInfo_tmpl<S> &a, 
                     Output_tmpl<S> &output, std::shared_ptr<Symmetry> Sym, const Params &P) {
  const auto Z = ranges::accumulate(diag, 0.0, [&Sym, &step](auto total, const auto &d) { const auto &[I, eig] = d;
    return total + Sym->mult(I) * ranges::accumulate(eig.value_zero, 0.0,
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

template<typename S>
CONSTFNC auto calc_trace_fdm_kept(const size_t ndx, const MatrixElements_tmpl<S> &n, const DensMatElements_tmpl<S> &rhoFDM, 
                                  const AllSteps_tmpl<S> &dm, std::shared_ptr<Symmetry> Sym) {
  typename traits<S>::t_matel tr{};
  for (const auto &[I, rhoI] : rhoFDM)
    tr += Sym->mult(I) * trace_contract(rhoI, n.at({I,I}), dm[ndx].at(I).kept); // over kept states ONLY
  return tr;
}

template<typename S>
void measure_singlet_fdm(const Step &step, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diag, const IterInfo_tmpl<S> &a, Output_tmpl<S> &output, 
                         const DensMatElements_tmpl<S> &rhoFDM, const AllSteps_tmpl<S> &dm, shared_ptr<Symmetry> Sym, const Params &P) {
  for (const auto &[name, m] : a.ops)  stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, dm, Sym);
  for (const auto &[name, m] : a.opsg) stats.fdmexpv[name] = calc_trace_fdm_kept(step.N(), m, rhoFDM, dm, Sym);
  output.customfdm->field_values(P.T);
}

// Calculate grand canonical partition function at current NRG energy shell. This is not the same as the true
// partition function of the full problem! Instead this is the Z_N that is used to initialize the density matrix,
// i.e. rho = 1/Z_N \sum_{l} exp{-beta E_l} |l;N> <l;N|.  grand_canonical_Z() is also used to calculate stats.Zft,
// that is used to compute the spectral function with the conventional approach, as well as stats.Zgt for G(T)
// calculations, stats.Zchit for chi(T) calculations.
template<typename S>
auto grand_canonical_Z(const Step &step, const DiagInfo_tmpl<S> &diag, std::shared_ptr<Symmetry> Sym, const double factor = 1.0) {
  double ZN{};
  for (const auto &[I, eig]: diag) 
    for (const auto &i : eig.kept()) // sum over all kept states
      ZN += Sym->mult(I) * exp(-eig.value_zero(i) * step.scT() * factor);
  my_assert(ZN >= 1.0);
  return ZN;
}

template<typename S> auto diagonal_exp(const Eigen_tmpl<S> &eig, const double factor) {
  const auto dim = eig.getnrstored();
  typename traits<S>::Matrix m(dim, dim, 0);
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
template<typename S>
auto init_rho(const Step &step, const DiagInfo_tmpl<S> &diag, std::shared_ptr<Symmetry> Sym) {
  DensMatElements_tmpl<S> rho;
  for (const auto &[I, eig]: diag)
    rho[I] = diagonal_exp(eig, step.scT()) / grand_canonical_Z(step, diag, Sym);
  check_trace_rho(rho, Sym);
  return rho;
}

// Determine the number of states to be retained. Returns Emax - the highest energy to still be retained.
template<typename S>
auto highest_retained_energy(const Step &step, const DiagInfo_tmpl<S> &diag, const Params &P) {
  auto energies = diag.sorted_energies();
  my_assert(energies.front() == 0.0); // check for the subtraction of Egs
  const auto totalnumber = energies.size();
  size_t nrkeep;
  if (P.keepenergy <= 0.0) {
    nrkeep = P.keep;
  } else {
    const auto keepenergy = P.keepenergy * step.unscale();
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
  template<typename S>
  truncate_stats(const DiagInfo_tmpl<S> &diag, shared_ptr<Symmetry> Sym) {
    nrall = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getdim(); });
    nrallmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getdim(); });
    nrkept = ranges::accumulate(diag, 0, [](int n, const auto &d) { const auto &[I, eig] = d; return n+eig.getnrkept(); });
    nrkeptmult = ranges::accumulate(diag, 0, [Sym](int n, const auto &d) { const auto &[I, eig] = d; return n+Sym->mult(I)*eig.getnrkept(); });
  }
  void report() {
    std::cout << nrgdump4(nrkept, nrkeptmult, nrall, nrallmult) << std::endl;
  }
};

struct NotEnough : public std::exception {};

// Compute the number of states to keep in each subspace. Returns true if an insufficient number of states has been
// obtained in the diagonalization and we need to compute more states.
template<typename S>
void truncate_prepare(const Step &step, DiagInfo_tmpl<S> &diag, std::shared_ptr<Symmetry> Sym, const Params &P) {
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
template<typename S>
void calc_ZnD(const AllSteps_tmpl<S> &dm, Stats_tmpl<S> &stats, std::shared_ptr<Symmetry> Sym, const double T) {
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
  std::cout << "ZZG=" << HIGHPREC(stats.ZZG) << std::endl;
  for (const auto N : dm.Nall()) {
    const double w  = pow(Sym->nr_combs(), int(dm.Nend - N - 1)) / stats.ZZG;
    stats.wnfactor[N] = w; // These ratios enter the terms for the spectral function.
    stats.wn[N] = w * mpf_get_d(stats.ZnDG[N]); // This is w_n defined after Eq. (8) in the WvD paper.
  }
  const auto sumwn = ranges::accumulate(stats.wn, 0.0);
  std::cout << "sumwn=" << sumwn << " sumwn-1=" << sumwn - 1.0 << std::endl;
  my_assert(num_equal(sumwn, 1.0));  // Check the sum-rule.
}

template<typename S>
void report_ZnD(Stats_tmpl<S> &stats, const Params &P) {
  for (const auto N : P.Nall())
    std::cout << "ZG[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDG[N])) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "ZN[" << N << "]=" << HIGHPREC(mpf_get_d(stats.ZnDN[N])) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "w[" << N << "]=" << HIGHPREC(stats.wn[N]) << std::endl;
  for (const auto N : P.Nall())
    std::cout << "wfactor[" << N << "]=" << HIGHPREC(stats.wnfactor[N]) << std::endl;
}

// TO DO: use Boost.Multiprecision instead of low-level GMP calls
// https://www.boost.org/doc/libs/1_72_0/libs/multiprecision/doc/html/index.html
template<typename S>
void fdm_thermodynamics(const AllSteps_tmpl<S> &dm, Stats_tmpl<S> &stats, std::shared_ptr<Symmetry> Sym, const double T)
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
  std::cout << std::endl;
  std::cout << "Z_fdm=" << HIGHPREC(stats.Z_fdm) << std::endl;
  std::cout << "F_fdm=" << HIGHPREC(stats.F_fdm) << std::endl;
  std::cout << "E_fdm=" << HIGHPREC(stats.E_fdm) << std::endl;
  std::cout << "C_fdm=" << HIGHPREC(stats.C_fdm) << std::endl;
  std::cout << "S_fdm=" << HIGHPREC(stats.S_fdm) << std::endl;
  std::cout << std::endl;
  stats.td_fdm.save_values();
}

template<typename F, typename S>
auto trace(F fnc, const double rescale_factor, const DiagInfo_tmpl<S> &diag, std::shared_ptr<Symmetry> Sym) {
    auto b = 0.0;
    for (const auto &[I, eig] : diag)
      b += Sym->mult(I) * ranges::accumulate(eig.value_zero, 0.0, [fnc, rescale_factor](auto acc, const auto x) { 
        const auto betaE = rescale_factor * x; return acc + fnc(betaE) * exp(-betaE); });
    return b;
  }

// We calculate thermodynamic quantities before truncation to make better use of the available states. Here we
// compute quantities which are defined for all symmetry types. Other calculations are performed by calculate_TD
// member functions defined in symmetry.cc.
template<typename S>
void calculate_TD(const Step &step, const DiagInfo_tmpl<S> &diag, Stats_tmpl<S> &stats, Output_tmpl<S> &output, 
                  std::shared_ptr<Symmetry> Sym, const double additional_factor = 1.0) {
  // Rescale factor for energies. The energies are expressed in units of omega_N, thus we need to appropriately
  // rescale them to calculate the Boltzmann weights at the temperature scale Teff (Teff=scale/betabar).
  const auto rescale_factor = step.TD_factor() * additional_factor;
  const auto Z = trace([](double x) { return 1; }, rescale_factor, diag, Sym); // partition function
  const auto E = trace([](double x) { return x; }, rescale_factor, diag, Sym); // Tr[beta H]
  const auto E2 = trace([](double x) { return pow(x,2); }, rescale_factor, diag, Sym); // Tr[(beta H)^2]
  stats.Z = Z;
  stats.td.T  = step.Teff();
  stats.td.E  = E/Z;               // beta <H>
  stats.td.E2 = E2/Z;              // beta^2 <H^2>
  stats.td.C  = E2/Z - pow(E/Z,2); // C/k_B=beta^2(<H^2>-<H>^2)
  stats.td.F  = -log(Z);           // F/(k_B T)=-ln(Z)
  stats.td.S  = E/Z+log(Z);        // S/k_B=beta<H>+ln(Z)
  Sym->calculate_TD(step, diag, stats, rescale_factor);  // symmetry-specific calculation routine
  stats.td.save_values();
}

template<typename S>
void calculate_spectral_and_expv(const Step &step, Stats_tmpl<S> &stats, Output_tmpl<S> &output, Oprecalc_tmpl<S> &oprecalc, 
                                 const DiagInfo_tmpl<S> &diag, const IterInfo_tmpl<S> &iterinfo, const AllSteps_tmpl<S> &dm, 
                                 std::shared_ptr<Symmetry> Sym, const Params &P) {
  // Zft is used in the spectral function calculations using the conventional approach. We calculate it here, in
  // order to avoid recalculations later on.
  stats.Zft = grand_canonical_Z(step, diag, Sym);
  if (string(P.specgt) != "" || string(P.speci1t) != "" || string(P.speci2t) != "")
    stats.Zgt = grand_canonical_Z(step, diag, Sym, 1.0/(P.gtp*step.scT()) ); // exp(-x*gtp)
  if (string(P.specchit) != "") 
    stats.Zchit = grand_canonical_Z(step, diag, Sym, 1.0/(P.chitp*step.scT()) ); // exp(-x*chitp)
  DensMatElements_tmpl<S> rho, rhoFDM;
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
template<typename S>
void perform_basic_measurements(const Step &step, const DiagInfo_tmpl<S> &diag, std::shared_ptr<Symmetry> Sym, 
                                Stats_tmpl<S> &stats, Output_tmpl<S> &output) {
  output.dump_all_energies(diag, step.ndx());
  calculate_TD(step, diag, stats, output, Sym);
  output.annotated.dump(step, diag, stats, Sym);
}

// Subspaces for the new iteration
template<typename S>
auto new_subspaces(const DiagInfo_tmpl<S> &diagprev, std::shared_ptr<Symmetry> Sym) {
  std::set<Invar> subspaces;
  for (const auto &I : diagprev.subspaces()) {
    const auto all = Sym->new_subspaces(I);
    const auto non_empty = all | ranges::views::filter([&Sym](const auto &In) { return Sym->Invar_allowed(In); }) | ranges::to<std::vector>();
    std::copy(non_empty.begin(), non_empty.end(), std::inserter(subspaces, subspaces.end()));
  }
  return subspaces;
}

template<typename S>
typename traits<S>::Matrix prepare_task_for_diag(const Step &step, const Invar &I, const Opch_tmpl<S> &opch, const Coef_tmpl<S> &coef, 
                                                 const DiagInfo_tmpl<S> &diagprev, std::shared_ptr<Symmetry> Sym, const Params &P) {
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

template<typename S>
auto diagonalisations_OpenMP(const Step &step, const Opch_tmpl<S> &opch, const Coef_tmpl<S> &coef, const DiagInfo_tmpl<S> &diagprev,
                             const std::vector<Invar> &tasks, const DiagParams &DP, std::shared_ptr<Symmetry> Sym, const Params &P) {
  DiagInfo_tmpl<S> diagnew;
  const auto nr = tasks.size();
  size_t itask = 0;
  // cppcheck-suppress unreadVariable symbolName=nth
  const int nth = P.diagth; // NOLINT
#pragma omp parallel for schedule(dynamic) num_threads(nth)
  for (itask = 0; itask < nr; itask++) {
    const Invar I  = tasks[itask];
    auto h = prepare_task_for_diag(step, I, opch, coef, diagprev, Sym, P); // non-const, consumed by diagonalise()
    const int thid = omp_get_thread_num();
#pragma omp critical
    { nrglog('(', "Diagonalizing " << I << " size=" << h.size1() << " (task " << itask + 1 << "/" << nr << ", thread " << thid << ")"); }
    Eigen e = diagonalise<S>(h, DP);
#pragma omp critical
    { diagnew[I] = e; }
  }
  return diagnew;
}

#ifdef NRG_MPI
enum TAG : int { TAG_EXIT = 1, TAG_DIAG_DBL, TAG_DIAG_CMPL, TAG_SYNC, TAG_MATRIX, TAG_INVAR, 
                 TAG_MATRIX_SIZE, TAG_MATRIX_LINE, TAG_EIGEN_INT, TAG_EIGEN_VEC };

void mpi_send_params(const DiagParams &DP) {
  mpilog("Sending diag parameters " << DP.diag << " " << DP.diagratio);
  for (auto i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_SYNC, 0);
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
    std::cout << "MPI communication error. rank=" << mpiw->rank() << std::endl;
    mpienv->abort(1);
  }
}

// NOTE: MPI is limited to message size of 2GB (or 4GB). For big problems we thus need to send objects line by line.
template<typename S> void mpi_send_matrix(const int dest, const ublas::matrix<S> &m) {
  const auto size1 = m.size1();
  mpiw->send(dest, TAG_MATRIX_SIZE, size1);
  const auto size2 = m.size2();
  mpiw->send(dest, TAG_MATRIX_SIZE, size2);
  mpilog("Sending matrix of size " << size1 << " x " << size2 << " line by line to " << dest);
  for (const auto i: range0(size1)) {
    ublas::vector<typename traits<S>::t_matel> vec = ublas::matrix_row<const ublas::matrix<S>>(m, i); // YYY
    mpiw->send(dest, TAG_MATRIX_LINE, vec);
  }
}

template<typename S> auto mpi_receive_matrix(const int source) {
  size_t size1;
  check_status(mpiw->recv(source, TAG_MATRIX_SIZE, size1));
  size_t size2;
  check_status(mpiw->recv(source, TAG_MATRIX_SIZE, size2));
  typename traits<S>::Matrix m(size1, size2);
  mpilog("Receiving matrix of size " << size1 << " x " << size2 << " line by line from " << source);
  for (const auto i: range0(size1)) {
    ublas::vector<typename traits<S>::t_matel> vec;
    check_status(mpiw->recv(source, TAG_MATRIX_LINE, vec));
    my_assert(vec.size() == size2);
    ublas::matrix_row<typename traits<S>::Matrix>(m, i) = vec;
  }
  return m;
}

template<typename S> void mpi_send_eigen(const int dest, const Eigen_tmpl<S> &eig) {
  mpilog("Sending eigen from " << mpiw->rank() << " to " << dest);
  mpiw->send(dest, TAG_EIGEN_VEC, eig.value_orig);
  mpi_send_matrix<S>(dest, eig.matrix);
}

template<typename S> auto mpi_receive_eigen(const int source) {
  mpilog("Receiving eigen from " << source << " on " << mpiw->rank());
  Eigen_tmpl<S> eig;
  check_status(mpiw->recv(source, TAG_EIGEN_VEC, eig.value_orig));
  eig.matrix = mpi_receive_matrix<S>(source);
  return eig;
}

// Read results from a slave process.
template<typename S> std::pair<Invar, Eigen_tmpl<S>> read_from(int source) {
  mpilog("Reading results from " << source);
  const auto eig = mpi_receive_eigen<S>(source);
  Invar Irecv;
  check_status(mpiw->recv(source, TAG_INVAR, Irecv));
  mpilog("Received results for subspace " << Irecv << " [nr=" << eig.getnrstored() << ", dim=" << eig.getdim() << "]");
  my_assert(eig.value_orig.size() == eig.matrix.size1());
  my_assert(eig.matrix.size1() <= eig.matrix.size2());
  return {Irecv, eig};
}

template<typename S>
DiagInfo diagonalisations_MPI(const Step &step, const Opch_tmpl<S> &opch, const Coef_tmpl<S> &coef, const DiagInfo_tmpl<S> &diagprev, 
                              const std::vector<Invar> &tasks, const DiagParams &DP, std::shared_ptr<Symmetry> Sym, const Params &P) {
  DiagInfo_tmpl<S> diagnew;
  mpi_send_params(DP); // Synchronise parameters
  std::list<Invar> todo; // List of all the tasks to handle
  std::copy(tasks.begin(), tasks.end(), std::back_inserter(todo)); // BBB: constr
  std::list<Invar> done; // List of finished tasks.
  // List of the available computation nodes (including the master,
  // which is always at the very beginnig of the deque).
  std::deque<int> nodes;
  for (const auto i: range0(mpiw->size())) nodes.push_back(i); // BBB: iota
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
    auto h = prepare_task_for_diag(step, I, opch, coef, diagprev, Sym, P); // non-const
    nrglog('M', "Scheduler: job " << I << " (dim=" << h.size1() << ")" << " on node " << i);
    if (i == 0) {
      // On master, diagonalize immediately.
      diagnew[I] = diagonalise<S>(h, DP);
      nodes.push_back(0);
      done.push_back(I);
    } else {
      mpiw->send(i, std::is_same_v<S, double> ? TAG_DIAG_DBL : TAG_DIAG_CMPL, 0);
      mpi_send_matrix<S>(i, h);
      mpiw->send(i, TAG_INVAR, I);
    }
    // Check for terminated jobs
    while (auto status = mpiw->iprobe(mpi::any_source, TAG_EIGEN_VEC)) {
      nrglog('M', "Receiveing results from " << status->source());
      const auto [Irecv, eig] = read_from<S>(status->source());
      diagnew[Irecv] = eig;
      done.push_back(Irecv);
      // The node is now available for new tasks!
      nodes.push_back(status->source());
    }
  }
  // Keep reading results sent from the slave processes until all tasks have been completed.
  while (done.size() != tasks.size()) {
    const auto status = mpiw->probe(mpi::any_source, TAG_EIGEN_VEC);
    const auto [Irecv, eig]  = read_from<S>(status.source());
    diagnew[Irecv] = eig;
    done.push_back(Irecv);
  }
  return diagnew;
}

// Handle a diagonalisation request
template<typename S> void slave_diag(const int master, const DiagParams &DP) {
  // 1. receive the matrix and the subspace identification
  auto m = mpi_receive_matrix<S>(master);
  Invar I;
  check_status(mpiw->recv(master, TAG_INVAR, I));
  // 2. preform the diagonalisation
  const auto eig = diagonalise<S>(m, DP);
  // 3. send back the results
  mpi_send_eigen<S>(master, eig);
  mpiw->send(master, TAG_INVAR, I);
}
#endif

// Build matrix H(ri;r'i') in each subspace and diagonalize it
template<typename S>
auto diagonalisations(const Step &step, const Opch_tmpl<S> &opch, const Coef_tmpl<S> &coef, const DiagInfo_tmpl<S> &diagprev, 
                      const std::vector<Invar> &tasks, const double diagratio, std::shared_ptr<Symmetry> Sym, const Params &P) {
  TIME("diag");
#ifdef NRG_MPI
  return diagonalisations_MPI<scalar>(step, opch, coef, diagprev, tasks, DiagParams(P, diagratio), Sym, P);
#else
  return diagonalisations_OpenMP(step, opch, coef, diagprev, tasks, DiagParams(P, diagratio), Sym, P);
#endif
}

// Determine the structure of matrices in the new NRG shell
template<typename S>
QSrmax::QSrmax(const DiagInfo_tmpl<S> &diagprev, std::shared_ptr<Symmetry> Sym) {
  for (const auto &I : new_subspaces(diagprev, Sym))
    (*this)[I] = Rmaxvals{I, Sym->ancestors(I), diagprev, Sym};
}

// Recalculate irreducible matrix elements for Wilson chains.
template<typename S>
void recalc_irreducible(const Step &step, const DiagInfo_tmpl<S> &diag, const QSrmax &qsrmax, Opch &opch, std::shared_ptr<Symmetry> Sym, const Params &P) {
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

template<typename S>
auto do_diag(const Step &step, IterInfo_tmpl<S> &iterinfo, const Coef_tmpl<S> &coef, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diagprev,
             QSrmax &qsrmax, std::shared_ptr<Symmetry> Sym, const Params &P) {
  step.infostring();
  Sym->show_coefficients(step, coef);
  auto tasks = qsrmax.task_list();
  double diagratio = P.diagratio; // non-const
  DiagInfo_tmpl<S> diag;
  while (true) {
    try {
      if (step.nrg()) {
        if (!(P.resume && int(step.ndx()) <= P.laststored))
          diag = diagonalisations(step, iterinfo.opch, coef, diagprev, tasks, diagratio, Sym, P); // compute in first run
        else
          diag = DiagInfo_tmpl<S>(step.ndx(), false); // or read from disk
      }
      if (step.dmnrg()) {
        diag = DiagInfo_tmpl<S>(step.ndx(), P.removefiles); // read from disk in second run
        diag.subtract_GS_energy(stats.GS_energy);
      }
      stats.Egs = diag.find_groundstate();
      if (step.nrg()) // should be done only once!
        diag.subtract_Egs(stats.Egs);
      const auto cluster_mapping = find_clusters(diag.sorted_energies(), P.fixeps);
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
template<typename S>
void calc_abs_energies(const Step &step, DiagInfo_tmpl<S> &diag, const Stats_tmpl<S> &stats) {
  for (auto &eig : diag.eigs()) {
    eig.absenergyN = eig.value_zero * step.scale();        // referenced to the lowest energy in current NRG step (not modified later on)
    eig.absenergy = eig.absenergyN;
    for (auto &x : eig.absenergy) x += stats.total_energy; // absolute energies (not modified later on)
    eig.absenergyG = eig.absenergy;                        // referenced to the absolute 0 (updated by shft_abs_energies())
  }
}

// Perform processing after a successful NRG step. Also called from doZBW() as a final step.
template<typename S>
void after_diag(const Step &step, IterInfo_tmpl<S> &iterinfo, Stats_tmpl<S> &stats, DiagInfo_tmpl<S> &diag, Output_tmpl<S> &output,
                QSrmax &qsrmax, AllSteps_tmpl<S> &dm, Oprecalc_tmpl<S> &oprecalc, std::shared_ptr<Symmetry> Sym, const Params &P) {
  stats.total_energy += stats.Egs * step.scale(); // stats.Egs has already been initialized
  std::cout << "Total energy=" << HIGHPREC(stats.total_energy) << "  Egs=" << HIGHPREC(stats.Egs) << std::endl;
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
template<typename S>
auto iterate(const Step &step, IterInfo_tmpl<S> &iterinfo, const Coef_tmpl<S> &coef, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diagprev,
             Output_tmpl<S> &output, AllSteps_tmpl<S> &dm, Oprecalc &oprecalc, std::shared_ptr<Symmetry> Sym, const Params &P) {
  QSrmax qsrmax{diagprev, Sym};
  auto diag = do_diag(step, iterinfo, coef, stats, diagprev, qsrmax, Sym, P);
  after_diag(step, iterinfo, stats, diag, output, qsrmax, dm, oprecalc, Sym, P);
  iterinfo.trim_matrices(diag);
  diag.clear_eigenvectors();
  time_mem::memory_time_brief_report();
  return diag;
}

// Perform calculations with quantities from 'data' file
template<typename S>
void docalc0(Step &step, const IterInfo_tmpl<S> &iterinfo, const DiagInfo_tmpl<S> &diag0, Stats_tmpl<S> &stats, Output_tmpl<S> &output, 
             Oprecalc_tmpl<S> &oprecalc, std::shared_ptr<Symmetry> Sym, const Params &P) {
  step.set(P.Ninit - 1); // in the usual case with Ninit=0, this will result in N=-1
  std::cout << endl << "Before NRG iteration";
  std::cout << " (N=" << step.N() << ")" << std::endl;
  perform_basic_measurements(step, diag0, Sym, stats, output);
  AllSteps_tmpl<S> empty_dm(0, 0);
  calculate_spectral_and_expv(step, stats, output, oprecalc, diag0, iterinfo, empty_dm, Sym, P);
  if (P.checksumrules) operator_sumrules(iterinfo, Sym);
}

// doZBW() takes the place of iterate() called from main_loop() in the case of zero-bandwidth calculation.
// It replaces do_diag() and calls after_diag() as the last step.
template<typename S>
auto nrg_ZBW(Step &step, IterInfo_tmpl<S> &iterinfo, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diag0, Output_tmpl<S> &output, 
             AllSteps_tmpl<S> &dm, Oprecalc_tmpl<S> &oprecalc, std::shared_ptr<Symmetry> Sym, const Params &P) {
  std::cout << std::endl << "Zero bandwidth calculation" << std::endl;
  step.set_ZBW();
  // --- begin do_diag() equivalent
  DiagInfo_tmpl<S> diag;
  if (step.nrg())
    diag = diag0;
  if (step.dmnrg()) {
    diag = DiagInfo_tmpl<S>(step.ndx(), P.removefiles);
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

template<typename S>
auto nrg_loop(Step &step, IterInfo_tmpl<S> &iterinfo, const Coef_tmpl<S> &coef, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diag0,
              Output &output, AllSteps_tmpl<S> &dm, Oprecalc_tmpl<S> &oprecalc, std::shared_ptr<Symmetry> Sym, const Params &P) {
  auto diag = diag0;
  for (step.init(); !step.end(); step.next())
    diag = iterate(step, iterinfo, coef, stats, diag, output, dm, oprecalc, Sym, P);
  step.set(step.lastndx());
  return diag;
}

template<typename S>
auto run_nrg(Step &step, IterInfo_tmpl<S> &iterinfo, const Coef_tmpl<S> &coef, Stats_tmpl<S> &stats, const DiagInfo_tmpl<S> &diag0,
             AllSteps_tmpl<S> &dm, std::shared_ptr<Symmetry> Sym, const Params &P) {
  diag0.states_report(Sym->multfnc());
  auto oprecalc = Oprecalc(step.runtype, iterinfo, Sym, P);
  auto output = Output(step.runtype, iterinfo, stats, P);
  // If calc0=true, a calculation of TD quantities is performed before starting the NRG iteration.
  if (step.nrg() && P.calc0 && !P.ZBW)
    docalc0(step, iterinfo, diag0, stats, output, oprecalc, Sym, P);
  auto diag = P.ZBW ? nrg_ZBW(step, iterinfo, stats, diag0, output, dm, oprecalc, Sym, P) 
    : nrg_loop(step, iterinfo, coef, stats, diag0, output, dm, oprecalc, Sym, P);
  fmt::print(fmt::emphasis::bold | fg(fmt::color::red), FMT_STRING("\nTotal energy: {:.18}\n"), stats.total_energy);
  stats.GS_energy = stats.total_energy;
  if (step.nrg() && P.dumpsubspaces) dm.dump_subspaces();
  fmt::print("\n** Iteration completed.\n\n");
  return diag;
}

void print_about_message() {
  fmt::print(fmt::emphasis::bold, "NRG Ljubljana - (c) rok.zitko@ijs.si\n");
  fmt::print(fmt::emphasis::bold, "Timestamp: {}\n",  __TIMESTAMP__);
  fmt::print(fmt::emphasis::bold, "Compiled on {} at {}\n\n", __DATE__, __TIME__);
}

template<typename S>
std::unique_ptr<Symmetry> get(const std::string &sym_string, const Params &P, Allfields &allfields)
{
  if (sym_string == "QS")     return std::make_unique<SymmetryQS_tmpl<S>>(P, allfields);
  if (sym_string == "QSZ")    return std::make_unique<SymmetryQSZ_tmpl<S>>(P, allfields);
#ifdef NRG_SYM_MORE
  if (sym_string == "ISO")    return std::make_unique<SymmetryISO_tmpl<S>>(P, allfields);
  if (sym_string == "ISO2")   return std::make_unique<SymmetryISO2_tmpl<S>>(P, allfields);
  if (sym_string == "ISOSZ")  return std::make_unique<SymmetryISOSZ_tmpl<S>>(P, allfields);
  if (sym_string == "SPSU2")  return std::make_unique<SymmetrySPSU2_tmpl<S>>(P, allfields);
  if (sym_string == "SPU1")   return std::make_unique<SymmetrySPU1_tmpl<S>>(P, allfields);
#endif
#ifdef NRG_SYM_ALL
  if (sym_string == "DBLISOSZ")  return std::make_unique<SymmetryDBLISOSZ_tmpl<S>>(P, allfields);
  if (sym_string == "DBLSU2")    return std::make_unique<SymmetryDBLSU2_tmpl<S>>(P, allfields);
  if (sym_string == "ISOLR")     return std::make_unique<SymmetryISOLR_tmpl<S>>(P, allfields);
  if (sym_string == "ISO2LR")    return std::make_unique<SymmetryISO2LR_tmpl<S>>(P, allfields);
  if (sym_string == "ISOSZLR")   return std::make_unique<SymmetryISOSZLR_tmpl<S>>(P, allfields);
  if (sym_string == "NONE")      return std::make_unique<SymmetryNONE_tmpl<S>>(P, allfields);
  if (sym_string == "P")         return std::make_unique<SymmetryP_tmpl<S>>(P, allfields);
  if (sym_string == "PP")        return std::make_unique<SymmetryPP_tmpl<S>>(P, allfields);
  if (sym_string == "QJ")        return std::make_unique<SymmetryQJ_tmpl<S>>(P, allfields);
  if (sym_string == "QSLR")      return std::make_unique<SymmetryQSLR_tmpl<S>>(P, allfields); 
  if (sym_string == "QST")       return std::make_unique<SymmetryQST_tmpl<S>>(P, allfields);
  if (sym_string == "QSTZ")      return std::make_unique<SymmetryQSTZ_tmpl<S>>(P, allfields);
  if (sym_string == "QSZLR")     return std::make_unique<SymmetryQSZLR_tmpl<S>>(P, allfields);
  if (sym_string == "QSZTZ")     return std::make_unique<SymmetryQSZTZ_tmpl<S>>(P, allfields);
  if (sym_string == "SL")        return std::make_unique<SymmetrySL_tmpl<S>>(P, allfields);
  if (sym_string == "SL3")       return std::make_unique<SymmetrySL3_tmpl<S>>(P, allfields);
  if (sym_string == "SPSU2LR")   return std::make_unique<SymmetrySPSU2LR_tmpl<S>>(P, allfields);
  if (sym_string == "SPSU2T")    return std::make_unique<SymmetrySPSU2T_tmpl<S>>(P, allfields);
  if (sym_string == "SPU1LR")    return std::make_unique<SymmetrySPU1LR_tmpl<S>>(P, allfields);
  if (sym_string == "SU2")       return std::make_unique<SymmetrySU2_tmpl<S>>(P, allfields);
  if (sym_string == "U1")        return std::make_unique<SymmetryU1_tmpl<S>>(P, allfields);
 #ifdef NRG_COMPLEX
  if (sym_string == "QSC3")      return std::make_unique<SymmetryQSC3_tmpl<S>>(P, allfields);
  if (sym_string == "SPSU2C3")   return std::make_unique<SymmetrySPSU2C3_tmpl<S>>(P, allfields);
 #endif
#endif 
  throw std::runtime_error("Unknown symmetry " + sym_string);
}

// Called immediately after parsing the information about the number of channels from the data file. This ensures
// that Invar can be parsed correctly.
template <typename S>
std::shared_ptr<Symmetry> set_symmetry(const Params &P, Stats_tmpl<S> &stats) {
  my_assert(P.channels > 0 && P.combs > 0); // must be set at this point
  std::cout << "SYMMETRY TYPE: " << P.symtype.value() << std::endl;
  auto Sym = get<S>(P.symtype.value(), P, stats.td.allfields);
  Sym->load();
  Sym->erase_first();
  return Sym;
}

template <typename S> class NRG_calculation {
private:
  // XXX: Workdir workdir;
  Params P;
  Stats_tmpl<S> stats;
public:
  NRG_calculation() : P("param", "param", workdir), stats(P) {}
  void go() {
    auto [diag0, iterinfo, coef, Sym] = read_data<S>(P, stats);
    Step step{P, RUNTYPE::NRG};
    AllSteps_tmpl<S> dm(P.Ninit, P.Nlen);
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
      if (std::string(P.stopafter) == "rho") exit1("*** Stopped after the DM calculation.");
      auto [diag0_dm, iterinfo_dm, coef_dm, Sym_dm] = read_data<S>(P, stats);
      Step step_dmnrg{P, RUNTYPE::DMNRG};
      run_nrg(step_dmnrg, iterinfo_dm, coef_dm, stats, diag0_dm, dm, Sym_dm, P);
      my_assert(num_equal(stats.GS_energy, stats.total_energy));
    }
  }
  ~NRG_calculation() {
    if (P.done) { std::ofstream D("DONE"); } // Indicate completion by creating a flag file
  }
};
  
// Master process does most of the i/o and passes calculations to the slaves.
void run_nrg_master() {
#ifdef NRG_REAL
  NRG_calculation<double> calc;
#endif
#ifdef NRG_COMPLEX
  NRG_calculation<std::complex<double>> calc;
#endif
  calc.go();
#ifdef NRG_MPI
  for (auto i = 1; i < mpiw->size(); i++) mpiw->send(i, TAG_EXIT, 0);
#endif
}

void run_nrg_slave() {
#ifdef NRG_MPI
  constexpr auto master = 0;
  DiagParams DP;
  for (;;) {
    if (mpiw->iprobe(master, mpi::any_tag)) { // message can be received.
      int task;
      const auto status = mpiw->recv(master, mpi::any_tag, task);
      mpilog("Slave " << mpiw->rank() << " received message with tag " << status.tag());
      check_status(status);
      switch (status.tag()) {
        case TAG_SYNC:
          DP = mpi_receive_params();
          break;
        case TAG_DIAG_DBL:
          slave_diag<double>(master, DP);
          break;
        case TAG_DIAG_CMPL:
          slave_diag<std::complex<double>>(master, DP);
          break;      
        case TAG_EXIT:
          return; // exit from run_slave()
        default: 
          std::cout << "MPI error: unknown tag on " << mpiw->rank() << std::endl; 
          break;
      }
    } else usleep(100); // sleep to reduce the load on the computer. (OpenMPI "feature" workaround)
  }
#endif
}
