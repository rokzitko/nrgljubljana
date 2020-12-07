#ifndef _eigen_hpp_
#define _eigen_hpp_

#include <vector>
#include <string>
#include <limits> // quiet_NaN

#include <boost/range/adaptor/map.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>

#include <range/v3/all.hpp>

#include "basicio.hpp"
#include "portabil.hpp"
#include "traits.hpp"
#include "invar.hpp"
#include "params.hpp"
#include "numerics.hpp"
#include "h5.hpp"

namespace NRG {

template<scalar S, typename t_eigen = eigen_traits<S>>
class Values {
  private:
   std::vector<t_eigen> v;
   double scale = std::numeric_limits<double>::quiet_NaN();
   double shift = std::numeric_limits<double>::quiet_NaN();
   double GS_energy = std::numeric_limits<double>::quiet_NaN();
  public:
 //  Values() = default;
 //  explicit Values(const double scale) : scale(scale) {}
   void resize(const size_t size) { v.resize(size); } // XXX for testing purposes
   auto rel(const size_t i) const { return v[i]; }
   auto abs(const size_t i) const { my_assert(std::isfinite(scale)); return rel(i) * scale; }
   auto rel_zero(const size_t i) const { my_assert(std::isfinite(shift)); return rel(i)-shift; }
   auto abs_zero(const size_t i) const { return (rel(i)-shift) * scale; }
   auto absG(const size_t i) const { return rel(i)*scale - GS_energy; }
   auto size() const { return v.size(); }
   auto lowest_rel() const { return v.front(); }
   auto all_rel() const { return v; }
   void set_scale(const double scale_) { scale = scale_; }
   void set_shift(const double shift_) { shift = shift_; }
   void set_GS_energy(const double GS_energy_) { GS_energy = GS_energy_; }
   void save(boost::archive::binary_oarchive &oa) const {
     oa << v << scale << shift << GS_energy;
   }
   void load(boost::archive::binary_iarchive &ia) {
     ia >> v >> scale >> shift >> GS_energy;
   }
   void copy(const std::vector<t_eigen> &in, const size_t M) {
     v.resize(M);
     std::copy(in.begin(), in.begin() + M, v.begin());
   }
   void copy(const std::vector<t_eigen> &in) { 
     v = in;
   }
  void move(const std::vector<t_eigen> &&in) { 
     v = std::move(in);
  }
  auto begin() {
    return v.begin();
  }
  auto end() {
    return v.end();
  }
};

// Result of a diagonalisation: eigenvalues and eigenvectors
template <scalar S, typename RVector = RVector_traits<S>, typename Matrix = Matrix_traits<S>> 
struct RawEigen {
public:
  RVector val;
  Matrix vec;
  RawEigen() = default;
  RawEigen(const size_t M, const size_t dim) {
    my_assert(M <= dim);
    val.resize(M);
    vec.resize(M, dim);
  }
  auto getnrcomputed() const { my_assert(val.size() == vec.size1()); return val.size(); } // nr eigenvalue/eigenvector pairs
  auto getdim() const { return vec.size2(); } // matrix dimension (length of eigenvectors)
  void dump_eigenvalues(const size_t max_nr = std::numeric_limits<size_t>::max(), std::ostream &F = std::cout) const {
    F << "eig= " << std::setprecision(std::numeric_limits<double>::max_digits10);
    ranges::for_each_n(val.begin(), std::min(val.size(), max_nr), [&F](const double x) { F << x << ' '; });
    F << std::endl;
  }
  void check_diag(const double NORMALIZATION_EPSILON = 1e-12, const double ORTHOGONALITY_EPSILON = 1e-12) const { // XXX: throw exception instead?
    const auto M   = getnrcomputed();
    const auto dim = getdim();
    // Check normalization
    for (const auto r : range0(M)) {
      assert_isfinite(val[r]);
      S sumabs{};
      for (const auto j : range0(dim)) sumabs += conj_me(vec(r, j)) * vec(r, j);
      my_assert(num_equal(abs(sumabs), 1.0, NORMALIZATION_EPSILON));
    }
    // Check orthogonality
    for (const auto r1 : range0(M))
      for (const auto r2 : boost::irange(r1 + 1, M)) {
        S skpdt{};
        for (const auto j : range0(dim)) skpdt += conj_me(vec(r1, j)) * vec(r2, j);
        my_assert(num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON));
      }
  }};

// High-level representation of eigenvalues/eigenvectors
template <scalar S, typename EVEC = evec_traits<S>, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>> 
class Eigen {
public:
  Values<S> values; // eigenvalues
  EVEC value_zero;  // eigenvalues with Egs subtracted
  EVEC value_corr;  // eigenvalues corrected for floating-point round-off errors
  Matrix matrix;    // eigenvectors
  Eigen() = default;
  explicit Eigen(const size_t M, const size_t dim) { // XXX for testing only
    my_assert(M <= dim);
    values.resize(M);
    value_zero.resize(M);
    matrix.resize(M, dim);
  }
  explicit Eigen(RawEigen<S> && raw) {
    values.move(std::move(raw.val));
    matrix = std::move(raw.vec);
  }
  [[nodiscard]] auto getnrcomputed() const { return values.size(); }     // number of computed eigenpairs
  [[nodiscard]] auto getdim() const { return matrix.size2(); }           // valid also after the split_in_blocks_Eigen() call
 private:
  long nrpost = -1;   // number of eigenpairs after truncation (-1: keep all)
  long nrstored = -1; // number of eigenpairs currently held in store 
  [[nodiscard]] auto getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }     // number of states after truncation
 public:
  [[nodiscard]] auto getnrall() const { return getnrcomputed(); }                              // all = all computed
  [[nodiscard]] auto getnrkept() const { return getnrpost(); }                                 // # of kept states
  [[nodiscard]] auto getnrdiscarded() const { return getnrcomputed()-getnrpost(); }            // # of discarded states
  [[nodiscard]] auto getnrstored() const  { return nrstored == -1 ? value_zero.size() : nrstored; }  // number of stored states
  [[nodiscard]] auto all() const { return range0(getnrcomputed()); }                           // iterator over all states
  [[nodiscard]] auto kept() const { return range0(getnrpost()); }                              // iterator over kept states
  [[nodiscard]] auto discarded() const { return boost::irange(getnrpost(), getnrcomputed()); } // iterator over discarded states
  [[nodiscard]] auto stored() const { return range0(getnrstored()); }                          // iterator over all stored states
  auto value_zero_kept() const { return ranges::subrange(value_zero.begin(), value_zero.begin() + getnrkept()); }
  // NOTE: "absolute" energy means that it is expressed in the absolute energy scale rather than SCALE(N).
  EVEC absenergy;      // absolute energies
  EVEC absenergyG;     // absolute energies (0 is the absolute ground state of the system) [SAVED TO FILE]
  EVEC absenergy_zero; // absolute energies (referenced to the lowest energy in the N-th step)
  // 'blocks' contains eigenvectors separated according to the invariant subspace from which they originate.
  // Required for using efficient BLAS routines when performing recalculations of the matrix elements.
  std::vector<Matrix> blocks;
  const Matrix & Ublock(const size_t i) const { // 1-based MMA index, called from recalc_f()
    my_assert(1 <= i && i <= blocks.size());
    return blocks[i-1];
  }
  // Truncate to nrpost states.
  void truncate_prepare(const size_t nrpost_) {
    nrpost = nrpost_;
    my_assert(nrpost <= getnrstored());
  }
  void truncate_perform() {
    for (auto &i : blocks) {
      my_assert(nrpost <= i.size1());
      i.resize(nrpost, i.size2());
    }
    nrstored = nrpost;
    value_zero.resize(nrpost); // ZZZ: necessary?? YES!
  }
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors form an identity matrix. This is used to
  // represent the spectral decomposition in the eigenbasis itself. Called when building DiagInfo from 'data' file.
  void diagonal(const EVEC &v) {
    value_zero = v; // YYY
    values.copy(v);
    values.set_shift(0.0);
    matrix   = ublas::identity_matrix<t_eigen>(v.size());
  }
  void subtract_Egs(const t_eigen Egs) {
    value_zero = values.all_rel(); // XXX
    for (auto &x : value_zero) x -= Egs; // XXX: subtract a scalar [fix after moving to Eigen]
    my_assert(value_zero[0] >= 0); // XXX
    values.set_shift(Egs);
    my_assert(values.rel_zero(0) == value_zero[0]); // XXX
  }
  void subtract_GS_energy(const t_eigen GS_energy) {
    for (auto &x : absenergyG) x -= GS_energy; // XXX
    my_assert(absenergyG[0] >= 0); // XXX
    values.set_GS_energy(GS_energy);
  }
  auto diagonal_exp(const double factor) const { // produce a diagonal matrix with exp(-factor*E) diagonal elements
    const auto dim = getnrstored();
    auto m = Zero_matrix<S>(dim);
    for (const auto i: range0(dim)) 
      m(i, i) = exp(-value_zero[i] * factor);
    return m;
  }
  template<typename F> auto trace(F fnc, const double factor) const { // Tr[fnc(factor*E) exp(-factor*E)]
    return ranges::accumulate(value_zero, 0.0, {}, [fnc, factor](const auto x) { return fnc(factor*x) * exp(-factor*x); });
  }
  void save(boost::archive::binary_oarchive &oa) const {
    values.save(oa);
    NRG::save(oa, matrix);
    oa << value_zero << nrpost << nrstored << absenergy << absenergyG << absenergy_zero;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    values.load(ia);
    NRG::load(ia, matrix);
    ia >> value_zero >> nrpost >> nrstored >> absenergy >> absenergyG >> absenergy_zero;
  }
  void h5save(H5Easy::File &fd, const std::string &name, const bool write_absG) const {
    H5Easy::dump(fd, name + "/value_orig",     values.all_rel());
    H5Easy::dump(fd, name + "/value_zero",     value_zero);
    H5Easy::dump(fd, name + "/absenergy",      absenergy);
    H5Easy::dump(fd, name + "/absenergy_zero", absenergy_zero);
    if (write_absG) 
      H5Easy::dump(fd, name + "/absenergyG",   absenergyG);
    h5_dump_matrix(fd, name + "/matrix", matrix);
    h5_dump_scalar(fd, name + "/nrkept", getnrkept());
    h5_dump_scalar(fd, name + "/nrstored", getnrstored()); // XXX: do we need this?
  }
};

// Full information after diagonalizations (eigenspectra in all subspaces)
template <scalar S, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>> 
class DiagInfo : public std::map<Invar, Eigen<S>> {
 public:
   explicit DiagInfo() = default;
   DiagInfo(std::istream &fdata, const size_t nsubs, const Params &P) {
     for (const auto i : range1(nsubs)) {
       const auto I = read_one<Invar>(fdata);
       auto energies = read_std_vector<t_eigen>(fdata);
       if (!P.data_has_rescaled_energies && !P.absolute)
         for (auto &x : energies) x /= P.SCALE(P.Ninit); // rescale to the suitable energy scale
       (*this)[I].diagonal(energies);
     }
     my_assert(this->size() == nsubs);
   }
   [[nodiscard]] auto subspaces() const { return *this | boost::adaptors::map_keys; }
   [[nodiscard]] auto eigs() const { return *this | boost::adaptors::map_values; }
   [[nodiscard]] auto eigs() { return *this | boost::adaptors::map_values; }
   [[nodiscard]] auto find_groundstate() const {
     const auto [Iground, eig] = *ranges::min_element(*this, {}, [](const auto &a) { return a.second.values.lowest_rel(); });
     const auto Egs = eig.values.lowest_rel();
     return Egs;
   }
   void subtract_Egs(const t_eigen Egs) {
     ranges::for_each(eigs(), [Egs](auto &eig)       { eig.subtract_Egs(Egs); });
   }
   t_eigen Egs_subtraction() {
     const auto Egs = find_groundstate();
     subtract_Egs(Egs);
     return Egs;
   }
   void subtract_GS_energy(const t_eigen GS_energy) {
     ranges::for_each(eigs(), [GS_energy](auto &eig) { eig.subtract_GS_energy(GS_energy); });
   }
   std::vector<t_eigen> sorted_energies_rel_zero() const { // YYY
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs())
       energies.insert(energies.end(), eig.value_zero.begin(), eig.value_zero.end());
     return energies | ranges::move | ranges::actions::sort;
   }
   std::vector<t_eigen> sorted_energies_rel() const { // YYY
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs()) {
       const auto &all = eig.values.all_rel();
       energies.insert(energies.end(), all.begin(), all.end());
     }
     return energies | ranges::move | ranges::actions::sort;
   }
   void dump_value_zero(std::ostream &F) const {
     for (const auto &[I, eig]: *this)
       F << "Subspace: " << I << std::endl << eig.value_zero << std::endl;
   }
   void truncate_perform() {
     ranges::for_each(eigs(), [](auto &eig){ eig.truncate_perform(); }); // Truncate subspace to appropriate size
   }
   [[nodiscard]] auto size_subspace(const Invar &I) const {
     const auto f = this->find(I);
     return f != this->cend() ? f->second.getnrstored() : 0;
   }
   [[nodiscard]] auto subs(const Invar &I1, const Invar &I2) const {
     return std::make_pair(this->at(I1), this->at(I2));
   }
   [[nodiscard]] auto dims(const Invar &I1, const Invar &I2) const { // Determines matrix sizes for operators (# stored)
     return std::make_pair(size_subspace(I1), size_subspace(I2));
   }
   void clear_eigenvectors() {
     ranges::for_each(eigs(), [](auto &eig){ ranges::fill(eig.blocks, Matrix()); });
   }
   // Total number of states (symmetry taken into account)
   template <typename MF> auto count_states(MF && mult) const {
     return ranges::accumulate(*this, 0, {}, [mult](const auto &x) { const auto &[I, eig] = x; return mult(I)*eig.getnrstored(); });
   }
   [[nodiscard]] auto count_subspaces() const {    // Count non-empty subspaces
     return ranges::count_if(eigs(), [](const auto &eig) { return eig.getnrstored()>0; });
   }
   template<typename F, typename M> auto trace(F fnc, const double factor, M mult) const { // Tr[fnc(factor*E) exp(-factor*E)]
     return ranges::accumulate(*this, 0.0, {}, [fnc, factor, mult](const auto &x) { const auto &[I, eig] = x; return mult(I) * eig.trace(fnc, factor); });
   }
   template <typename MF>
     void states_report(MF && mult) const {
       fmt::print("Number of invariant subspaces: {}\n", count_subspaces());
       for (const auto &[I, eig]: *this) 
         if (eig.getnrstored()) 
           fmt::print("({}) {} states: {}\n", I.str(), eig.getnrstored(), eig.values.all_rel());
       fmt::print("Number of states (multiplicity taken into account): {}\n\n", count_states(mult));
     }
   void save(const size_t N, const Params &P) const {
     const std::string fn = P.workdir->unitaryfn(N);
     std::ofstream MATRIXF(fn, std::ios::binary | std::ios::out);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
     boost::archive::binary_oarchive oa(MATRIXF);
     oa << this->size();
     for(const auto &[I, eig]: *this) {
       oa << I;
       eig.save(oa);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn)); // Check after each write.
     }
   }
   void load(const size_t N, const Params &P, const bool remove_files = false) {
     const std::string fn = P.workdir->unitaryfn(N);
     std::ifstream MATRIXF(fn, std::ios::binary | std::ios::in);
     if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
     boost::archive::binary_iarchive ia(MATRIXF);
     const auto nr = read_one<size_t>(ia); // Number of subspaces
     for (const auto cnt : range0(nr)) {
       const auto inv = read_one<Invar>(ia);
       (*this)[inv].load(ia);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
     }
     if (remove_files) NRG::remove(fn);
   }
   void h5save(H5Easy::File &fd, const std::string &name, const bool write_absG) const {
     for (const auto &[I, eig]: *this)
       eig.h5save(fd, name + "/" + I.name(), write_absG);
   }
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) { 
     load(N, P, remove_files); 
   } // called from do_diag()
};

} // namespace

#endif
