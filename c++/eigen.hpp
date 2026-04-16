#ifndef _eigen_hpp_
#define _eigen_hpp_

#include <vector>
#include <string>
#include <limits> // quiet_NaN
#include <stdexcept>

#include <boost/range/adaptor/map.hpp>

#include <range/v3/all.hpp>

#include "basicio.hpp"
#include "portabil.hpp"
#include "traits.hpp"
#include "invar.hpp"
#include "params.hpp"
#include "step.hpp"
#include "h5.hpp"
#include "numerics.hpp"

#include <fmt/format.h>

namespace NRG {

template<typename T>
std::vector<T> negate_copy(const std::vector<T>& v) {
  std::vector<T> out(v.size());
  std::transform(v.begin(), v.end(), out.begin(),
                 [](const T& x) { return -x; });
  return out;
}

template<typename T>
std::vector<T> abs_copy(const std::vector<T>& v) {
  std::vector<T> out(v.size());
  std::transform(v.begin(), v.end(), out.begin(),
                 [](const T& x) { return std::abs(x); });
  return out;
}

template<typename T>
void negate_inplace(std::vector<T>& v) {
  for (auto& x : v)
    x = -x;
}

template<typename T>
void abs_inplace(std::vector<T>& v) {
  for (auto& x : v)
    x = std::abs(x);
}

template<typename T>
void shift_inplace(std::vector<T>& v, const T shift) {
  for (auto& x : v)
    x -= shift; // note that we do a subtraction!
}

// Storage container for eigenvalues. Vector v contains the raw eigenvalues as computed in the Hamiltonian
// diagonalisation. If scale, shift, T_shift and/or GS_energy parameters are defined, then one has also access to various
// derived quantities. 'Relative' means in the units of the current NRG shell. 'Absolute' means in the units of
// half-bandwidth (or some other general scale). abs_T is the total energy of the state. abs_G is the total energy
// of the state referenced to the absolute many-body ground state of the whole system. 'Corrected' means the eigenvalue
// corrected for the floating-point round-off errors by fixing the small splittings.
template<scalar S, typename t_eigen = eigen_traits<S>>
class Values {
  private:
   std::vector<t_eigen> v;
   double scale = std::numeric_limits<double>::quiet_NaN();
   double shift = std::numeric_limits<double>::quiet_NaN();
   double T_shift = std::numeric_limits<double>::quiet_NaN();
   double abs_GS_energy = std::numeric_limits<double>::quiet_NaN();
   std::vector<t_eigen> corrected; // these are referenced to Egs in current step (in the same way as rel_zero)
   std::vector<t_eigen> c; // criterion, value of "cost" function for truncation; constrols the sort order of states!!
   double c_shift = std::numeric_limits<double>::quiet_NaN(); // shift for criterion 'c'
  public:
   std::vector<t_eigen> & ref_v() { return v; }
   std::vector<t_eigen> & ref_corrected() { return corrected; }
   std::vector<t_eigen> & ref_c() { return c; }
   void resize(const size_t size) { 
     v.resize(size); 
     c.resize(size); 
   }
   [[nodiscard]] auto raw(const size_t i) const { return v[i]; }
   [[nodiscard]] auto rel(const size_t i) const { return v[i]; }
   [[nodiscard]] auto abs(const size_t i) const { assert(std::isfinite(scale)); return rel(i) * scale; }
   [[nodiscard]] auto rel_zero(const size_t i) const { assert(std::isfinite(shift)); return rel(i)-shift; } // shift subtracted
   [[nodiscard]] auto abs_zero(const size_t i) const { return (rel(i)-shift) * scale; }
   [[nodiscard]] auto abs_T(const size_t i) const { assert(std::isfinite(T_shift)); return abs_zero(i) + T_shift; } // T_shift added
   [[nodiscard]] auto abs_G(const size_t i) const { assert(std::isfinite(abs_GS_energy)); return abs_T(i) - abs_GS_energy; } // abs_GS_energy subtracted
   [[nodiscard]] auto corr(const size_t i) const { return corrected[i]; }
   [[nodiscard]] auto size() const noexcept { return v.size(); }
   [[nodiscard]] auto lowest_rel() const { return v.front(); }
   [[nodiscard]] auto lowest_crit() const { return c.front(); }
   [[nodiscard]] auto highest_crit() const { return c.back(); }
   void report(bool verbose = false) const {
     assert(v.size() > 0);
     std::cout << "v[0]=" << v.front() << " corrected[0]=" << corrected.front() << " crit[0]=" << c.front() << std::endl;
     if (verbose)
       std::cout << "scale=" << scale << " shift=" << shift << " T_shift=" << T_shift << " c_shift=" << c_shift << std::endl;
   }
   [[nodiscard]] const auto & all_rel() const noexcept { return v; }
   [[nodiscard]] const auto & all_crit() const noexcept { return c; }
   [[nodiscard]] auto all_rel_zero() const noexcept {
     assert(v.size() == 0 || std::isfinite(shift));
     return ranges::views::transform(v, [this](const auto x){ return x-shift; });
   }
   [[nodiscard]] auto all_abs_zero() const noexcept {
     assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale; });
   }
   [[nodiscard]] auto all_abs_T() const noexcept {
     assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale) && std::isfinite(T_shift)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale + T_shift; });
   }
   [[nodiscard]] auto all_abs_G() const noexcept {
     assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale) && std::isfinite(T_shift) && std::isfinite(abs_GS_energy)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale + T_shift - abs_GS_energy; });
   }
   [[nodiscard]] const auto & all_corr() const noexcept {
     return corrected;
   }
   void set(std::vector<t_eigen> in) { v = std::move(in); }
   void set_scale(const double scale_) { scale = scale_; }
   void set_shift(const double shift_) { shift = shift_; }
   void set_c_shift(const double c_shift_) { c_shift = c_shift_; }
   void set_T_shift(const double T_shift_) { T_shift = T_shift_; }
   void set_abs_GS_energy(const double abs_GS_energy_) { abs_GS_energy = abs_GS_energy_; }
   void set_corr(std::vector<t_eigen> in) { corrected = std::move(in); }
   void do_shift(const double shift) { shift_inplace(v,shift); }
   void do_c_shift(const double shift) { shift_inplace(c, shift); }
   void do_corr_shift(const double shift) { shift_inplace(corrected, shift); }
   [[nodiscard]] auto crit(const size_t i) const { return c[i]; }
   void set_crit(const size_t i, const t_eigen x) { c[i] = x; }
   void crit_copy_raw() { c = v; }
   void crit_copy_corr() { c = corrected; }
   void crit_negate() { negate_inplace(c); }
   void crit_abs() { abs_inplace(c); }
   [[nodiscard]] auto has_abs() const noexcept { return std::isfinite(scale); }
   [[nodiscard]] auto has_zero() const noexcept { return std::isfinite(shift); }
   [[nodiscard]] auto has_abs_zero() const noexcept { return has_abs() && has_zero(); }
   [[nodiscard]] auto has_abs_T() const noexcept { return has_abs_zero() && std::isfinite(T_shift); }
   [[nodiscard]] auto has_abs_G() const noexcept { return has_abs_T() && std::isfinite(abs_GS_energy); }
   [[nodiscard]] auto has_corr() const noexcept { return corrected.size() > 0; }
   [[nodiscard]] auto has_crit() const noexcept { return c.size() > 0; }
   void save(boost::archive::binary_oarchive &oa) const {
     oa << v << scale << shift << T_shift << abs_GS_energy << corrected << c;
   }
   void load(boost::archive::binary_iarchive &ia) {
     ia >> v >> scale >> shift >> T_shift >> abs_GS_energy >> corrected >> c;
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
     h5_dump_vector(fd, name + "/value_orig", all_rel());
     if (has_zero())     h5_dump_vector(fd, name + "/value_zero",     all_rel_zero() | ranges::to_vector);
     if (has_abs_zero()) h5_dump_vector(fd, name + "/absenergy_zero", all_abs_zero() | ranges::to_vector);
     if (has_abs_T())    h5_dump_vector(fd, name + "/absenergy",      all_abs_T()    | ranges::to_vector);
     if (has_abs_G())    h5_dump_vector(fd, name + "/absenergyG",     all_abs_G()    | ranges::to_vector);
     if (has_corr())     h5_dump_vector(fd, name + "/value_corr",     all_corr());
     if (has_crit())     h5_dump_vector(fd, name + "/crit",           all_crit());
   }
};

template <scalar S, typename t_matel = matel_traits<S>, typename Matrix = Matrix_traits<S>>
class Vectors {
  private:
    Matrix m;
    size_t _dim = 0;
  public:
    Matrix & ref_m() { return m; }
    [[nodiscard]] auto M() const noexcept { return size1(m); }
    [[nodiscard]] auto dim() const noexcept { return _dim; }
    void set(Matrix m_) {
      m = std::move(m_);
      _dim = size2(m);
      assert(is_unitary<S>(m));
      assert(M() <= dim());
    }
    [[nodiscard]] const auto & get() const noexcept { return m; }
    [[nodiscard]] const auto & operator()() const noexcept { return m; }
    void standard_basis(const size_t size) { m = id_matrix<t_matel>(size); }
    auto submatrix_const(const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2) const {
      return NRG::submatrix_const(m, r1, r2);
    }
    void resize(const size_t new_size1, const size_t new_size2) { // called from Eigen(int,int) constructor
      m.resize(new_size1, new_size2); // non-conserving resize
      _dim = new_size2;
    }
    void shrink() {
      m = Matrix(0, 0);
    }
    void save(boost::archive::binary_oarchive &oa) const {
      NRG::save(oa, m);
      oa << _dim;
    }
    void load(boost::archive::binary_iarchive &ia) {
      m = NRG::load<S>(ia);
      ia >> _dim;
    }
    void h5save(H5Easy::File &fd, const std::string &name) const {
      h5_dump_matrix(fd, name + "/matrix", m);
    }
   void dump(std::ostream &F) const {
     for (size_t i = 0; i < M(); i++) {
       F << "vec(" << i << ")=[";
       double sum = 0.0;
       for (size_t j = 0; j < dim(); j++) {
         F << m(i,j) << (j != dim()-1 ? ", " : "");
         sum += pow(abs(m(i,j)),2);
       }
       const double diff = sum-1.0;
       F << "] norm-1=" << diff << std::endl;
     }
   }
};

// Eigenvectors separated according to the invariant subspace from which they originate.
// Required for using efficient BLAS routines when performing recalculations of the matrix elements.
template <scalar S, typename Matrix = Matrix_traits<S>>
class Blocks {
private:
  std::vector<Matrix> blocks;
public:
  void resize(const size_t nr) {
    blocks.resize(nr);
  }
  void set(const size_t i, Matrix m) {
    assert(i < blocks.size());
    blocks[i] = std::move(m);
  }
  [[nodiscard]] bool is_unitary() const noexcept {
    return NRG::is_unitary_blocks<S>(blocks);
  }
  const Matrix & get(const size_t i) const {
    assert(i < blocks.size());
    return blocks[i];
  }
  const Matrix & operator()(const size_t i) const { // 1-based MMA index, called from recalc_f()
    assert(1 <= i && i <= blocks.size());
    return blocks[i-1];
  }
  void truncate(const size_t nr) {
    for (auto &i : blocks) {
      assert(nr <= nrvec(i));
      NRG::resize(i, nr, dim(i)); // conserving matrix resize
    }
  }
  void save(boost::archive::binary_oarchive &oa) const {
    oa << blocks.size();
    for (const auto &b: blocks) NRG::save(oa, b);
  }
  void load(boost::archive::binary_iarchive &ia) {
    resize(read_one<size_t>(ia));
    for (auto &b: blocks) b = NRG::load<S>(ia);
  }
  void clear() { ranges::fill(blocks, Matrix()); }
};

// Result of a diagonalisation: eigenvalues and eigenvectors
template <scalar S, typename RVector = RVector_traits<S>, typename Matrix = Matrix_traits<S>>
struct RawEigen {
public:
  RVector val;
  Matrix vec;
  RawEigen() = default;
  RawEigen(const size_t M, const size_t dim) {
    assert(M <= dim);
    val.resize(M);
    vec.resize(M, dim); // non-conserving matrix resize
  }
  [[nodiscard]] auto getnrcomputed() const noexcept { assert(val.size() == nrvec(vec)); return val.size(); } // nr eigenvalue/eigenvector pairs
  [[nodiscard]] auto getdim() const noexcept { return dim(vec); } // matrix dimension (length of eigenvectors)
  void dump_eigenvalues(const size_t max_nr = std::numeric_limits<size_t>::max(), std::ostream &F = std::cout) const {
    F << "eig= " << std::setprecision(std::numeric_limits<double>::max_digits10);
    ranges::for_each_n(val.begin(), std::min(val.size(), max_nr), [&F](const double x) { F << x << ' '; });
    F << std::endl;
  }
  void check_diag() const { NRG::check_diag<S>(val, vec); }
};

template <typename T, typename S>
void perform_sort_by_c(std::vector<T>& c,
                       std::vector<T>& v,
                       std::vector<T>& corrected,
                       Eigen::Matrix<S, -1, -1, Eigen::RowMajor>& m) {
  const std::size_t n = c.size();
  assert(v.size() == n);
  assert(corrected.size() == n);
  assert(static_cast<std::size_t>(m.rows()) == n);
  // Build permutation: p[i] = original index of the i-th smallest element of c
  std::vector<std::size_t> p(n);
  std::iota(p.begin(), p.end(), 0);
  std::sort(p.begin(), p.end(),
            [&](std::size_t a, std::size_t b) {
              return c[a] < c[b];
            });
  // Apply permutation to c and v
  std::vector<T> c_sorted;
  std::vector<T> v_sorted;
  std::vector<T> corrected_sorted;
  c_sorted.reserve(n);
  v_sorted.reserve(n);
  corrected_sorted.reserve(n);
  for (std::size_t i = 0; i < n; ++i) {
    c_sorted.push_back(std::move(c[p[i]]));
    v_sorted.push_back(std::move(v[p[i]]));
    corrected_sorted.push_back(std::move(corrected[p[i]]));
  }
  // Apply permutation to rows of m
  Eigen::Matrix<S, -1, -1, Eigen::RowMajor> m_sorted(m.rows(), m.cols());
  for (std::size_t i = 0; i < n; ++i)
    m_sorted.row(static_cast<Eigen::Index>(i)) = m.row(static_cast<Eigen::Index>(p[i]));
  c = std::move(c_sorted);
  v = std::move(v_sorted);
  corrected = std::move(corrected_sorted);
  m = std::move(m_sorted);
}

// High-level representation of eigenvalues/eigenvectors
template <scalar S, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>>
class Eigen {
public:
  Values<S> values;   // eigenvalues
  Vectors<S> vectors; // eigenvectors
  Blocks<S> U;        // eigenvectors in blocks
  Eigen() = default;
  explicit Eigen(const size_t M, const size_t dim) {
    assert(M <= dim);
    values.resize(M);
    vectors.resize(M, dim); // non-conserving matrix-resize
  }
  explicit Eigen(RawEigen<S> && raw, const Step &step) {
    values.set(std::move(raw.val));
    values.crit_copy_raw(); // conventionally, cost = energy for truncation criterion
    vectors.set(std::move(raw.vec));
    last = step.last();
  }
  // Called when building DiagInfo from 'data' file.
  explicit Eigen(const std::vector<t_eigen> &v, const double scale) {
    diagonal(v);
    values.set_corr(v); // required for matrix construction in the first step!
    values.set_shift(0.0); // required in the first step!
    values.set_scale(scale);
    values.crit_copy_raw();
    last = false; // true only for ZBW, but support for that case has been removed in 2026
  }
  [[nodiscard]] auto getnrcomputed() const noexcept { return values.size(); } // number of computed eigenpairs
  [[nodiscard]] auto getdim() const noexcept { return vectors.dim(); }        // valid also after the split_in_blocks_Eigen() call
 private:
  long nrpost = -1;   // number of eigenpairs after truncation (-1: keep all)
  long nrstored = -1; // number of eigenpairs currently held in store
  bool last = false; // eigenspectrum from the last step of the NRG iteration
  [[nodiscard]] auto getnrpost() const noexcept { return nrpost == -1 ? getnrcomputed() : nrpost; }     // number of states after truncation
  [[nodiscard]] auto boundary() const noexcept { return last ? 0ul : getnrkept(); } // for FDM
 public:
  [[nodiscard]] auto getnrall() const noexcept { return getnrcomputed(); }                              // all = all computed
  [[nodiscard]] auto getnrkept() const noexcept { return getnrpost(); }                                 // # of kept states
  [[nodiscard]] auto getnrdiscarded() const noexcept { return getnrcomputed()-getnrpost(); }            // # of discarded states
  [[nodiscard]] auto getnrstored() const noexcept { return nrstored == -1 ? values.size() : nrstored; }  // number of stored states
  [[nodiscard]] auto all() const noexcept { return range0(getnrcomputed()); }                           // iterator over all states
  [[nodiscard]] auto kept() const noexcept { return range0(getnrpost()); }                              // iterator over kept states
  [[nodiscard]] auto discarded() const noexcept { return boost::irange(getnrpost(), getnrcomputed()); } // iterator over discarded states
  [[nodiscard]] auto stored() const noexcept { return range0(getnrstored()); }                          // iterator over all stored states
  // Ranges for FDM algorithm with different semantics of D/K states for the last step
  [[nodiscard]] auto Drange() const noexcept { return boost::irange(boundary(), getnrall()); }
  [[nodiscard]] auto Krange() const noexcept { return boost::irange(0ul, boundary()); }
  [[nodiscard]] auto value_corr_kept() const noexcept { return ranges::subrange(values.all_corr().begin(), values.all_corr().begin() + getnrkept()); }
  [[nodiscard]] auto value_corr_msr() const noexcept { return ranges::subrange(values.all_corr().begin(), values.all_corr().begin() + getnrstored()); } // range used in measurements (all or kept, depending on the moment of call)
  // Reordering
  void do_sort_by_c() {
    perform_sort_by_c(values.ref_c(), values.ref_v(), values.ref_corrected(), vectors.ref_m());
  }
  // Truncate to nrpost states.
  void truncate_prepare(const size_t nrpost_) {
    nrpost = nrpost_;
    assert(nrpost <= getnrcomputed());
  }
  void truncate_perform() {
    nrstored = nrpost;
    U.truncate(nrpost);
  }
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors form an identity matrix. This is used to
  // represent the spectral decomposition in the eigenbasis itself.
  void diagonal(const std::vector<t_eigen> &v) {
    values.set(v);
    vectors.standard_basis(v.size());
  }
  void set_shift_Egs(const t_eigen Egs) {
    values.set_shift(Egs);
  }
  void shift(const t_eigen Egs, const t_eigen Clw) {
    values.do_shift(Egs);
    values.do_corr_shift(Egs);
    values.do_c_shift(Clw);
  }
  void subtract_GS_energy(const t_eigen GS_energy) {
    values.set_abs_GS_energy(GS_energy);
  }
  [[nodiscard]] auto diagonal_exp(const double factor) const noexcept { // produce a diagonal matrix with exp(-factor*E) diagonal elements, used in init_rho()
    const auto dim = getnrstored();
    auto m = zero_matrix<S>(dim);
    for (const auto i: range0(dim))
      m(i, i) = exp(-values.corr(i) * factor); // corrected eigenvalues!
    return m;
  }
  template<typename F>
  [[nodiscard]] auto trace(F fnc, const double factor) const noexcept { // Tr[fnc(factor*E) exp(-factor*E)]
    return ranges::accumulate(values.all_rel_zero(), 0.0, {}, [fnc, factor](const auto x) { return fnc(factor*x) * exp(-factor*x); });
  }
  void clear_eigenvectors() {
    vectors.shrink();
    U.clear();
  }
  void save(boost::archive::binary_oarchive &oa) const {
    values.save(oa);
    vectors.save(oa);
    U.save(oa);
    oa << nrpost << nrstored << last;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    values.load(ia);
    vectors.load(ia);
    U.load(ia);
    ia >> nrpost >> nrstored >> last;
  }
  void h5save(H5Easy::File &fd, const std::string &name, const bool save_vectors = true) const {
    values.h5save(fd, name);
    if (save_vectors) vectors.h5save(fd, name);
    h5_dump_scalar(fd, name + "/nrkept", getnrkept());
  }
};

// Full information after diagonalizations (eigenspectra in all subspaces)
template <scalar S, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>>
class DiagInfo : public std::map<Invar, Eigen<S>> {
 public:
   explicit DiagInfo() = default;
   DiagInfo(std::istream &fdata, const size_t nsubs, const Params &P) {
     skip_comments(fdata);
     for ([[maybe_unused]] const auto i : range1(nsubs)) {
       const auto I = read_one<Invar>(fdata);
       auto energies = read_std_vector<t_eigen>(fdata);
       if (!(P.data_has_rescaled_energies || P.absolute))
         for (auto &x : energies) x /= P.SCALE(P.Ninit); // rescale to the suitable energy scale
       (*this)[I] = Eigen<S>(energies, P.absolute ? 1.0 : P.SCALE(P.Ninit));
     }
     my_assert(this->size() == nsubs);
   }
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) {  // called from do_diag()
       load(N, P, remove_files);
   }
   [[nodiscard]] auto subspaces() const noexcept { return *this | boost::adaptors::map_keys; }
   [[nodiscard]] auto eigs() const noexcept { return *this | boost::adaptors::map_values; }
   [[nodiscard]] auto eigs() noexcept { return *this | boost::adaptors::map_values; }
   [[nodiscard]] auto find_Egs() const {
     const auto [Iground, eig] = *ranges::min_element(*this, {}, [](const auto &a) { return a.second.values.lowest_rel(); });
     const auto Egs = eig.values.lowest_rel();
     return Egs;
   }
   [[nodiscard]] auto find_Clw() const {
     const auto [Ilowest, eig] = *ranges::min_element(*this, {}, [](const auto &a) { return a.second.values.lowest_crit(); });
     const auto Clw = eig.values.lowest_crit();
     return Clw;
   }
   void set_shift_Egs(const t_eigen Egs) {
     ranges::for_each(eigs(), [Egs](auto &eig) { eig.set_shift_Egs(Egs); });
   }
   // shifts corr(ected) and criterion; raw values (v) remain untouched
   void shift(const t_eigen Egs, const t_eigen Clw) {
     ranges::for_each(eigs(), [Egs, Clw](auto &eig) { eig.shift(Egs, Clw); });
   }
   void subtract_GS_energy(const t_eigen GS_energy) {
     ranges::for_each(eigs(), [GS_energy](auto &eig) { eig.subtract_GS_energy(GS_energy); });
   }
   [[nodiscard]] std::vector<t_eigen> sorted_energies_rel_zero() const {
     std::vector<t_eigen> all;
     for (const auto &eig: eigs())
       all.insert(all.end(), eig.values.all_rel_zero().begin(), eig.values.all_rel_zero().end());
     return all | ranges::move | ranges::actions::sort;
   }
   [[nodiscard]] std::vector<t_eigen> sorted_criterion_values() const {
     std::vector<t_eigen> all;
     for (const auto &eig: eigs())
       all.insert(all.end(), eig.values.all_crit().begin(), eig.values.all_crit().end());
     return all | ranges::move | ranges::actions::sort;
   }
   void sort_by_c() {
     for (auto &eig: eigs())
       eig.do_sort_by_c();
   }
   void negate_c() {
     for (auto &eig: eigs())
       eig.values.crit_negate();
   }
   void abs_c() {
     for (auto &eig: eigs())
       eig.values.crit_abs();
   }
   void report(const bool verbose = false) const {
     std::cout << "DiagInfo report" << std::endl;
     for (auto &[I, eig]: *this) {
       std::cout << "I=" << I << std::endl;
       eig.values.report(true);
     }
   }
   void truncate_perform() {
     ranges::for_each(eigs(), &Eigen<S>::truncate_perform); // Truncate subspace to appropriate size
   }
   [[nodiscard]] auto size_subspace(const Invar &I) const {
     const auto f = this->find(I);
     const auto val = f != this->cend() ? f->second.getnrstored() : 0;
     my_assert(0 <= val && val <= 1000000); // bug trap
     return val;
   }
   [[nodiscard]] auto subs(const Invar &I1, const Invar &I2) const {
     my_assert(this->find(I1) != this->cend() && this->find(I2) != this->cend()); // bug trap: check for existance
     return std::make_pair(this->at(I1), this->at(I2));
   }
   [[nodiscard]] auto dims(const Invar &I1, const Invar &I2) const { // Determines matrix sizes for operators (# stored)
     return std::make_pair(size_subspace(I1), size_subspace(I2));
   }
   void clear_eigenvectors() {
     ranges::for_each(eigs(), &Eigen<S>::clear_eigenvectors);
   }
   // Total number of states (symmetry taken into account)
   template <typename MF>
   [[nodiscard]] auto count_states(MF && mult) const noexcept {
     return ranges::accumulate(*this, 0, {}, [mult](const auto &x) { const auto &[I, eig] = x; return mult(I)*eig.getnrstored(); });
   }
   [[nodiscard]] auto count_subspaces() const noexcept {    // Count non-empty subspaces
     return ranges::count_if(eigs(), [](const auto &eig) { return eig.getnrstored()>0; });
   }
   template<typename F, typename M>
   [[nodiscard]] auto trace(F fnc, const double factor, M mult) const noexcept { // Tr[fnc(factor*E) exp(-factor*E)]
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
     for ([[maybe_unused]] const auto cnt : range0(nr)) {
       const auto inv = read_one<Invar>(ia);
       (*this)[inv].load(ia);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
     }
     if (remove_files) NRG::remove(fn);
   }
   void h5save(H5Easy::File &fd, const std::string &name, const bool save_vectors = true) const {
     for (const auto &[I, eig]: *this) eig.h5save(fd, name + "/" + I.name(), save_vectors);
   }
};

inline auto rescaled(std::vector<double> v, double a) { std::ranges::for_each(v, [a](double& x){ x *= a; }); return v; }

// This dumps all energies. These are the 'relative' (scaled), but *not* shifted to zero!
template<scalar S>
void dump_all_energies(const DiagInfo<S> &diag, const double rescaled_by, std::ostream &F, const Params &P)  {
  for (const auto &[I, eig]: diag) {
    F << "Subspace: " << I << std::endl;
    F << rescaled(eig.values.all_rel(), rescaled_by) << std::endl;
    if (P.dumpcorr)
      F << "corr=" << rescaled(eig.values.all_corr(), rescaled_by) << std::endl;
    if (P.dumpcrit)
      F << "crit=" << rescaled(eig.values.all_crit(), rescaled_by) << std::endl;
  }
}

template<scalar S>
void dump_all_states(const DiagInfo<S> &diag, const double rescaled_by, std::ostream &F, const Params &P)  {
  for (const auto &[I, eig]: diag) {
    F << "Subspace: " << I << std::endl;
    F << "Energies (rel): " << rescaled(eig.values.all_rel(), rescaled_by) << std::endl;
    F << "Vectors:" << std::endl;
    eig.vectors.dump(F);
    F << std::endl; // empty line
  }
}

} // namespace

#endif
