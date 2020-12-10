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
   std::vector<t_eigen> corrected;
  public:
   void resize(const size_t size) { v.resize(size); }
   auto raw(const size_t i) const { return v[i]; }
   auto rel(const size_t i) const { return v[i]; }
   auto abs(const size_t i) const { my_assert(std::isfinite(scale)); return rel(i) * scale; }
   auto rel_zero(const size_t i) const { my_assert(std::isfinite(shift)); return rel(i)-shift; } // shift subtracted
   auto abs_zero(const size_t i) const { return (rel(i)-shift) * scale; }
   auto abs_T(const size_t i) const { my_assert(std::isfinite(T_shift)); return abs_zero(i) + T_shift; } // T_shift added
   auto abs_G(const size_t i) const { my_assert(std::isfinite(abs_GS_energy)); return abs_T(i) - abs_GS_energy; } // abs_GS_energy subtracted
   auto corr(const size_t i) const { return corrected[i]; }
   auto size() const { return v.size(); }
   auto lowest_rel() const { return v.front(); }
   auto highest_corr() const { return corrected.back(); }
   const auto & all_rel() const { return v; }
   auto all_rel_zero() const {
     my_assert(v.size() == 0 || std::isfinite(shift));
     return ranges::views::transform(v, [this](const auto x){ return x-shift; });
   }
   auto all_abs_zero() const {
     my_assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale; });
   }
   auto all_abs_T() const {
     my_assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale) && std::isfinite(T_shift)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale + T_shift; });
   }
   auto all_abs_G() const {
     my_assert(v.size() == 0 || (std::isfinite(shift) && std::isfinite(scale) && std::isfinite(T_shift) && std::isfinite(abs_GS_energy)));
     return ranges::views::transform(v, [this](const auto x){ return (x-shift) * scale + T_shift - abs_GS_energy; });
   }
   const auto & all_corr() const {
     return corrected;
   }
   void set(std::vector<t_eigen> in) { v = std::move(in); }
   void set_scale(const double scale_) { scale = scale_; }
   void set_shift(const double shift_) { shift = shift_; }
   void set_T_shift(const double T_shift_) { T_shift = T_shift_; }
   void set_abs_GS_energy(const double abs_GS_energy_) { abs_GS_energy = abs_GS_energy_; }  
   void set_corr(std::vector<t_eigen> in) { corrected = std::move(in); }
   auto has_abs() const { return std::isfinite(scale); }
   auto has_zero() const { return std::isfinite(shift); }
   auto has_abs_zero() const { return has_abs() && has_zero(); }
   auto has_abs_T() const { return has_abs_zero() && std::isfinite(T_shift); }
   auto has_abs_G() const { return has_abs_T() && std::isfinite(abs_GS_energy); }
   auto has_corr() const { return corrected.size() > 0; }
   void save(boost::archive::binary_oarchive &oa) const {
     oa << v << scale << shift << T_shift << abs_GS_energy << corrected;
   }
   void load(boost::archive::binary_iarchive &ia) {
     ia >> v >> scale >> shift >> T_shift >> abs_GS_energy >> corrected;
   }
   void h5save(H5Easy::File &fd, const std::string &name) const {
    h5_dump_vector(fd, name + "/value_orig", all_rel());
    if (has_zero())     h5_dump_vector(fd, name + "/value_zero",     all_rel_zero() | ranges::to_vector);
    if (has_abs_zero()) h5_dump_vector(fd, name + "/absenergy_zero", all_abs_zero() | ranges::to_vector);
    if (has_abs_T())    h5_dump_vector(fd, name + "/absenergy",      all_abs_T()    | ranges::to_vector);
    if (has_abs_G())    h5_dump_vector(fd, name + "/absenergyG",     all_abs_G()    | ranges::to_vector);
    if (has_corr())     h5_dump_vector(fd, name + "/value_corr",     all_corr());
  }
};

template <scalar S, typename t_matel = matel_traits<S>, typename Matrix = Matrix_traits<S>>
class Vectors {
  private:
    Matrix m;
  public:
    auto M() const { return m.size1(); }
    auto dim() const { return m.size2(); }
    void set(Matrix m_) { m = std::move(m_); my_assert(M() <= dim()); }
    const auto & get() const { return m; }
    void resize(const size_t size1, const size_t size2) { m.resize(size1, size2); }
    const auto & operator()() const { return m; }
    void standard_basis(const size_t size) {
      m = ublas::identity_matrix<t_matel>(size);
    }
    auto submatrix(const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2) const {
      return NRG::submatrix(m, r1, r2);
    }
    void shrink() {
      const auto d = dim();
      resize(0, d); // We keep the information about the dimensionality!!
      my_assert(M() == 0 && dim() == d);
    }
    void save(boost::archive::binary_oarchive &oa) const {
      NRG::save(oa, m);
    }
    void load(boost::archive::binary_iarchive &ia) {
      NRG::load(ia, m);
    }
    void h5save(H5Easy::File &fd, const std::string &name) const {
      h5_dump_matrix(fd, name + "/matrix", m);
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
    my_assert(i < blocks.size());
    blocks[i] = std::move(m);
  }
  const Matrix & get(const size_t i) const {
    my_assert(i < blocks.size());
    return blocks[i];
  }
  const Matrix & operator()(const size_t i) const { // 1-based MMA index, called from recalc_f()
    my_assert(1 <= i && i <= blocks.size());
    return blocks[i-1];
  }
  void truncate(const size_t nr) {
    for (auto &i : blocks) {
      my_assert(nr <= i.size1());
      i.resize(nr, i.size2());
    }
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
    for (const auto r1 : range0(M)) {
      for (const auto r2 : boost::irange(r1 + 1, M)) {
        S skpdt{};
        for (const auto j : range0(dim)) skpdt += conj_me(vec(r1, j)) * vec(r2, j);
        my_assert(num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON));
      }
    }
  }
};

// High-level representation of eigenvalues/eigenvectors
template <scalar S, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>> 
class Eigen {
public:
  Values<S> values;   // eigenvalues
  Vectors<S> vectors; // eigenvectors
  Blocks<S> U;        // eigenvectors in blocks
  Eigen() = default;
  explicit Eigen(const size_t M, const size_t dim) {
    my_assert(M <= dim);
    values.resize(M);
    vectors.resize(M, dim);
  }
  explicit Eigen(RawEigen<S> && raw) {
    values.set(std::move(raw.val));
    vectors.set(std::move(raw.vec));
  }
  [[nodiscard]] auto getnrcomputed() const { return values.size(); } // number of computed eigenpairs
  [[nodiscard]] auto getdim() const { return vectors.dim(); }        // valid also after the split_in_blocks_Eigen() call
 private:
  long nrpost = -1;   // number of eigenpairs after truncation (-1: keep all)
  long nrstored = -1; // number of eigenpairs currently held in store 
  [[nodiscard]] auto getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }     // number of states after truncation
 public:
  [[nodiscard]] auto getnrall() const { return getnrcomputed(); }                              // all = all computed
  [[nodiscard]] auto getnrkept() const { return getnrpost(); }                                 // # of kept states
  [[nodiscard]] auto getnrdiscarded() const { return getnrcomputed()-getnrpost(); }            // # of discarded states
  [[nodiscard]] auto getnrstored() const  { return nrstored == -1 ? values.size() : nrstored; }  // number of stored states
  [[nodiscard]] auto all() const { return range0(getnrcomputed()); }                           // iterator over all states
  [[nodiscard]] auto kept() const { return range0(getnrpost()); }                              // iterator over kept states
  [[nodiscard]] auto discarded() const { return boost::irange(getnrpost(), getnrcomputed()); } // iterator over discarded states
  [[nodiscard]] auto stored() const { return range0(getnrstored()); }                          // iterator over all stored states
  auto value_corr_kept() const { return ranges::subrange(values.all_corr().begin(), values.all_corr().begin() + getnrkept()); }
  auto value_corr_msr() const { return ranges::subrange(values.all_corr().begin(), values.all_corr().begin() + getnrstored()); } // range used in measurements (all or kept, depending on the moment of call)
  // Truncate to nrpost states.
  void truncate_prepare(const size_t nrpost_) {
    nrpost = nrpost_;
    my_assert(nrpost <= getnrcomputed());
  }
  void truncate_perform() {
    nrstored = nrpost;
    U.truncate(nrpost);
  }
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors form an identity matrix. This is used to
  // represent the spectral decomposition in the eigenbasis itself. Called when building DiagInfo from 'data' file.
  void diagonal(const std::vector<t_eigen> &v) {
    values.set(v);
    values.set_corr(v); // required, for matrix construction in the first step!
    values.set_shift(0.0); // required in the first step!
    vectors.standard_basis(v.size());
  }
  void subtract_Egs(const t_eigen Egs) {
    values.set_shift(Egs); 
  }
  void subtract_GS_energy(const t_eigen GS_energy) {
    values.set_abs_GS_energy(GS_energy);
  }
  auto diagonal_exp(const double factor) const { // produce a diagonal matrix with exp(-factor*E) diagonal elements, used in init_rho()
    const auto dim = getnrstored();
    auto m = Zero_matrix<S>(dim);
    for (const auto i: range0(dim)) 
      m(i, i) = exp(-values.corr(i) * factor); // corrected eigenvalues!
    return m;
  }
  template<typename F> auto trace(F fnc, const double factor) const { // Tr[fnc(factor*E) exp(-factor*E)]
    return ranges::accumulate(values.all_rel_zero(), 0.0, {}, [fnc, factor](const auto x) { return fnc(factor*x) * exp(-factor*x); });
  }
  void clear_eigenvectors() {
    vectors.shrink();
    U.clear();
  }
  void save(boost::archive::binary_oarchive &oa) const {
    values.save(oa);
    vectors.save(oa);
    oa << nrpost << nrstored;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    values.load(ia);
    vectors.load(ia);
    ia >> nrpost >> nrstored;
  }
  void h5save(H5Easy::File &fd, const std::string &name) const {
    values.h5save(fd, name);
    vectors.h5save(fd, name);
    h5_dump_scalar(fd, name + "/nrkept", getnrkept());
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
     ranges::for_each(eigs(), [Egs](auto &eig) { eig.subtract_Egs(Egs); });
   }
   t_eigen Egs_subtraction() {
     const auto Egs = find_groundstate();
     subtract_Egs(Egs);
     return Egs;
   }
   void subtract_GS_energy(const t_eigen GS_energy) {
     ranges::for_each(eigs(), [GS_energy](auto &eig) { eig.subtract_GS_energy(GS_energy); });
   }
   std::vector<t_eigen> sorted_energies_rel_zero() const {
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs()) 
       energies.insert(energies.end(), eig.values.all_rel_zero().begin(), eig.values.all_rel_zero().end());
     return energies | ranges::move | ranges::actions::sort;
   }
   std::vector<t_eigen> sorted_energies_corr() const {
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs()) 
       energies.insert(energies.end(), eig.values.all_corr().begin(), eig.values.all_corr().end());
     return energies | ranges::move | ranges::actions::sort;
   }
   void dump_energies(std::ostream &F) const {
     for (const auto &[I, eig]: *this)
       F << "Subspace: " << I << std::endl << eig.values.all_rel() << std::endl;
   }
   void truncate_perform() {
     ranges::for_each(eigs(), &Eigen<S>::truncate_perform); // Truncate subspace to appropriate size
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
     ranges::for_each(eigs(), &Eigen<S>::clear_eigenvectors);
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
   void h5save(H5Easy::File &fd, const std::string &name) const {
     for (const auto &[I, eig]: *this) eig.h5save(fd, name + "/" + I.name());
   }
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) { 
     load(N, P, remove_files); 
   } // called from do_diag()
};

} // namespace

#endif
