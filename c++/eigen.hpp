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
   [[nodiscard]] auto highest_corr() const { return corrected.back(); }
   [[nodiscard]] const auto & all_rel() const noexcept { return v; }
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
   void set_T_shift(const double T_shift_) { T_shift = T_shift_; }
   void set_abs_GS_energy(const double abs_GS_energy_) { abs_GS_energy = abs_GS_energy_; }
   void set_corr(std::vector<t_eigen> in) { corrected = std::move(in); }
   [[nodiscard]] auto has_abs() const noexcept { return std::isfinite(scale); }
   [[nodiscard]] auto has_zero() const noexcept { return std::isfinite(shift); }
   [[nodiscard]] auto has_abs_zero() const noexcept { return has_abs() && has_zero(); }
   [[nodiscard]] auto has_abs_T() const noexcept { return has_abs_zero() && std::isfinite(T_shift); }
   [[nodiscard]] auto has_abs_G() const noexcept { return has_abs_T() && std::isfinite(abs_GS_energy); }
   [[nodiscard]] auto has_corr() const noexcept { return corrected.size() > 0; }
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
    [[nodiscard]] auto M() const noexcept { return size1(m); }
    [[nodiscard]] auto dim() const noexcept { return size2(m); }
    void set(Matrix m_) {
      m = std::move(m_);
      assert(is_unitary<S>(m));
      assert(M() <= dim());
    }
    [[nodiscard]] const auto & get() const noexcept { return m; }
    [[nodiscard]] const auto & operator()() const noexcept { return m; }
    void standard_basis(const size_t size) { m = id_matrix<t_matel>(size); }
    auto submatrix_const(const std::pair<size_t,size_t> &r1, const std::pair<size_t,size_t> &r2) const {
      return NRG::submatrix_const(m, r1, r2);
    }
    void resize(const size_t new_size1, const size_t new_size2) {
      m.resize(new_size1, new_size2); // non-conserving resize
    }
    void shrink() {
      const auto d = dim();
      resize(0, d); // Shrink to zero size, but keep the information about the dimensionality!!
      assert(M() == 0 && dim() == d);
    }
    void save(boost::archive::binary_oarchive &oa) const {
      NRG::save(oa, m);
    }
    void load(boost::archive::binary_iarchive &ia) {
      m = NRG::load<S>(ia);
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
    vectors.set(std::move(raw.vec));
    last = step.last();
  }
  // Called when building DiagInfo from 'data' file.
  explicit Eigen(const std::vector<t_eigen> &v, const double scale, const bool last_step) {
    diagonal(v);
    values.set_corr(v); // required for matrix construction in the first step!
    values.set_shift(0.0); // required in the first step!
    values.set_scale(scale);
    last = last_step;
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
  void subtract_Egs(const t_eigen Egs) {
    values.set_shift(Egs);
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
    oa << nrpost << nrstored << last;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    values.load(ia);
    vectors.load(ia);
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
     for (const auto i : range1(nsubs)) {
       const auto I = read_one<Invar>(fdata);
       auto energies = read_std_vector<t_eigen>(fdata);
       if (!(P.data_has_rescaled_energies || P.absolute))
         for (auto &x : energies) x /= P.SCALE(P.Ninit); // rescale to the suitable energy scale
       (*this)[I] = Eigen<S>(energies, P.absolute ? 1.0 : P.SCALE(P.Ninit), P.ZBW());
     }
     my_assert(this->size() == nsubs);
   }
   [[nodiscard]] auto subspaces() const noexcept { return *this | boost::adaptors::map_keys; }
   [[nodiscard]] auto eigs() const noexcept { return *this | boost::adaptors::map_values; }
   [[nodiscard]] auto eigs() noexcept { return *this | boost::adaptors::map_values; }
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
   [[nodiscard]] std::vector<t_eigen> sorted_energies_rel_zero() const {
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs())
       energies.insert(energies.end(), eig.values.all_rel_zero().begin(), eig.values.all_rel_zero().end());
     return energies | ranges::move | ranges::actions::sort;
   }
   [[nodiscard]] std::vector<t_eigen> sorted_energies_corr() const {
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
     for (const auto cnt : range0(nr)) {
       const auto inv = read_one<Invar>(ia);
       (*this)[inv].load(ia);
       if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
     }
     if (remove_files) NRG::remove(fn);
   }
   void h5save(H5Easy::File &fd, const std::string &name, const bool save_vectors = true) const {
     for (const auto &[I, eig]: *this) eig.h5save(fd, name + "/" + I.name(), save_vectors);
   }
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) {
     load(N, P, remove_files);
   } // called from do_diag()
};

} // namespace

#endif
