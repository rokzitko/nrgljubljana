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
class Values : public std::vector<t_eigen> {
  private:
   double scale = std::numeric_limits<double>::quiet_NaN();
   double shift = std::numeric_limits<double>::quiet_NaN();
   double GS_energy = std::numeric_limits<double>::quiet_NaN();
  public:
 //  Values() {}
 //  explicit Values(const double scale) : scale(scale) {}
   auto rel(const size_t i) const { return (*this)[i]; }
   auto abs(const size_t i) const { return rel(i) * scale; }
   auto rel_zero(const size_t i) const { return rel(i)-shift; }
   auto abs_zero(const size_t i) const { return (rel(i)-shift) * scale; }
   auto absG(const size_t i) const { return rel(i)*scale - GS_energy; }
   auto lowest_rel() const { return this->front(); }
   auto getnrcomputed() const { return this->size(); }
   void set_scale(const double scale_) { scale = scale_; }
   void set_shift(const double shift_) { shift = shift_; }
   void set_GS_energy(const double GS_energy_) { GS_energy = GS_energy_; }
   void save(boost::archive::binary_oarchive &oa) const {
     std::vector<t_eigen> tmp(this->begin(), this->end());
     oa << tmp << scale << shift << GS_energy;
   }
   void load(boost::archive::binary_iarchive &ia) {
     std::vector<t_eigen> tmp;
     ia >> tmp >> scale >> shift >> GS_energy;
     this->resize(tmp.size());
     std::copy(tmp.begin(), tmp.end(), this->begin());
   }
   void copy(const std::vector<t_eigen> &v, const size_t M) {
     this->resize(M);
     std::copy(v.begin(), v.begin() + M, this->begin());
   }
   void copy(const std::vector<t_eigen> &v) { this->copy(v, v.size()); }
};

// Result of a diagonalisation: eigenvalues and eigenvectors
template <scalar S, typename EVEC = evec_traits<S>, typename Matrix = Matrix_traits<S>, typename t_eigen = eigen_traits<S>> 
class Eigen {
public:
  Values<S> values; // eigenvalues
  EVEC value_orig;  // eigenvalues as computed YYY
  EVEC value_zero;  // eigenvalues with Egs subtracted
  Matrix matrix;    // eigenvectors
  Eigen() = default;
  Eigen(const size_t nr, const size_t dim) {
    my_assert(nr <= dim);
    values.resize(nr); // YYY
    value_orig.resize(nr);
    value_zero.resize(nr);
    matrix.resize(nr, dim);
  }
  [[nodiscard]] auto getnrcomputed() const { return value_orig.size(); } // number of computed eigenpairs
  [[nodiscard]] auto getdim() const { return matrix.size2(); }           // valid also after the split_in_blocks_Eigen() call
 private:
  long nrpost = -1;  // number of eigenpairs after truncation (-1: keep all)
  [[nodiscard]] auto getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }     // number of states after truncation
 public:
  [[nodiscard]] auto getnrall() const { return getnrcomputed(); }                              // all = all computed
  [[nodiscard]] auto getnrkept() const { return getnrpost(); }                                 // # of kept states
  [[nodiscard]] auto getnrdiscarded() const { return getnrcomputed()-getnrpost(); }            // # of discarded states
  [[nodiscard]] auto getnrstored() const  { return value_zero.size(); }                        // number of stored states
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
    value_zero.resize(nrpost);
  }
  // Initialize the data structures with eigenvalues 'v'. The eigenvectors form an identity matrix. This is used to
  // represent the spectral decomposition in the eigenbasis itself.
  void diagonal(const EVEC &v) {
    value_orig = value_zero = v; // YYY
    values.copy(v);
    matrix   = ublas::identity_matrix<t_eigen>(v.size());
  }
  void subtract_Egs(const t_eigen Egs) {
    value_zero = value_orig; // XXX
    for (auto &x : value_zero) x -= Egs; // XXX: subtract a scalar [fix after moving to Eigen]
    my_assert(value_zero[0] >= 0);
  }
  void subtract_GS_energy(const t_eigen GS_energy) {
    for (auto &x : absenergyG) x -= GS_energy; // XXX
    my_assert(absenergyG[0] >= 0);
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
    oa << value_orig;
    values.save(oa);
    NRG::save(oa, matrix);
    oa << value_zero << nrpost << absenergy << absenergyG << absenergy_zero;
  }  
  void load(boost::archive::binary_iarchive &ia) {
    ia >> value_orig;
    values.load(ia);
    NRG::load(ia, matrix);
    ia >> value_zero >> nrpost >> absenergy >> absenergyG >> absenergy_zero;
  }
  void h5save(H5Easy::File &fd, const std::string &name, const bool write_absG) const {
    H5Easy::dump(fd, name + "/value_orig",     value_orig);
    H5Easy::dump(fd, name + "/value_zero",     value_zero);
    H5Easy::dump(fd, name + "/absenergy",      absenergy);
    H5Easy::dump(fd, name + "/absenergy_zero", absenergy_zero);
    if (write_absG) 
      H5Easy::dump(fd, name + "/absenergyG",   absenergyG);
    h5_dump_matrix(fd, name + "/matrix", matrix);
    h5_dump_scalar(fd, name + "/nrkept", getnrkept());
  }
  void dump_eigenvalues(const size_t max_nr = std::numeric_limits<size_t>::max(), std::ostream &F = std::cout) const {
    F << "eig= " << std::setprecision(std::numeric_limits<double>::max_digits10);
    ranges::for_each_n(values.begin(), std::min(getnrcomputed(), max_nr),
                     [&F](const double x) { F << x << ' '; });
    F << std::endl;
  }
  void checkdiag(const double NORMALIZATION_EPSILON = 1e-12, const double ORTHOGONALITY_EPSILON = 1e-12) {
    const auto M   = getnrcomputed(); // number of eigenpairs
    const auto dim = getdim();        // dimension of the eigenvector
    // Check normalization
    for (const auto r : range0(M)) { // ZZZ
      assert_isfinite(values[r]);
      S sumabs{};
      for (const auto j : range0(dim)) sumabs += conj_me(matrix(r, j)) * matrix(r, j);
      my_assert(num_equal(abs(sumabs), 1.0, NORMALIZATION_EPSILON));
    }
    // Check orthogonality
    for (const auto r1 : range0(M))
      for (const auto r2 : boost::irange(r1 + 1, M)) {
        S skpdt{};
        for (const auto j : range0(dim)) skpdt += conj_me(matrix(r1, j)) * matrix(r2, j);
        my_assert(num_equal(abs(skpdt), 0.0, ORTHOGONALITY_EPSILON));
      }
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
   std::vector<t_eigen> sorted_energies() const {
     std::vector<t_eigen> energies;
     for (const auto &eig: eigs())
       energies.insert(energies.end(), eig.value_zero.begin(), eig.value_zero.end());
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
           fmt::print("({}) {} states: {}\n", I.str(), eig.getnrstored(), eig.value_orig);
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
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) { load(N, P, remove_files); } // called from do_diag()
};

} // namespace

#endif
