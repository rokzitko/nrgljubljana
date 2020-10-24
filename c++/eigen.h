#ifndef _eigen_h_
#define _eigen_h_

#include <vector>
#include <string>
#include "portabil.h"
#include "traits.h"
#include "invar.h"
#include "params.h"

// Result of a diagonalisation: eigenvalues and eigenvectorse
template <typename S> class Eigen {
public:
  using t_eigen = typename traits<S>::t_eigen;
  using EVEC = ublas::vector<t_eigen>;
  using Matrix = typename traits<S>::Matrix;
  EVEC value_orig; // eigenvalues as computed
  Matrix matrix;   // eigenvectors
  Eigen() {}
  Eigen(const size_t nr, const size_t dim) {
    my_assert(nr <= dim);
    value_orig.resize(nr);
    matrix.resize(nr, dim);
  }
  auto getnrcomputed() const { return value_orig.size(); } // number of computed eigenpairs
  auto getdim() const { return matrix.size2(); }           // valid also after the split_in_blocks_Eigen() call
 private:
  long nrpost = -1;  // number of eigenpairs after truncation (-1: keep all)
 public:
  EVEC value_zero;                                                               // eigenvalues with Egs subtracted
  auto getnrpost() const { return nrpost == -1 ? getnrcomputed() : nrpost; }     // number of states after truncation
  auto getnrstored() const  { return value_zero.size(); }                        // number of stored states
  auto getnrall() const { return getnrcomputed(); }                              // all = all computed
  auto getnrkept() const { return getnrpost(); }                                 // # of kept states
  auto getnrdiscarded() const { return getnrcomputed()-getnrpost(); }            // # of discarded states
  auto all() const { return range0(getnrcomputed()); }                           // iterator over all states
  auto kept() const { return range0(getnrpost()); }                              // iterator over kept states
  auto discarded() const { return boost::irange(getnrpost(), getnrcomputed()); } // iterator over discarded states
  auto stored() const { return range0(getnrstored()); }                          // iterator over all stored states
  // NOTE: "absolute" energy means that it is expressed in the absolute energy scale rather than SCALE(N).
  EVEC absenergy;      // absolute energies
  EVEC absenergyG;     // absolute energies (0 is the absolute ground state of the system) [SAVED TO FILE]
  EVEC absenergyN;     // absolute energies (referenced to the lowest energy in the N-th step)
  // 'blocks' contains eigenvectors separated according to the invariant subspace from which they originate.
  // Required for using efficient BLAS routines when performing recalculations of the matrix elements.
  std::vector<Matrix> blocks;
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
  auto diagonal_exp(const double factor) const { // produce a diagonal matrix with exp(-factor*E) diagonal elements
    const auto dim = getnrstored();
    auto m = Zero_matrix<S>(dim);
    for (const auto i: range0(dim)) 
      m(i, i) = exp(-value_zero(i) * factor);
    return m;
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
template<typename S>
class DiagInfo : public std::map<Invar, Eigen<S>> {
 public:
   using t_eigen = typename traits<S>::t_eigen;
   using Matrix  = typename traits<S>::Matrix;
   explicit DiagInfo() {}
   DiagInfo(std::ifstream &fdata, const size_t nsubs, const Params &P) {
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
   auto size_subspace(const Invar &I) const {
     const auto f = this->find(I);
     return f != this->cend() ? f->second.getnrstored() : 0;
   }
   auto subs(const Invar &I1, const Invar &I2) const {
     return std::make_pair(this->at(I1), this->at(I2));
   }
   auto dims(const Invar &I1, const Invar &I2) const { // Determines matrix sizes for operators (# stored)
     return std::make_pair(size_subspace(I1), size_subspace(I2));
   }
   void clear_eigenvectors() {
     for (auto &eig : this->eigs())
       for (auto &m : eig.blocks) 
         m = Matrix();
   }
   // Total number of states (symmetry taken into account)
   template <typename MF> auto count_states(MF && mult) const {
     return ranges::accumulate(*this, 0, [mult](auto n, const auto &x) { const auto &[I, eig] = x; return n + mult(I)*eig.getnrstored(); });
   }
   auto count_subspaces() const {    // Count non-empty subspaces
     return ranges::count_if(this->eigs(), [](const auto &eig) { return eig.getnrstored()>0; });
   }
   template<typename F, typename M> auto trace(F fnc, const double factor, M mult) const { // Tr[fnc exp(-factor*E)]
     auto b = 0.0;
     for (const auto &[I, eig] : *this)
       b += mult(I) * ranges::accumulate(eig.value_zero, 0.0, [fnc, factor](auto acc, const auto x) { 
         const auto betaE = factor * x; return acc + fnc(betaE) * exp(-betaE); });
     return b;
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
     const std::string fn = P.workdir.unitaryfn(N);
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
     const std::string fn = P.workdir.unitaryfn(N);
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
   explicit DiagInfo(const size_t N, const Params &P, const bool remove_files = false) { load(N, P, remove_files); } // called from do_diag()
};

#endif
