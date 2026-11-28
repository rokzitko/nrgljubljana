// symmetry.cc - Classes representing various symmetry types
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _symmetry_hpp_
#define _symmetry_hpp_

#include <cmath>
#include <complex>
#include <cstddef>
#include <functional>
#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <limits>
#include <optional>
#include <ostream>
#include <utility>

#include "operators.hpp"
#include "params.hpp"
#include "traits.hpp"
#include "operators.hpp"
#include "invar.hpp"
#include "step.hpp"
#include "eigen.hpp"
#include "subspaces.hpp"
#include "stats.hpp"
#include "coef.hpp"
#include "debug.hpp"

namespace NRG {

using std::string_literals::operator""s;

using cmpl = std::complex<double>;

// Check if the triangle inequality is satisfied (i.e. if Clebsch-Gordan coefficient can be different from zero).
// This is important, for example, for triplet operators, which are zero when evaluated between two singlet states.
// Arguments ss1, ss2, ss3 are spin multiplicities. Returns true if the inequality is satisfied, false otherwise.
inline auto su2_triangle_inequality(const int ss1, const int ss2, const int ss3) {
  return (std::abs(ss1-ss2) <= ss3-1) && (std::abs(ss2-ss3) <= ss1-1) && (std::abs(ss3-ss1) <= ss2-1);
}

inline auto u1_equality(const int q1, const int q2, const int q3) { return q1 == q2 + q3; }     // Equality for U(1) symmetry
inline auto z2_equality(const int p1, const int p2, const int p3) { return p1 == p2 * p3; }
inline auto c3_equality(const int p1, const int p2, const int p3) { return p1 == (p2+p3) % 3; } // C_3 quantum number: Equality modulo 3

// Factors due to the SU(2) multiplet structure in the traces over the invariant subspaces. Arguments are the
// multiplicities (2j+1) of the subspaces Ip and I1. See spectral_density_clebschgordan.nb and DMNRG_clebschgordan.nb.
// Doublet operator (spectral functions), j1 = jp +- 1/2:
inline double su2_doublet_factor(const int mp, const int m1) {
  my_assert(std::abs(m1 - mp) == 1);
  return m1 == mp + 1 ? S(mp) + 1.0 : S(mp);
}
// Triplet operator (spin susceptibilities), j1 = jp, jp +- 1:
inline double su2_triplet_factor(const int mp, const int m1) {
  my_assert(std::abs(m1 - mp) == 2 || m1 == mp);
  return switch3(m1, mp + 2, 1. + (mp - 1) / 3., mp, mp / 3., mp - 2, (-2. + mp) / 3.);
}
// Orbital-triplet operator (orbital angular momentum T, multiplicities 2T+1): as su2_triplet_factor(), except that two
// orbital singlets do not couple.
inline double orbital_triplet_factor(const int mp, const int m1) { return m1 == 1 && mp == 1 ? 0.0 : su2_triplet_factor(mp, m1); }

// Doublet operator with respect to one of two SU(2) labels qn1, qn2 (multiplicities); the other label is a spectator
// and contributes its multiplicity.
inline double two_su2_doublet_factor(const Invar &Ip, const Invar &I1, const std::string &qn1, const std::string &qn2) {
  const int m1p = Ip.get(qn1), m11 = I1.get(qn1), m2p = Ip.get(qn2), m21 = I1.get(qn2);
  if (std::abs(m11 - m1p) == 1) {
    my_assert(m2p == m21 && m2p >= 1);
    return su2_doublet_factor(m1p, m11) * m2p;
  }
  my_assert(m1p == m11 && m1p >= 1);
  return su2_doublet_factor(m2p, m21) * m1p;
}

// A thermodynamic quantity computed by Symmetry::calculate_TD() as a thermal average over the invariant subspaces:
// the column 'name' of the file td, and the value of the quantity in a subspace, averaged over the states of an SU(2)
// multiplet. With variance=true, the column holds <x^2> - <x>^2 instead of <x>.
struct TDTrace {
  std::string name;
  std::function<double(const Invar &)> value;
  bool variance = false;
};

namespace td {
// Abelian quantum number qn: its mean value <qn>, the mean of its square <qn^2>, and its variance
inline TDTrace mean(std::string name, std::string qn) { return {std::move(name), [qn](const Invar &I) { return double(I.get(qn)); }}; }
inline TDTrace square(std::string name, std::string qn) {
  return {std::move(name), [qn](const Invar &I) { const auto q = I.get(qn); return double(q * q); }};
}
inline TDTrace variance(std::string name, std::string qn) { return {std::move(name), [qn](const Invar &I) { return double(I.get(qn)); }, true}; }
// Spin projection, stored as 2S_z+1 (label SSZ): <S_z> and <S_z^2>
inline TDTrace sz(std::string name, std::string qn = "SSZ") { return {std::move(name), [qn](const Invar &I) { return SZ(I.get(qn)); }}; }
inline TDTrace sz_square(std::string name, std::string qn = "SSZ") {
  return {std::move(name), [qn](const Invar &I) { const auto sz = SZ(I.get(qn)); return sz * sz; }};
}
// SU(2) multiplet with multiplicity m = 2J+1 (label qn): <J_z^2> = J(J+1)/3 = (m^2-1)/12, times 'scale'
inline TDTrace su2_z_square(std::string name, std::string qn, const double scale = 1.0) {
  return {std::move(name), [qn, scale](const Invar &I) { const auto m = I.get(qn); return scale * (m * m - 1) / 12.; }};
}
} // namespace td

template<scalar S>
struct Recalc_f {
  size_t i1; // subspace indexes
  size_t ip;
  coef_traits<S> factor;
};

// Structure which holds subspace information and factor for each of nonzero irreducible matrix elements. cf.
// Hofstetter PhD p. 120. <Q+1 S+-1/2 .. i1 ||f^\dag|| Q S .. ip>_N = factor < IN1 .. ||f^\dag|| INp ..>_{N_1}
template<scalar S>
struct Recalc {
  size_t i1{}; // combination of states
  size_t ip{};
  Invar IN1; // subspace in N-1 stage
  Invar INp;
  coef_traits<S> factor{}; // additional multiplicative factor
};

template<scalar S> using RecalcTable = std::vector<Recalc<S>>;
template<scalar S> using RecalcFTable = std::vector<Recalc_f<S>>;

template<scalar S>
class Symmetry {
 protected:
   const Params &P;
   const std::vector<TDTrace> td_traces; // thermodynamic quantities in td, see calculate_TD()
   // In and QN contain information about how the invariant subspaces at consecutive iteration steps are combined. In
   // is the array of quantum number DIFFERENCES used in the construction of the basis. QN is the array of the
   // conserved quantum numbers corresponding to the states being added. For example, in case of SU(2) symmetry, In
   // will include S_z, while QN will include S. In other words, QN involves those quantum numbers that we need to
   // retain in the calculation, while In involves those quantum numbers that "drop out" of the problem due to the
   // symmetry.
   std::vector<Invar> In, QN;
   const Invar InvarSinglet; // QNs for singlet operator
   const Invar Invar_f;      // QNs for f operator
 public:
   virtual void load() = 0; // load In, QN
   void erase_first() { // drop the first element in In, QN to convert to 0-based vectors; call after load()
     In.erase(In.begin());
     QN.erase(QN.begin());
   }
   Symmetry(const Params &P_, std::vector<TDTrace> td_traces_, const Invar & InvarSinglet_ = {}, const Invar & Invar_f_ = {}) :
     P(P_), td_traces(std::move(td_traces_)), In(P_.combs+1), QN(P_.combs+1), InvarSinglet(InvarSinglet_), Invar_f(Invar_f_) {}
   Symmetry(const Symmetry &) = delete;
   Symmetry(Symmetry &&) = delete;
   Symmetry &operator=(const Symmetry &) = delete;
   Symmetry &operator=(Symmetry &&) = delete;
   virtual ~Symmetry() {}
   [[nodiscard]] const auto &input_subspaces() const { return In; }
   auto QN_subspace(const size_t i) const { my_assert(i < P.combs); return QN[i]; }
   auto ancestor(const Invar &I, const size_t i) const {
     my_assert(i < P.combs);
     Invar anc = I;
     anc.combine(In[i]);
     return anc; // I.combine(In[i]) == In[i].combine(I)
   }
   [[nodiscard]] size_t nr_combs() const {
     my_assert(P.combs == In.size());
     my_assert(P.combs == QN.size());
     return P.combs;
   }
   auto combs() const { return range0(nr_combs()); }
   auto ancestors(const Invar &I) const {
     auto input = input_subspaces();
     for (const auto i: combs())
       input[i].combine(I); // In is the list of differences wrt I
     return input;
   }
   auto new_subspaces(const Invar &I) const {
     auto input = input_subspaces();
     for (const auto i: combs()) {
       input[i].inverse();
       input[i].combine(I);
     }
     return input;
   }
   auto get_td_fields() const {
     std::vector<std::string> names;
     for (const auto &t : td_traces) names.push_back(t.name);
     return names;
   }
   // For some symmetry types with two-channels we distinguish between even and odd parity with respect to the
   // channel-interchange operation.
   virtual bool islr() const { return false; }
   // Ditto for 3 channels: C_3 symmetry.
   virtual bool isc3() const { return false; }
   // For some symmetry types, we may distinguish between spin-up and spin-down quantities (in particular spin-up and
   // spin-down spectral functions).
   virtual bool isfield() const { return false; }
   // Superconducting symmetry types include the on-site pairing term (Wilson-chain coefficients delta) and possibly
   // the anomalous hopping term (coefficients kappa) of the added site in the Hamiltonian. All other symmetry types
   // ignore these coefficients, so nonzero values are rejected when the input is read (see read-input.hpp).
   //
   // Phase convention: the ISOSPINX tables hold the (real) matrix elements of the isospin x component between shell
   // states i < j, and only the upper triangle of the Hamiltonian is filled. Since the state with the lower index has
   // fewer electrons, the coefficient of these entries is conj(delta). This is consistent with the pairing term of the
   // initial sites set up by nrginit and is the same for all superconducting symmetry types (test/c++sym_conventions)
   // except SPSU2T, which supports only real delta (complex_pairing_supported() is false).
   [[nodiscard]] virtual bool uses_pairing() const { return false; }
   [[nodiscard]] virtual bool uses_anomalous_hopping() const { return false; }
   [[nodiscard]] virtual bool complex_pairing_supported() const { return true; }
   // Multiplicity of the states in the invariant subspace
   [[nodiscard]] virtual size_t mult(const Invar &) const { return 1; };
   auto multfnc() const { return [this](const Invar &I) { return this->mult(I); }; }
   auto calculate_Z(const Invar &I, const Eigen<S> &eig, const double rescale_factor) const {
     return mult(I) * ranges::accumulate(eig.value_corr_msr(), 0.0, {}, [rf=rescale_factor](const auto &x) { return std::exp(-rf*x); });
   }
   // Does the combination of subspaces I1 and I2 contribute to the spectral function corresponding to spin SPIN?
   [[nodiscard]] virtual bool check_SPIN([[maybe_unused]] const Invar &I1, [[maybe_unused]] const Invar &I2, [[maybe_unused]] const int &SPIN) const { return true; }
   // Is the triangle inequality satisfied (i.e. can Clebsch-Gordan coefficient be different from zero)?
   [[nodiscard]] virtual bool triangle_inequality([[maybe_unused]] const Invar &I1, [[maybe_unused]] const Invar &I2, [[maybe_unused]] const Invar &I3) const { return true; }
   // Is an invariant subspace with given quantum numbers allowed?
   [[nodiscard]] virtual bool Invar_allowed([[maybe_unused]] const Invar &I) const { return true; }

   // Project the states before taking measurements. String p defines what kind of projection is performed.
   virtual bool project_subspace([[maybe_unused]] const Invar &I, [[maybe_unused]] const std::string &p) const {
     return true; // by default retain all (i.e., no projection)
   }
   virtual DiagInfo<S> project(const DiagInfo<S> &diag, const std::string &p) const {
     my_assert(p != ""s); // for performance reasons, this should never be called if the projection is not specified explicitly
     DiagInfo<S> proj;
     if (p == "trivial"s) {
       proj = diag; // no projection
     } else {
       for (const auto &[I, eig]: diag) if (project_subspace(I, p)) proj[I] = eig;
     }
     return proj;
   }

   using Matrix  = Matrix_traits<S>;
   using t_matel = matel_traits<S>;
   using t_coef  = coef_traits<S>;

   void offdiag_function_impl(const Step &step, const size_t i, const size_t j, const size_t ch, const size_t fnr, const t_coef factor,
                              Matrix &h, const SubspaceDimensions &qq, const InvarVec &In, const Opch<S> &opch) const;
   void diag_function_impl(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                           Matrix &h, const SubspaceDimensions &qq, const double f) const;
   void diag_function(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                      Matrix &h, const SubspaceDimensions &qq) const;
   void diag_function_half(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                           Matrix &h, const SubspaceDimensions &qq) const;
   void diag_offdiag_function(const Step &step, const size_t i, const size_t j, const t_coef factor,
                              Matrix &h, const SubspaceDimensions &qq) const;

   virtual void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In,
                            const Opch<S> &opch, const Coef<S> &coef) const = 0;

   // Called from recalc_dynamicsusceptibility().  This is the factor due
   // to the spin degeneracy when calculating the trace of Sz.Sz.
   [[nodiscard]] virtual double dynamicsusceptibility_factor([[maybe_unused]] const Invar &Ip, [[maybe_unused]] const Invar &I1) const { return 1.0; }

   // Called from recalc_dynamic_orb_susceptibility().  This is the factor due
   // to the orbital moment degeneracy when calculating the trace of Tz.Tz.
   [[nodiscard]] virtual double dynamic_orb_susceptibility_factor([[maybe_unused]] const Invar &Ip, [[maybe_unused]] const Invar &I1) const { return 1.0; }

   // Called from calc_specdens().
   // See spectral_density_clebschgordan.nb and DMNRG_clebschgordan.nb.
   [[nodiscard]] virtual double specdens_factor([[maybe_unused]] const Invar &Ip, [[maybe_unused]] const Invar &I1) const { return 1.0; }
   [[nodiscard]] virtual double specdensquad_factor([[maybe_unused]]  const Invar &Ip, [[maybe_unused]]  const Invar &I1) const { return 1.0; }

   // Thermal averages of the quantities listed in td_traces. 'factor' rescales the energies in the Boltzmann weights.
   virtual void calculate_TD([[maybe_unused]] const Step &step, const DiagInfo<S> &diag, Stats<S> &stats, const double factor) const {
     std::vector<bucket> tr(td_traces.size()), tr2(td_traces.size()); // Tr[x], Tr[x^2]
     for (const auto &[I, eig] : diag) {
       const double sumZ = calculate_Z(I, eig, factor);
       for (size_t k = 0; k < td_traces.size(); k++) {
         const auto x = td_traces[k].value(I);
         tr[k] += sumZ * x;
         if (td_traces[k].variance) tr2[k] += sumZ * x * x;
       }
     }
     for (size_t k = 0; k < td_traces.size(); k++) {
       const double mean = tr[k] / stats.Z;
       stats.td.set(td_traces[k].name, td_traces[k].variance ? tr2[k] / stats.Z - mean * mean : mean);
     }
   }

   virtual Opch<S> recalc_irreduc([[maybe_unused]] const Step &step, [[maybe_unused]] const DiagInfo<S> &diag) const { my_assert_not_reached(); }
   virtual OpchChannel<S> recalc_irreduc_substeps([[maybe_unused]] const Step &step, [[maybe_unused]] const DiagInfo<S> &diag,
                                                       [[maybe_unused]] int M) const { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_doublet([[maybe_unused]] const DiagInfo<S> &diag, [[maybe_unused]] const SubspaceStructure &substruct,
                                                 [[maybe_unused]] const MatrixElements<S> &cold) const { my_assert_not_reached(); }
   // Optional restriction of the doublet recalculation to the transfers of quantum numbers present in the seed of an
   // operator (parameter doublet_transfer_mask, see doublet_transfer.hpp): the mask of a seed (nullopt if the symmetry
   // type does not support the restriction), and the recalculation restricted to the transfers in a mask.
   [[nodiscard]] virtual std::optional<unsigned> doublet_transfer_mask([[maybe_unused]] const std::string &name,
                                                                       [[maybe_unused]] const MatrixElements<S> &seed) const { return std::nullopt; }
   virtual MatrixElements<S> recalc_doublet_masked(const DiagInfo<S> &diag, const SubspaceStructure &substruct,
                                                   const MatrixElements<S> &cold, [[maybe_unused]] unsigned mask) const {
     return recalc_doublet(diag, substruct, cold);
   }
   virtual MatrixElements<S> recalc_triplet([[maybe_unused]] const DiagInfo<S> &diag, [[maybe_unused]] const SubspaceStructure &substruct,
                                                 [[maybe_unused]] const MatrixElements<S> &cold) const { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_orb_triplet([[maybe_unused]] const DiagInfo<S> &diag, [[maybe_unused]] const SubspaceStructure &substruct,
                                                     [[maybe_unused]] const MatrixElements<S> &cold) const { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_quadruplet([[maybe_unused]] const DiagInfo<S> &diag, [[maybe_unused]] const SubspaceStructure &substruct,
                                                    [[maybe_unused]] const MatrixElements<S> &cold) const { my_assert_not_reached(); }
   virtual void recalc_global([[maybe_unused]] const Step &step, [[maybe_unused]] const DiagInfo<S> &diag, [[maybe_unused]] std::string name,
                              [[maybe_unused]] MatrixElements<S> &cnew) const { my_assert_not_reached(); }

   // Recalculates irreducible matrix elements of a singlet operator, as well as odd-parity spin-singlet operator (for
   //  parity -1). Generic implementation, valid for all symmetry types.
   MatrixElements<S> recalc_singlet(const DiagInfo<S> &diag, const SubspaceStructure &substruct, const MatrixElements<S> &nold, const int parity) const {
     MatrixElements<S> nnew;
     my_assert(islr() ? parity == 1 || parity == -1 : parity == 1);
     for (const auto &I : diag.subspaces()) {
       const Invar I1 = I;
       const Invar Ip = parity == -1 ? I.InvertParity() : I;
       std::vector<Recalc<S>> recalc_table;
       for (const auto i: combs()) {
         const auto anc = ancestor(I, i);
         recalc_table.push_back({i+1, i+1, anc, parity == -1 ? anc.InvertParity() : anc, 1.0});
       }
       const auto Iop = parity == -1 ? InvarSinglet.InvertParity() : InvarSinglet;
       auto nn = recalc_general(diag, substruct, nold, I1, Ip, recalc_table, Iop);
       if (nn) nnew[Twoinvar(I1,Ip)] = *nn;
     }
     return nnew;
   }

   virtual void show_coefficients(const Step &step, const Coef<S> &coef) const {
     if (P.silent) return;
     std::cout << std::setprecision(std::numeric_limits<double>::max_digits10);
     if (!P.substeps) {
       for (size_t i = 0; i < P.coefchannels; i++) {
         const auto N = step.N();
         std::cout << "[" << i + 1 << "]"
           << " xi(" << N << ")=" << coef.xi(N, i) << " xi_scaled(" << N << ")=" << coef.xi(N, i)/step.scale()
             << " zeta(" << N+1 << ")=" << coef.zeta(N+1, i) << std::endl;
       }
     } else {
       const auto [N, M] = step.NM();
       for (auto i = 0; i < P.coeffactor; i++) {
         const auto index = M + P.channels * i;
         std::cout << "[" << index << "]"
           << " xi(" << N << ")=" << coef.xi(N, index) << " zeta(" << N+1 << ")=" << coef.zeta(N+1, index) << std::endl;
       }
     }
   }

   // only overriden for symtypes QST and SPSU2T
   virtual bool recalc_f_coupled([[maybe_unused]] const Invar &I1, [[maybe_unused]] const Invar &I2, [[maybe_unused]] const Invar &If) const { return true; } // used in recalc_f()

   [[nodiscard]] bool have_stored_subspaces(const DiagInfo<S> &diag, const Invar &I1, const Invar &Ip) const {
     const auto it1 = diag.find(I1);
     if (it1 == diag.end() || !it1->second.getnrstored()) return false;
     const auto itp = diag.find(Ip);
     return itp != diag.end() && itp->second.getnrstored();
   }

   [[nodiscard]] bool have_stored_coupled_subspaces(const DiagInfo<S> &diag, const Invar &I1, const Invar &Ip) const {
     return recalc_f_coupled(I1, Ip, Invar_f) && have_stored_subspaces(diag, I1, Ip);
   }

   template<typename T>
      auto recalc_f(const DiagInfo<S> &diag, const Invar &I1, const Invar &Ip, const T &table) const;

   // One recalculation table of coefnew (rows {i1, ip, IN1, INp, factor}, see Recalc) for the matrix element
   // <I1||O||Ip> of an operator with quantum numbers Iop, stored in cnew. 'table' returns the rows; it is only called
   // if both subspaces have stored states, since the factors may be undefined otherwise.
   template<typename F>
   void recalc_table(const DiagInfo<S> &diag, const SubspaceStructure &substruct, const MatrixElements<S> &cold, const Invar &I1,
                     const Invar &Ip, const Invar &Iop, MatrixElements<S> &cnew, F &&table) const {
     nrglog('f', "RECALC(I1=" << I1 << ", Ip=" << Ip << ", Iop=" << Iop << ")");
     if (!have_stored_subspaces(diag, I1, Ip)) return;
     if (auto cn = recalc_general(diag, substruct, cold, I1, Ip, table(), Iop)) cnew[Twoinvar(I1, Ip)] = std::move(*cn);
   }
   // Ditto for the irreducible matrix elements <I1||f||Ip> (rows {i1, ip, factor}, see Recalc_f), stored in f
   template<typename F>
   void recalc_f_table(const DiagInfo<S> &diag, const Invar &I1, const Invar &Ip, MatrixElements<S> &f, F &&table) const {
     nrglog('f', "RECALC_F(I1=" << I1 << ", Ip=" << Ip << ")");
     if (have_stored_coupled_subspaces(diag, I1, Ip)) f[Twoinvar(I1, Ip)] = recalc_f(diag, I1, Ip, table());
   }

   template<typename T>
     std::optional<Matrix_traits<S>> recalc_general(const DiagInfo<S> &diag, const SubspaceStructure &substruct, const MatrixElements<S> &cold,
                         const Invar &I1, const Invar &Ip, const T &table, const Invar &Iop) const;

   void recalc1_global(const DiagInfo<S> &diag, const Invar &I,
                       Matrix &m, const size_t i1, const size_t ip, const t_coef value) const;

   auto CorrelatorFactorFnc() const   { return [this]([[maybe_unused]] const Invar &Ip, const Invar &I1) { return this->mult(I1); }; }
   auto SpecdensFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->specdens_factor(Ip, I1); }; }
   auto SpecdensquadFactorFnc() const { return [this](const Invar &Ip, const Invar &I1) { return this->specdensquad_factor(Ip, I1); }; }
   auto SpinSuscFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->dynamicsusceptibility_factor(Ip, I1); }; }
   auto OrbSuscFactorFnc() const      { return [this](const Invar &Ip, const Invar &I1) { return this->dynamic_orb_susceptibility_factor(Ip, I1); }; }
   auto TrivialCheckSpinFnc() const   { return []([[maybe_unused]] const Invar &Ip, [[maybe_unused]] const Invar &I1, [[maybe_unused]] int SPIN) { return true; }; }
   auto SpecdensCheckSpinFnc() const  { return [this](const Invar &I1, const Invar &Ip, int SPIN) { return this->check_SPIN(I1, Ip, SPIN); }; }
};

// Add DECL declaration in each symmetry class
#define DECL                                                                                                      \
  void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &I, const InvarVec &In, \
             const Opch<SC> &opch, const Coef<SC> &coef) const override;                                          \
  Opch<SC> recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const override

// Optional declaration
#define HAS_SUBSTEPS OpchChannel<SC> recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, int M) const override
#define HAS_DOUBLET MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, \
                                                      const SubspaceStructure &substruct, \
                                                      const MatrixElements<SC> &cold) const override
#define HAS_TRIPLET MatrixElements<SC> recalc_triplet(const DiagInfo<SC> &diag, \
                                                      const SubspaceStructure &substruct, \
                                                      const MatrixElements<SC> &cold) const override
#define HAS_ORB_TRIPLET MatrixElements<SC> recalc_orb_triplet(const DiagInfo<SC> &diag, \
                                                              const SubspaceStructure &substruct, \
                                                              const MatrixElements<SC> &cold) const override
#define HAS_QUADRUPLET MatrixElements<SC> recalc_quadruplet(const DiagInfo<SC> &diag, \
                                                            const SubspaceStructure &substruct, \
                                                            const MatrixElements<SC> &cold) const override
#define HAS_GLOBAL void recalc_global(const Step &step, const DiagInfo<SC> &diag, \
                                      std::string name, MatrixElements<SC> &cnew) const override

// The spin projection of a doublet operator is the difference of the labels SSZ (2S_z+1) of the two subspaces.
inline bool sz_check_spin(const Invar &I1, const Invar &Ip, const int SPIN) { return I1.get("SSZ") - Ip.get("SSZ") == SPIN; }

template<scalar S>
class SymField : public Symmetry<S> {
 public:
   using Symmetry<S>::Symmetry;
   bool isfield() const final override { return true; }
   bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override { return sz_check_spin(I1, Ip, SPIN); }
};

template<scalar S>
class SymLR : public Symmetry<S> {
 public:
   using Symmetry<S>::Symmetry;
   bool islr() const final override { return true; }
};

template<scalar S>
class SymC3 : public Symmetry<S> {
 public:
   using Symmetry<S>::Symmetry;
   bool isc3() const final override { return true; }
};

template<scalar S>
class SymFieldLR : public Symmetry<S> {
  public:
   using Symmetry<S>::Symmetry;
   bool isfield() const final override { return true; }
   bool islr() const final override { return true; }
   bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override { return sz_check_spin(I1, Ip, SPIN); }
};

// Helper functions
inline void check_abs_diff(const Invar &Ip, const Invar &I1, const std::string &what, const int diff) {
  const auto a = Ip.get(what);
  const auto b = I1.get(what);
  my_assert(std::abs(b - a) == diff);
}

inline void check_diff(const Invar &Ip, const Invar &I1, const std::string &what, const int diff) {
  const auto a = Ip.get(what);
  const auto b = I1.get(what);
  my_assert(b - a == diff);
}

} // namespace

#endif
