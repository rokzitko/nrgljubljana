// symmetry.cc - Classes representing various symmetry types
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _symmetry_h_
#define _symmetry_h_

#include "nrg-general.h"

// Check if the triangle inequality is satisfied (i.e. if Clebsch-Gordan coefficient can be different from zero).
// This is important, for example, for triplet operators, which are zero when evaluated between two singlet states.
// Arguments ss1, ss2, ss3 are spin multiplicities. Returns true if the inequality is satisfied, false otherwise.
inline auto su2_triangle_inequality(const int ss1, const int ss2, const int ss3) {
  return (abs(ss1-ss2) <= ss3-1) && (abs(ss2-ss3) <= ss1-1) && (abs(ss3-ss1) <= ss2-1);
}

inline auto u1_equality(const int q1, const int q2, const int q3) { return q1 == q2 + q3; }     // Equality for U(1) symmetry
inline auto z2_equality(const int p1, const int p2, const int p3) { return p1 == p2 * p3; }
inline auto c3_equality(const int p1, const int p2, const int p3) { return p1 == (p2+p3) % 3; } // C_3 quantum number: Equality modulo 3

template<typename S>
inline auto newopch(const Params &P)
{
  Opch<S> opch(P.channels);
  for (auto &oc: opch) {
    oc.resize(P.perchannel);
    for (auto &o: oc)
      o.clear(); // set all ublas matrix elements to zero
  }
  return opch;
}

template<typename S>
struct Recalc_f {
  size_t i1; // subspace indexes
  size_t ip;
  typename traits<S>::t_coef factor;
};

// Structure which holds subspace information and factor for each of nonzero irreducible matrix elements. cf.
// Hofstetter PhD p. 120. <Q+1 S+-1/2 .. i1 ||f^\dag|| Q S .. ip>_N = factor < IN1 .. ||f^\dag|| INp ..>_{N_1}
template<typename S>
struct Recalc {
  size_t i1{}; // combination of states
  size_t ip{};
  Invar IN1; // subspace in N-1 stage
  Invar INp;
  typename traits<S>::t_coef factor{}; // additional multiplicative factor
};

template<typename S>
class Symmetry {
 protected:
   const Params &P;
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
   explicit Symmetry(const Params &P_, const Invar InvarSinglet = {}, const Invar Invar_f = {}) : 
     P(P_), In(P.combs+1), QN(P.combs+1), InvarSinglet(InvarSinglet), Invar_f(Invar_f) {}
   auto input_subspaces() const { return In; }
   auto QN_subspace(const size_t i) const { my_assert(i < P.combs); return QN[i]; }
   auto ancestor(const Invar &I, const size_t i) const {
     my_assert(i < P.combs);
     const auto input = input_subspaces();
     Invar anc = I;
     anc.combine(input[i]);
     return anc; // I.combine(input[i]) == input[i].combine(I)
   }
   size_t nr_combs() const {
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
   // For some symmetry types with two-channels we distinguish between even and odd parity with respect to the
   // channel-interchange operation.
   virtual bool islr() { return false; }
   // Ditto for 3 channels: C_3 symmetry.
   virtual bool isc3() { return false; }
   // For some symmetry types, we may distinguish between spin-up and spin-down quantities (in particular spin-up and
   // spin-down spectral functions).
   virtual bool isfield() { return false; }
   // Multiplicity of the states in the invariant subspace
   virtual size_t mult(const Invar &) const { return 1; };
   auto multfnc() const { return [this](const Invar &I) { return this->mult(I); }; }
   auto calculate_Z(const Invar &I, const Eigen<S> &eig, const double rescale_factor) const {
     return mult(I) * ranges::accumulate(eig.value_zero, 0.0, [rf=rescale_factor](auto sum, const auto &x) { return sum+exp(-rf*x); });
   }
   // Does the combination of subspaces I1 and I2 contribute to the spectral function corresponding to spin SPIN?
   virtual bool check_SPIN(const Invar &I1, const Invar &I2, const int &SPIN) const { return true; }
   // Is the triangle inequality satisfied (i.e. can Clebsch-Gordan coefficient be different from zero)?
   virtual bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const { return true; }
   // Is an invariant subspace with given quantum numbers allowed?
   virtual bool Invar_allowed(const Invar &I) const { return true; }

   using Matrix = typename traits<S>::Matrix;
   using t_matel = typename traits<S>::t_matel;
   using t_coef = typename traits<S>::t_coef;

   void offdiag_function_impl(const Step &step, const size_t i, const size_t j, const size_t ch, const size_t fnr, const t_coef factor,
                              Matrix &h, const Rmaxvals &qq, const InvarVec &In, const Opch<S> &opch) const;
   void diag_function_impl(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta,
                           Matrix &h, const Rmaxvals &qq, const double f) const;
   void diag_function(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta, 
                      Matrix &h, const Rmaxvals &qq) const;
   void diag_function_half(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta,
                           Matrix &h, const Rmaxvals &qq) const;
   void diag_offdiag_function(const Step &step, const size_t i, const size_t j, const size_t chin, const t_coef factor,
                              Matrix &h, const Rmaxvals &qq) const;

   virtual void make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In, 
                            const Opch<S> &opch, const Coef<S> &coef) = 0;

   // Called from recalc_dynamicsusceptibility().  This is the factor due
   // to the spin degeneracy when calculating the trace of Sz.Sz.
   virtual double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   // Called from recalc_dynamic_orb_susceptibility().  This is the factor due
   // to the orbital moment degeneracy when calculating the trace of Tz.Tz.
   virtual double dynamic_orb_susceptibility_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   // Called from calc_specdens().
   // See spectral_density_clebschgordan.nb and DMNRG_clebschgordan.nb.
   virtual double specdens_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }
   virtual double specdensquad_factor(const Invar &Ip, const Invar &I1) const { return 1.0; }

   virtual void calculate_TD(const Step &step, const DiagInfo<S> &diag, const Stats<S> &stats, const double factor) = 0;

   virtual Opch<S> recalc_irreduc(const Step &step, const DiagInfo<S> &diag, const QSrmax &qsrmax) { my_assert_not_reached(); }
   virtual OpchChannel<S> recalc_irreduc_substeps(const Step &step, const DiagInfo<S> &diag, 
                                                       const QSrmax &qsrmax, int M) { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_doublet(const DiagInfo<S> &diag, const QSrmax &qsrmax, 
                                                 const MatrixElements<S> &cold) { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_triplet(const DiagInfo<S> &diag, const QSrmax &qsrmax, 
                                                 const MatrixElements<S> &cold) { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_orb_triplet(const DiagInfo<S> &diag, const QSrmax &qsrmax, 
                                                     const MatrixElements<S> &cold) { my_assert_not_reached(); }
   virtual MatrixElements<S> recalc_quadruplet(const DiagInfo<S> &diag, const QSrmax &qsrmax, 
                                                    const MatrixElements<S> &cold) { my_assert_not_reached(); }
   virtual void recalc_global(const Step &step, const DiagInfo<S> &diag, const QSrmax &qsrmax, std::string name, 
                              MatrixElements<S> &cnew) { my_assert_not_reached(); }

   // Recalculates irreducible matrix elements of a singlet operator, as well as odd-parity spin-singlet operator (for
   //  parity -1). Generic implementation, valid for all symmetry types.
   MatrixElements<S> recalc_singlet(const DiagInfo<S> &diag, const QSrmax &qsrmax, const MatrixElements<S> &nold, const int parity) {
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
       nnew[Twoinvar(I1,Ip)] = recalc_general(diag, qsrmax, nold, I1, Ip, recalc_table, Iop);
     }
     return nnew;
   }

   virtual void show_coefficients(const Step &step, const Coef<S> &coef) {
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

   virtual bool recalc_f_coupled(const Invar &I1, const Invar &I2, const Invar &If) { return true; } // used in recalc_f()

   template<typename T>
     auto recalc_f(const DiagInfo<S> &diag, const QSrmax &qsrmax, const Invar &I1,
                   const Invar &Ip, const T &table);

   template<typename T>
     auto recalc_general(const DiagInfo<S> &diag, const QSrmax &qsrmax, const MatrixElements<S> &cold,
                         const Invar &I1, const Invar &Ip, const T &table, const Invar &Iop) const;

   void recalc1_global(const DiagInfo<S> &diag, const QSrmax &qsrmax, const Invar &I,
                       Matrix &m, const size_t i1, const size_t ip, const t_coef value) const;

   auto CorrelatorFactorFnc() const   { return [this](const Invar &Ip, const Invar &I1) { return this->mult(I1); }; }
   auto SpecdensFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->specdens_factor(Ip, I1); }; }
   auto SpecdensquadFactorFnc() const { return [this](const Invar &Ip, const Invar &I1) { return this->specdensquad_factor(Ip, I1); }; }
   auto SpinSuscFactorFnc() const     { return [this](const Invar &Ip, const Invar &I1) { return this->dynamicsusceptibility_factor(Ip, I1); }; }
   auto OrbSuscFactorFnc() const      { return [this](const Invar &Ip, const Invar &I1) { return this->dynamic_orb_susceptibility_factor(Ip, I1); }; }
   auto TrivialCheckSpinFnc() const   { return [this](const Invar &Ip, const Invar &I1, int SPIN) { return true; }; }
   auto SpecdensCheckSpinFnc() const  { return [this](const Invar &I1, const Invar &Ip, int SPIN) { return this->check_SPIN(I1, Ip, SPIN); }; }
};

// Add DECL declaration in each symmetry class
#define DECL                                                                                                 \
  void make_matrix(Matrix &h, const Step &step, const Rmaxvals &qq, const Invar &I, const InvarVec &In,      \
             const Opch<SC> &opch, const Coef<SC> &coef) override;                                           \
  Opch<SC> recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax) override

// Optional declaration
#define HAS_SUBSTEPS OpchChannel<SC> recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax, int M) override
#define HAS_DOUBLET MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, \
                                                      const MatrixElements<SC> &cold) override
#define HAS_TRIPLET MatrixElements<SC> recalc_triplet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, \
                                                      const MatrixElements<SC> &cold) override
#define HAS_ORB_TRIPLET MatrixElements<SC> recalc_orb_triplet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, \
                                                              const MatrixElements<SC> &cold) override
#define HAS_QUADRUPLET MatrixElements<SC> recalc_quadruplet(const DiagInfo<SC> &diag, const QSrmax &qsrmax, \
                                                            const MatrixElements<SC> &cold) override
#define HAS_GLOBAL void recalc_global(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax, \
                                      std::string name, MatrixElements<SC> &cnew) override

template<typename S>
class SymField : public Symmetry<S> {
 public:
   template<typename ... Args> explicit SymField(Args && ... args) : Symmetry<S>(std::forward<Args>(args)...) {}
   bool isfield() override { return true; }
};

template<typename S>
class SymLR : public Symmetry<S> {
 public:
   template<typename ... Args> explicit SymLR(Args && ... args) : Symmetry<S>(std::forward<Args>(args)...) {}
   bool islr() override { return true; }
};

template<typename S>
class SymC3 : public Symmetry<S> {
 public:
   template<typename ... Args> explicit SymC3(Args && ... args) : Symmetry<S>(std::forward<Args>(args)...) {}
   bool isc3() override { return true; }
};

template<typename S>
class SymFieldLR : public Symmetry<S> {
  public:
   template<typename ... Args> explicit SymFieldLR(Args && ... args) : Symmetry<S>(std::forward<Args>(args)...) {}
   bool isfield() override { return true; }
   bool islr() override { return true; }
};

// Helper functions
inline void check_abs_diff(const Invar &Ip, const Invar &I1, const std::string &what, int diff) {
  const auto a = Ip.get(what);
  const auto b = I1.get(what);
  my_assert(abs(b - a) == diff);
}

inline void check_diff(const Invar &Ip, const Invar &I1, const std::string &what, int diff) {
  const auto a = Ip.get(what);
  const auto b = I1.get(what);
  my_assert(b - a == diff);
}

#endif
