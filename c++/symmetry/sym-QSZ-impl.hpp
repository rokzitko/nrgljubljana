#include <array>
#include <cmath>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

#include "abelian-common.hpp"

namespace NRG {

using std::pow;

// U(1) charge and U(1) spin symmetry (conserved Q and S_z). Any number of channels; the coefficient tables are
// generated at run time (abelian_shell.hpp). One <||f||> block per channel holds both f^dag_up and f^dag_down.
//
// The class also implements QSZ4 (sym-QSZ4-impl.hpp), where only 2S_z modulo 4 is conserved: with modulus = 4 the
// spin label SSZ is taken modulo 4 (see SymmetrySPU1 for the same construction).
template<typename SC>
class SymmetryQSZ : public SymField<SC> {
 protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
   using t_coef = coef_traits<SC>;

   AbelianShell::Kind kind = AbelianShell::Kind::QSZ;
   int modulus = 0; // 0: SSZ is additive (QSZ); 4: SSZ is taken modulo 4 (QSZ4)

   [[nodiscard]] Invar qssz(const int q, const int x) const { return Invar(q, modulus ? ((x % modulus) + modulus) % modulus : x); }

   std::optional<AbelianShell::Shell> shell;
   // Hopping f^dag_N a_{N+1} for both spins and for each spin separately (polarized chains). All (i,j), see
   // Symmetry::offdiag_function_impl().
   std::vector<std::vector<AbelianShell::Entry<double>>> offdiag;
   std::vector<std::array<std::vector<AbelianShell::Entry<double>>, 2>> offdiag_spin;
   std::vector<std::array<std::vector<Recalc_f<SC>>, 2>> ftables; // [ch][spin]: f^dag_sigma
   // Channel mixing for P.rungs (two channels): f^dag_{N,ch} a_{N+1,1-ch}, and on-site hopping between the channels
   std::vector<std::vector<AbelianShell::Entry<double>>> offdiag_mix;
   std::vector<AbelianShell::Entry<AbelianShell::cmpl>> runghop;

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryQSZ(const Params &P, std::vector<TDTrace> td_traces) : SymField<SC>(P, std::move(td_traces), Invar(0,0), Invar(1,2)) {}
   SymmetryQSZ(const Params &P) : SymmetryQSZ(P, std::vector<TDTrace>{td::sz_square("<Sz^2>"), td::sz("<Sz>"), td::mean("<Q>", "Q"), td::square("<Q^2>", "Q")}) {
     initInvar({
       {"Q", additive},  // charge
       {"SSZ", additive} // spin projection
     });
   }

   // The spin projection of the operator is defined by the difference in Sz of both the invariant subspaces.
   bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
     int sszop = I1.get("SSZ") - Ip.get("SSZ");
     if (!modulus) return sszop == SPIN;
     sszop = ((sszop % modulus) + modulus) % modulus;
     return (sszop == 1 && SPIN == 1) || (sszop == 3 && SPIN == -1);
   }

   bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override {
     auto I = I2;
     I.combine(I3);
     return I == I1;
   }

   void load() override {
     using namespace AbelianShell;
     shell.emplace(kind, shell_channels(P));
     load_shell(*shell, In, QN, P);
     const auto K = shell->channels();
     offdiag.assign(K, {});
     offdiag_spin.assign(K, {});
     ftables.assign(K, {});
     for (std::size_t ch = 0; ch < K; ++ch) {
       offdiag[ch] = shell->hop_cr(c(ch, Spin::up) + c(ch, Spin::down));
       for (const auto s : both_spins) {
         offdiag_spin[ch][static_cast<std::size_t>(s)] = shell->hop_cr(c(ch, s));
         ftables[ch][static_cast<std::size_t>(s)] = f_table<SC>(*shell, cdag(ch, s));
       }
     }
     offdiag_mix.clear();
     runghop.clear();
     if (K == 2) {
       for (std::size_t ch = 0; ch < 2; ++ch) offdiag_mix.push_back(shell->hop_cr(c(1 - ch, Spin::up) + c(1 - ch, Spin::down)));
       runghop = shell->entries(hop(0, 1));
     }
   }

   // "Non-polarized" means that the coefficients xi do not depend on spin; a global magnetic field is supported
   // through P.globalB. With substeps, the shell has one channel, which is channel M of the problem.
   void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &, const InvarVec &In_,
                    const Opch<SC> &opch, const Coef<SC> &coef) const override {
     using namespace AbelianShell;
     if (P.substeps || P.polarized) my_assert(!P.rungs); // not implemented
     const auto [N, M] = P.substeps ? step.NM() : std::pair{step.N(), std::size_t{0}};
     for (std::size_t sch = 0; sch < shell->channels(); ++sch) {
       const auto ch = P.substeps ? M : sch; // channel of the coefficients and of <||f||>
       if (!P.polarized) {
         add_hopping<SC>(*this, step, offdiag[sch], ch, 0, coef.xi(N, ch), h, qq, In_, opch);
         add_number<SC>(*this, step, *shell, sch, std::nullopt, coef.zeta(N + 1, ch), h, qq);
       } else {
         add_hopping<SC>(*this, step, offdiag_spin[sch][0], ch, 0, coef.xiUP(N, ch), h, qq, In_, opch);
         add_hopping<SC>(*this, step, offdiag_spin[sch][1], ch, 0, coef.xiDOWN(N, ch), h, qq, In_, opch);
         add_number<SC>(*this, step, *shell, sch, Spin::up, coef.zetaUP(N + 1, ch), h, qq);
         add_number<SC>(*this, step, *shell, sch, Spin::down, coef.zetaDOWN(N + 1, ch), h, qq);
       }
     }
     if (P.rungs) {
       my_assert(shell->channels() == 2);
       for (std::size_t ch = 0; ch < 2; ++ch) add_hopping<SC>(*this, step, offdiag_mix[ch], ch, 0, coef.xiR(N, ch), h, qq, In_, opch);
       add_onsite<SC>(*this, step, runghop, coef.zetaR(N + 1, 0), h, qq);
     }
   }

   // <||f||> of channel 'sch' of the shell, stored as channel ch
   void recalc_f_channel(const DiagInfo<SC> &diag, const std::size_t sch, const std::size_t ch, Opch<SC> &opch) const {
     for (const auto &[Ip, eig] : diag) {
       // NOTE: q,ssz only couples to q+1,ssz+-1 in general, even for several channels.
       const int qp   = Ip.get("Q");
       const int sszp = Ip.get("SSZ");
       AbelianShell::recalc_f_block<SC>(*this, diag, qssz(qp + 1, sszp + 1), Ip, ftables[sch][0], opch[ch][0]); // f^dag_up
       AbelianShell::recalc_f_block<SC>(*this, diag, qssz(qp + 1, sszp - 1), Ip, ftables[sch][1], opch[ch][0]); // f^dag_down
     }
   }

   Opch<SC> recalc_irreduc(const Step &, const DiagInfo<SC> &diag) const override {
     Opch<SC> opch(P);
     for (std::size_t ch = 0; ch < shell->channels(); ++ch) recalc_f_channel(diag, ch, ch, opch);
     return opch;
   }

   OpchChannel<SC> recalc_irreduc_substeps(const Step &, const DiagInfo<SC> &diag, int M) const override {
     Opch<SC> opch(P);
     recalc_f_channel(diag, 0, M, opch);
     return opch[M];
   }

   // NOTE: p is ket (right side), 1 is bra (left side). OP is sandwiched in between. Thus Q[p] + Q[op] = Q[1]. In
   // the (Q,S_z) basis, spin up and spin down are not equivalent: two spin-diagonal spectral densities (and also the
   // spin-flip spectral density) can be defined.
   MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
     MatrixElements<SC> cnew;
     for (const auto &[I1, eig] : diag) {
       const int q1   = I1.get("Q");
       const int ssz1 = I1.get("SSZ");
       AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, qssz(q1 - 1, ssz1 + 1), qssz(1, -1), true, cnew);
       AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, qssz(q1 - 1, ssz1 - 1), qssz(1, +1), true, cnew);
     }
     return cnew;
   }

   MatrixElements<SC> recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
     MatrixElements<SC> cnew;
     for (const auto &[I1, eig] : diag) {
       const int q1   = I1.get("Q");
       const int ssz1 = I1.get("SSZ");
       AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, qssz(q1, ssz1), qssz(0, 0), false, cnew);
       AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, qssz(q1, ssz1 + 2), qssz(0, -2), false, cnew);
       // with modulus 4, Delta SSZ = +2 and -2 coincide
       if (modulus != 4) AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, qssz(q1, ssz1 - 2), qssz(0, +2), false, cnew);
     }
     return cnew;
   }

   // Q1u and Q1d refer to channel 0: n_{0,sigma} - 1/2.
   void recalc_global(const Step &, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const override {
     using namespace AbelianShell;
     const t_coef one = 1.0;
     const auto add = [&](const Operator &op) { add_site_operator<SC>(*this, *shell, diag, op, one, cnew); };
     if (name == "SZtot") return add(sum_channels(*shell, spinz));
     if (name == "Q1u") return add(number(0, Spin::up) + (-0.5) * unity());
     if (name == "Q1d") return add(number(0, Spin::down) + (-0.5) * unity());
     my_assert_not_reached();
   }


   void show_coefficients(const Step &step, const Coef<SC> &coef) const override {
     Symmetry<SC>::show_coefficients(step, coef);
     if (P.silent) return;
     if (P.rungs)
       for (unsigned int i = 0; i < P.channels; i++) {
         std::cout << "[" << i + 1 << "]"
           << " xi_rung(" << step.N() << ")=" << coef.xiR(step.N(), i) << " zeta_rung(" << step.N() + 1 << ")=" << coef.zetaR(step.N() + 1, i) << std::endl;
       }
   }
};

}
