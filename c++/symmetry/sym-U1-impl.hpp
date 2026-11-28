#include <array>
#include <complex>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "abelian-common.hpp"

namespace NRG {

// U(1) charge symmetry. Any number of channels; the coefficient tables are generated at run time (abelian_shell.hpp).
// Two <||f||> blocks per channel, as written by nrginit: fnr = 0 for f^dag_up, fnr = 1 for f^dag_down.
template<typename SC>
class SymmetryU1 : public Symmetry<SC> {
 private:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
   using t_coef = coef_traits<SC>;

   static constexpr std::size_t fnr_up = 0;
   static constexpr std::size_t fnr_down = 1;

   std::optional<AbelianShell::Shell> shell;
   // Hopping f^dag_{N,sigma} a_{N+1,sigma'}: [ch][sigma'] (all (i,j), see Symmetry::offdiag_function_impl())
   std::vector<std::array<std::vector<AbelianShell::Entry<double>>, 2>> offdiag;
   std::vector<std::vector<AbelianShell::Entry<AbelianShell::cmpl>>> spinflip; // a^dag_down a_up, upper triangle
   std::vector<std::array<std::vector<Recalc_f<SC>>, 2>> ftables;             // [ch][fnr]

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   SymmetryU1(const Params &P) : Symmetry<SC>(P, std::vector{"<Q>", "<Q^2>"}, Invar(0)) {
     initInvar({
       {"Q", additive} // charge
     });
   }

  bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) const override { return u1_equality(I1.get("Q"), I2.get("Q"), I3.get("Q")); }

  void load() override {
    using namespace AbelianShell;
    shell.emplace(Kind::U1, shell_channels(P));
    load_shell(*shell, In, QN, P);
    const auto K = shell->channels();
    offdiag.assign(K, {});
    spinflip.assign(K, {});
    ftables.assign(K, {});
    for (std::size_t ch = 0; ch < K; ++ch) {
      offdiag[ch][0] = shell->hop_cr(c(ch, Spin::up));
      offdiag[ch][1] = shell->hop_cr(c(ch, Spin::down));
      spinflip[ch] = shell->entries(product(cdag(ch, Spin::down), c(ch, Spin::up)), true);
      ftables[ch][fnr_up] = f_table<SC>(*shell, cdag(ch, Spin::up));
      ftables[ch][fnr_down] = f_table<SC>(*shell, cdag(ch, Spin::down));
    }
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
    for (const auto &[I, eig]: diag) {
      const int q    = I.get("Q");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trQ += sumZ * q;
      trQ2 += sumZ * q * q;
    }
    stats.td.set("<Q>",   trQ / stats.Z);
    stats.td.set("<Q^2>", trQ2 / stats.Z);
  }

  // Non-polarized, polarized (spin-dependent coefficients), and full 2x2 spin structure (pol2x2) of the chain. With
  // substeps (non-polarized only), the shell has one channel, which is channel M of the problem.
  void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &, const InvarVec &In_, const Opch<SC> &opch, const Coef<SC> &coef) const override {
    using namespace AbelianShell;
    const bool spin_dependent = P.polarized || P.pol2x2;
    if (P.substeps) my_assert(!spin_dependent && P.coeffactor == 1); // not implemented
    const auto [N, M] = P.substeps ? step.NM() : std::pair{step.N(), std::size_t{0}};
    for (std::size_t sch = 0; sch < shell->channels(); ++sch) {
      const auto ch = P.substeps ? M : sch; // channel of the coefficients and of <||f||>
      const t_coef xi_up = spin_dependent ? coef.xiUP(N, ch) : coef.xi(N, ch);
      const t_coef xi_down = spin_dependent ? coef.xiDOWN(N, ch) : coef.xi(N, ch);
      add_hopping<SC>(*this, step, offdiag[sch][0], ch, fnr_up, xi_up, h, qq, In_, opch);     // f^dag_up a_up
      add_hopping<SC>(*this, step, offdiag[sch][1], ch, fnr_down, xi_down, h, qq, In_, opch); // f^dag_down a_down
      if (P.pol2x2) {
        add_hopping<SC>(*this, step, offdiag[sch][1], ch, fnr_up, coef.xiUPDO(N, ch), h, qq, In_, opch);   // f^dag_up a_down
        add_hopping<SC>(*this, step, offdiag[sch][0], ch, fnr_down, coef.xiDOUP(N, ch), h, qq, In_, opch); // f^dag_down a_up
      }
      if (!spin_dependent) {
        add_number<SC>(*this, step, *shell, sch, std::nullopt, coef.zeta(N + 1, ch), h, qq);
      } else {
        add_number<SC>(*this, step, *shell, sch, Spin::up, coef.zetaUP(N + 1, ch), h, qq);
        add_number<SC>(*this, step, *shell, sch, Spin::down, coef.zetaDOWN(N + 1, ch), h, qq);
      }
      // a^dag_down a_up on the added site; the Hermitian conjugate is implied (upper triangle)
      if (P.pol2x2) add_onsite<SC>(*this, step, spinflip[sch], coef.zetaDOUP(N + 1, ch), h, qq);
    }
  }

  // <||f||> of channel 'sch' of the shell, stored as channel ch
  void recalc_f_channel(const DiagInfo<SC> &diag, const std::size_t sch, const std::size_t ch, Opch<SC> &opch) const {
    for (const auto &[Ip, eig] : diag) {
      const Invar I1(Ip.get("Q") + 1);
      for (const auto n : {fnr_up, fnr_down}) AbelianShell::recalc_f_block<SC>(*this, diag, I1, Ip, ftables[sch][n], opch[ch][n]);
    }
  }

  Opch<SC> recalc_irreduc(const Step &, const DiagInfo<SC> &diag) const override {
    my_assert(!P.substeps);
    Opch<SC> opch(P);
    for (std::size_t ch = 0; ch < shell->channels(); ++ch) recalc_f_channel(diag, ch, ch, opch);
    return opch;
  }

  // The <||f||> of the other channels are recalculated as doublets (see recalc_irreducible()).
  OpchChannel<SC> recalc_irreduc_substeps(const Step &, const DiagInfo<SC> &diag, int M) const override {
    my_assert(P.substeps);
    Opch<SC> opch(P);
    recalc_f_channel(diag, 0, M, opch);
    return opch[M];
  }

  bool supports_substeps() const override { return true; }

  MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
    MatrixElements<SC> cnew;
    for (const auto &[I1, eig] : diag)
      AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, Invar(I1.get("Q") - 1), Invar(1), true, cnew);
    return cnew;
  }

  void recalc_global(const Step &, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const override {
    using namespace AbelianShell;
    const t_coef one = 1.0;
    const auto add = [&](const Operator &op) { add_site_operator<SC>(*this, *shell, diag, op, one, cnew); };
    if (name == "SZtot") return add(sum_channels(*shell, spinz));
    if (name == "SXtot") return add(sum_channels(*shell, spinx));
    if constexpr (std::is_same_v<SC, std::complex<double>>)
      if (name == "SYtot") return add(sum_channels(*shell, spiny));
    my_assert_not_reached();
  }
};

}
