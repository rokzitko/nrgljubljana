// Common implementation of the symmetry types NONE and P. Both use four <||f||> blocks per channel, the creation
// and annihilation operators of both spins of the last site, numbered fnr = 0 (CR, DO), 1 (CR, UP), 2 (AN, DO),
// 3 (AN, UP), as in the data files produced by nrginit. The coefficient tables are generated at run time for any
// number of channels (abelian_shell.hpp).

#ifndef _abelian_cran_hpp_
#define _abelian_cran_hpp_

#include <array>
#include <complex>
#include <cstddef>
#include <optional>
#include <string>
#include <type_traits>
#include <vector>

#include "abelian-common.hpp"

namespace NRG {

template<typename SC>
class SymmetryCRAN : public Symmetry<SC> {
 protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
   using Matrix = typename traits<SC>::Matrix;
   using t_coef = coef_traits<SC>;

   AbelianShell::Kind kind;
   Invar Iferm;                  // quantum numbers of a fermion operator: I1 = Ip.combine(Iferm)
   bool annihilators;            // P: the hopping and <||f||> recalculation also use the annihilation operators
   std::optional<AbelianShell::Shell> shell;
   std::vector<std::array<std::vector<AbelianShell::Entry<double>>, 2>> cr, an; // [ch][spin], upper triangle
   std::vector<std::vector<AbelianShell::Entry<AbelianShell::cmpl>>> isospinx; // [ch], upper triangle
   std::vector<std::array<std::vector<Recalc_f<SC>>, 4>> ftables;             // [ch][fnr]

   static constexpr std::size_t fnr(const bool dagger, const AbelianShell::Spin s) {
     return (dagger ? 0 : 2) + (s == AbelianShell::Spin::down ? 0 : 1);
   }

 public:
   SymmetryCRAN(const Params &P_, const AbelianShell::Kind kind_, const Invar &Isinglet_, const Invar &Iferm_, const bool annihilators_)
     : Symmetry<SC>(P_, std::vector<std::string>{}, Isinglet_), kind(kind_), Iferm(Iferm_), annihilators(annihilators_) {}

   void load() override {
     using namespace AbelianShell;
     shell.emplace(kind, shell_channels(P));
     load_shell(*shell, In, QN, P);
     const auto K = shell->channels();
     cr.assign(K, {});
     an.assign(K, {});
     isospinx.assign(K, {});
     ftables.assign(K, {});
     for (std::size_t ch = 0; ch < K; ++ch) {
       for (const auto s : both_spins) {
         const auto sp = static_cast<std::size_t>(s);
         cr[ch][sp] = shell->hop_cr(c(ch, s), true);
         if (annihilators) an[ch][sp] = shell->hop_an(cdag(ch, s), true);
         ftables[ch][fnr(true, s)] = f_table<SC>(*shell, cdag(ch, s));
         if (annihilators) ftables[ch][fnr(false, s)] = f_table<SC>(*shell, c(ch, s));
       }
       isospinx[ch] = shell->entries(AbelianShell::isospinx(ch), true);
     }
   }

   void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &, const InvarVec &In_,
                    const Opch<SC> &opch, const Coef<SC> &coef) const override {
     using namespace AbelianShell;
     const auto N = step.N();
     for (std::size_t ch = 0; ch < shell->channels(); ++ch) {
       for (const auto s : both_spins) {
         const auto sp = static_cast<std::size_t>(s);
         const t_coef xi = !P.polarized ? coef.xi(N, ch) : (s == Spin::up ? coef.xiUP(N, ch) : coef.xiDOWN(N, ch));
         add_hopping<SC>(*this, step, cr[ch][sp], ch, fnr(true, s), xi, h, qq, In_, opch);   // f^dag_N a_{N+1}
         add_hopping<SC>(*this, step, an[ch][sp], ch, fnr(false, s), xi, h, qq, In_, opch);  // a^dag_{N+1} f_N
       }
       if (!P.polarized) {
         add_number<SC>(*this, step, *shell, ch, std::nullopt, coef.zeta(N + 1, ch), h, qq);
       } else {
         add_number<SC>(*this, step, *shell, ch, Spin::up, coef.zetaUP(N + 1, ch), h, qq);
         add_number<SC>(*this, step, *shell, ch, Spin::down, coef.zetaDOWN(N + 1, ch), h, qq);
       }
       // Pairing Delta f^dag_up f^dag_down + h.c. = 2 Re(Delta) I_x - 2 Im(Delta) I_y; conj(Delta), see Symmetry::uses_pairing().
       add_onsite<SC>(*this, step, isospinx[ch], t_coef(2.0) * conj_me(coef.delta(N + 1, ch)), h, qq);
     }
   }

   Opch<SC> recalc_irreduc(const Step &, const DiagInfo<SC> &diag) const override {
     Opch<SC> opch(P);
     for (const auto &[Ip, eig] : diag) {
       auto I1 = Ip;
       I1.combine(Iferm);
       for (std::size_t ch = 0; ch < shell->channels(); ++ch)
         for (std::size_t n = 0; n < 4; ++n)
           if (n < 2 || annihilators) AbelianShell::recalc_f_block<SC>(*this, diag, I1, Ip, ftables[ch][n], opch[ch][n]);
     }
     return opch;
   }

   MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
     MatrixElements<SC> cnew;
     for (const auto &[I1, eig] : diag) {
       auto Ip = I1;
       auto inv = Iferm;
       inv.inverse();
       Ip.combine(inv);
       AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, Ip, Iferm, true, cnew);
     }
     return cnew;
   }

   void recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const override {
     using namespace AbelianShell;
     const t_coef one = 1.0;
     const t_coef iso = psgn(step.getnn() + 1); // transverse isospin components alternate in sign along the chain
     const auto add = [&](const Operator &op, const t_coef factor) { add_site_operator<SC>(*this, *shell, diag, op, factor, cnew); };
     const auto all = [this](auto op) { return sum_channels(*shell, op); };
     if (name == "SZtot") return add(all(spinz), one);
     if (name == "SXtot") return add(all(spinx), one);
     if (name == "Qtot") return add(all(charge), one);
     if (name == "Iztot") return add(all(isospinz), one);
     if (name == "Ixtot") return add(all(AbelianShell::isospinx), iso);
     if (name == "Iptot") return add(all(isospinplus), iso);
     if (name == "Imtot") return add(all(isospinminus), iso);
     if constexpr (std::is_same_v<SC, std::complex<double>>) {
       if (name == "SYtot") return add(all(spiny), one);
       if (name == "Iytot") return add(all(isospiny), iso);
     }
     my_assert_not_reached();
   }

   bool uses_pairing() const override { return true; }
};

} // namespace NRG

#endif
