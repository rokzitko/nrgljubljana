#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <string>
#include <vector>

#include "abelian-common.hpp"

namespace NRG {

using std::pow;

// U(1) spin symmetry (conserved S_z) for problems with superconducting pairing. Any number of channels; the
// coefficient tables are generated at run time (abelian_shell.hpp). One <||f||> block per channel holds both f^dag_up
// (Delta SSZ = +1) and f^dag_down (Delta SSZ = -1).
//
// The class also implements SPZ4 (sym-SPZ4-impl.hpp), where only 2S_z modulo 4 is conserved: with modulus = 4 the spin
// label SSZ is taken modulo 4, f^dag_up and f^dag_down change it by +1 and +3, and the spin-flip components of a
// triplet (Delta SSZ = +2 and -2) coincide.
template<typename SC>
class SymmetrySPU1 : public SymField<SC> {
 protected:
   using Symmetry<SC>::P;
   using Symmetry<SC>::In;
   using Symmetry<SC>::QN;
   using t_coef = coef_traits<SC>;

   AbelianShell::Kind kind = AbelianShell::Kind::SPU1;
   int modulus = 0; // 0: SSZ is additive (SPU1); 4: SSZ is taken modulo 4 (SPZ4)

   std::optional<AbelianShell::Shell> shell;
   // Hopping f^dag_N a_{N+1}: both spins, and each spin separately (polarized chains); anomalous hopping
   // f^dag_{N,up} a^dag_{N+1,down} - f^dag_{N,down} a^dag_{N+1,up}. All (i,j), see Symmetry::offdiag_function_impl().
   std::vector<std::vector<AbelianShell::Entry<double>>> offdiag, anomalous;
   std::vector<std::array<std::vector<AbelianShell::Entry<double>>, 2>> offdiag_spin;
   std::vector<std::vector<AbelianShell::Entry<AbelianShell::cmpl>>> isospinx; // upper triangle
   std::vector<std::array<std::vector<Recalc_f<SC>>, 2>> ftables;             // [ch][spin]: f^dag_sigma

   // Spin label shifted by d
   [[nodiscard]] Invar ssz(const int x) const { return Invar(modulus ? ((x % modulus) + modulus) % modulus : x); }
   [[nodiscard]] int ssz_difference(const Invar &I1, const Invar &I2) const {
     const int d = I1.get("SSZ") - I2.get("SSZ");
     return modulus ? ((d % modulus) + modulus) % modulus : d;
   }

 public:
   using Matrix = typename traits<SC>::Matrix;
   using t_matel = typename traits<SC>::t_matel;
   template<typename T>
   SymmetrySPU1(const Params &P, const T td_fields) : SymField<SC>(P, td_fields, Invar(0)) {}
   SymmetrySPU1(const Params &P) : SymmetrySPU1(P, std::vector{"<Sz^2>", "<Sz>"}) {
     initInvar({
        {"SSZ", additive} // spin projection
     });
   }

  // The spin projection of the operator is defined by the difference in Sz of both the invariant subspaces.
  bool check_SPIN(const Invar &I1, const Invar &Ip, const int &SPIN) const override {
    const int sszop = ssz_difference(I1, Ip);
    return modulus ? (sszop == 1 && SPIN == 1) || (sszop == 3 && SPIN == -1) : sszop == SPIN;
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
    anomalous.assign(K, {});
    offdiag_spin.assign(K, {});
    isospinx.assign(K, {});
    ftables.assign(K, {});
    for (std::size_t ch = 0; ch < K; ++ch) {
      offdiag[ch] = shell->hop_cr(c(ch, Spin::up) + c(ch, Spin::down));
      anomalous[ch] = shell->hop_cr(cdag(ch, Spin::down) + (-1.0) * cdag(ch, Spin::up));
      for (const auto s : both_spins) {
        offdiag_spin[ch][static_cast<std::size_t>(s)] = shell->hop_cr(c(ch, s));
        ftables[ch][static_cast<std::size_t>(s)] = f_table<SC>(*shell, cdag(ch, s));
      }
      isospinx[ch] = shell->entries(AbelianShell::isospinx(ch), true);
    }
  }

  void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
    bucket trSZ, trSZ2; // Tr[S_z], Tr[S_z^2]
    for (const auto &[I, eig]: diag) {
      const int ssz  = I.get("SSZ");
      const double sumZ = this->calculate_Z(I, eig, factor);
      trSZ += sumZ * SZ(ssz);
      trSZ2 += sumZ * pow(SZ(ssz),2);
    }
    stats.td.set("<Sz^2>", trSZ2 / stats.Z);
    stats.td.set("<Sz>",   trSZ / stats.Z);
  }

  // With substeps, the shell has one channel, which is channel M of the problem.
  void make_matrix(Matrix &h, const Step &step, const SubspaceDimensions &qq, const Invar &, const InvarVec &In_, const Opch<SC> &opch, const Coef<SC> &coef) const override {
    using namespace AbelianShell;
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
      add_hopping<SC>(*this, step, anomalous[sch], ch, 0, coef.kappa(N, ch), h, qq, In_, opch);
      // conj(Delta), see Symmetry::uses_pairing()
      add_onsite<SC>(*this, step, isospinx[sch], t_coef(2.0) * conj_me(coef.delta(N + 1, ch)), h, qq);
    }
  }

  // <||f||> of channel 'sch' of the shell, stored as channel ch
  void recalc_f_channel(const DiagInfo<SC> &diag, const std::size_t sch, const std::size_t ch, Opch<SC> &opch) const {
    for (const auto &[Ip, eig] : diag) {
      const int sszp = Ip.get("SSZ");
      AbelianShell::recalc_f_block<SC>(*this, diag, ssz(sszp + 1), Ip, ftables[sch][0], opch[ch][0]); // f^dag_up
      AbelianShell::recalc_f_block<SC>(*this, diag, ssz(sszp - 1), Ip, ftables[sch][1], opch[ch][0]); // f^dag_down
    }
  }

  Opch<SC> recalc_irreduc(const Step &, const DiagInfo<SC> &diag) const override {
    my_assert(!P.substeps);
    Opch<SC> opch(P);
    for (std::size_t ch = 0; ch < shell->channels(); ++ch) recalc_f_channel(diag, ch, ch, opch);
    return opch;
  }

  OpchChannel<SC> recalc_irreduc_substeps(const Step &, const DiagInfo<SC> &diag, int M) const override {
    my_assert(P.substeps);
    Opch<SC> opch(P);
    recalc_f_channel(diag, 0, M, opch);
    return opch[M];
  }

  // Doublets: fermionic operators changing SSZ by +1 or -1.
  MatrixElements<SC> recalc_doublet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
    MatrixElements<SC> cnew;
    for (const auto &[I1, eig] : diag) {
      const int ssz1 = I1.get("SSZ");
      AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, ssz(ssz1 + 1), ssz(-1), true, cnew);
      AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, ssz(ssz1 - 1), ssz(+1), true, cnew);
    }
    return cnew;
  }

  // Triplets: bosonic operators changing SSZ by 0, +2, or -2 (with modulus 4, +2 and -2 coincide).
  MatrixElements<SC> recalc_triplet(const DiagInfo<SC> &diag, const SubspaceStructure &substruct, const MatrixElements<SC> &cold) const override {
    MatrixElements<SC> cnew;
    for (const auto &[I1, eig] : diag) {
      const int ssz1 = I1.get("SSZ");
      AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, ssz(ssz1), ssz(0), false, cnew);
      AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, ssz(ssz1 + 2), ssz(-2), false, cnew);
      if (modulus != 4) AbelianShell::recalc_tensor<SC>(*this, *shell, diag, substruct, cold, I1, ssz(ssz1 - 2), ssz(+2), false, cnew);
    }
    return cnew;
  }

  // Q1, Q2, ... refer to channels 0 and 1; Qdiff = Q2 - Q1.
  void recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const override {
    using namespace AbelianShell;
    const t_coef one = 1.0;
    const t_coef iso = psgn(step.getnn() + 1); // transverse isospin components alternate in sign along the chain
    const auto add = [&](const Operator &op, const t_coef factor) { add_site_operator<SC>(*this, *shell, diag, op, factor, cnew); };
    const auto all = [this](auto op) { return sum_channels(*shell, op); };
    const auto channel = [this, &name](const std::size_t ch) {
      if (ch >= shell->channels()) throw std::invalid_argument(fmt::format("Operator {} requires at least {} channels.", name, ch + 1));
      return ch;
    };
    if (name == "Qtot") return add(all(charge), one);
    if (name == "Iztot") return add(all(isospinz), one);
    if (name == "Ixtot") return add(all(AbelianShell::isospinx), iso);
    if (name == "Iptot") return add(all(isospinplus), iso);
    if (name == "Imtot") return add(all(isospinminus), iso);
    if (name == "Qdiff") return add(number(channel(1)) + (-1.0) * number(channel(0)), one);
    if (name == "Q1") return add(number(channel(0)), one);
    if (name == "Q2") return add(number(channel(1)), one);
    if (name == "Q1UP") return add(number(channel(0), Spin::up), one);
    if (name == "Q2UP") return add(number(channel(1), Spin::up), one);
    if (name == "Q1DO") return add(number(channel(0), Spin::down), one);
    if (name == "Q2DO") return add(number(channel(1), Spin::down), one);
    my_assert_not_reached();
  }

  bool uses_pairing() const override { return true; }
  bool uses_anomalous_hopping() const override { return true; }
  bool supports_substeps() const override { return true; }
};

}
