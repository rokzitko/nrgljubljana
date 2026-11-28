// Helpers shared by the abelian symmetry types (NONE, P, U1, QSZ, SPU1), whose coefficient tables are generated at
// run time by AbelianShell::Shell (c++/abelian_shell.hpp) for any number of channels.

#ifndef _abelian_common_hpp_
#define _abelian_common_hpp_

#include <complex>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

#include "../abelian_shell.hpp"

namespace NRG::AbelianShell {

template<scalar S> [[nodiscard]] coef_traits<S> to_coef(const cmpl v) {
  if constexpr (std::is_same_v<coef_traits<S>, double>) {
    my_assert(v.imag() == 0.0);
    return v.real();
  } else {
    return v;
  }
}

// Number of channels of the shell added in one NRG step: with substeps, one channel at a time.
[[nodiscard]] inline std::size_t shell_channels(const Params &P) { return P.substeps ? 1 : P.channels; }

// Fill In and QN (1-based, P.combs+1 entries; Symmetry::erase_first() follows) from the shell.
inline void load_shell(const Shell &shell, std::vector<Invar> &In, std::vector<Invar> &QN, const Params &P) {
  if (shell.size() != P.combs)
    throw std::invalid_argument(fmt::format("Shell has {} states, expected combs={}.", shell.size(), P.combs));
  for (std::size_t i = 1; i <= shell.size(); ++i) {
    In[i] = shell.in(i);
    QN[i] = shell.qn(i);
  }
}

// Hopping blocks: h(i,j) += value * coefficient * <In[i]||f(ch, fnr)||In[j]>; for i > j the block is added in
// conjugate-transposed form (see Symmetry::offdiag_function_impl()).
template<scalar S>
void add_hopping(const Symmetry<S> &sym, const Step &step, const std::vector<Entry<double>> &entries, const std::size_t ch,
                 const std::size_t fnr, const coef_traits<S> coefficient, Matrix_traits<S> &h, const SubspaceDimensions &qq,
                 const InvarVec &In, const Opch<S> &opch) {
  for (const auto &[i, j, value] : entries)
    if (qq.offdiag_contributes(i, j)) sym.offdiag_function_impl(step, i, j, ch, fnr, coef_traits<S>(value) * coefficient, h, qq, In, opch);
}

// On-site terms coupling blocks i <= j with the same ancestor: h(i,j) += value * coefficient * 1.
template<scalar S>
void add_onsite(const Symmetry<S> &sym, const Step &step, const std::vector<Entry<cmpl>> &entries, const coef_traits<S> coefficient,
                Matrix_traits<S> &h, const SubspaceDimensions &qq) {
  for (const auto &[i, j, value] : entries) sym.diag_offdiag_function(step, i, j, to_coef<S>(value) * coefficient, h, qq);
}

// Diagonal energy shift zeta*(n - 1) for each channel (half=false) or zeta*(n_sigma - 1/2) for each spin (half=true).
template<scalar S>
void add_number(const Symmetry<S> &sym, const Step &step, const Shell &shell, const std::size_t ch, const std::optional<Spin> spin,
                const coef_traits<S> zeta, Matrix_traits<S> &h, const SubspaceDimensions &qq) {
  for (std::size_t i = 1; i <= shell.size(); ++i) {
    if (spin)
      sym.diag_function_half(step, i, shell.electrons(i, ch, *spin), zeta, h, qq);
    else
      sym.diag_function(step, i, shell.electrons(i, ch), zeta, h, qq);
  }
}

// Table <i|op|j> of an operator of the added site for recalc_f().
template<scalar S> [[nodiscard]] std::vector<Recalc_f<S>> f_table(const Shell &shell, const Operator &op) {
  std::vector<Recalc_f<S>> t;
  for (const auto &[i, j, value] : shell.entries(op)) t.push_back({i, j, to_coef<S>(value)});
  return t;
}

// <I1||f||Ip> of the added site, stored in 'target' if both subspaces have stored states.
template<scalar S>
void recalc_f_block(const Symmetry<S> &sym, const DiagInfo<S> &diag, const Invar &I1, const Invar &Ip,
                    const std::vector<Recalc_f<S>> &table, MatrixElements<S> &target) {
  if (sym.have_stored_coupled_subspaces(diag, I1, Ip)) target[{I1, Ip}] = sym.recalc_f(diag, I1, Ip, table);
}

// Recalculation of an operator of the previous sites with quantum numbers Iop between subspaces I1 and Ip
// (Ip combined with Iop gives I1). The added-site state is unchanged; fermionic operators acquire (-1)^n_i.
template<scalar S>
void recalc_tensor(const Symmetry<S> &sym, const Shell &shell, const DiagInfo<S> &diag, const SubspaceStructure &substruct,
                   const MatrixElements<S> &cold, const Invar &I1, const Invar &Ip, const Invar &Iop, const bool fermionic,
                   MatrixElements<S> &cnew) {
  if (!sym.have_stored_subspaces(diag, I1, Ip)) return;
  std::vector<Recalc<S>> table;
  table.reserve(shell.size());
  for (const auto &[i, ip, factor] : shell.recalc_rows(fermionic))
    table.push_back({i, ip, sym.ancestor(I1, i - 1), sym.ancestor(Ip, ip - 1), coef_traits<S>(factor)});
  if (auto cn = sym.recalc_general(diag, substruct, cold, I1, Ip, table, Iop)) cnew[{I1, Ip}] = std::move(*cn);
}

// Contribution of an operator of the added site, which conserves the quantum numbers, to a global operator:
// <I||O||I> += factor * sum_ij <i|O|j> U(i) U(j)^dagger. The previous-site part has been added by recalc_singlet().
template<scalar S>
void add_site_operator(const Symmetry<S> &sym, const Shell &shell, const DiagInfo<S> &diag, const Operator &op,
                       const coef_traits<S> factor, MatrixElements<S> &cnew) {
  const auto entries = shell.entries(op);
  for (const auto &[I, eig] : diag) {
    auto &cn = cnew[{I, I}];
    for (const auto &[i, j, value] : entries) sym.recalc1_global(diag, I, cn, i, j, factor * to_coef<S>(value));
  }
}

// Sum of a channel operator over all channels of the shell.
template<typename F> [[nodiscard]] Operator sum_channels(const Shell &shell, F op) {
  Operator r;
  for (std::size_t ch = 0; ch < shell.channels(); ++ch) r = r + op(ch);
  return r;
}

} // namespace NRG::AbelianShell

#endif
