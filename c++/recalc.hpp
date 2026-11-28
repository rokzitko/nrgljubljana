#ifndef _recalc_hpp_
#define _recalc_hpp_

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <iterator>
#include <iostream>
#include <optional>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "debug.hpp"
#include "traits.hpp"
#include "invar.hpp"
#include "eigen.hpp"
#include "subspaces.hpp"
#include "operators.hpp"
#include "symmetry.hpp"
#include "numerics.hpp"

namespace NRG {

// We split the matrices of eigenvectors into blocks according to the partition into "ancestor subspaces". At the price
// of some copying, this increases memory locality of data and thus improves numerical performence of matrix-matrix multiplies
// in the recalculation of matrix elements. Note that the original (matrix) data is discarded after the splitting.
template<scalar S, typename Matrix = Matrix_traits<S>>
inline void split_in_blocks_Eigen(Eigen<S> &e, const SubspaceDimensions &sub, const bool discard = true) {
  const auto combs = sub.combs();
  e.U.resize(combs);
  const auto nr = e.getnrstored();
  my_assert(0 < nr && nr <= e.getdim());
  for (const auto block: range0(combs)) {
    Matrix U = e.vectors.submatrix_const({0, nr}, sub.part(block));
    e.U.set(block, std::move(U));
  }
  assert(e.U.is_unitary());
  if (discard)
    e.vectors.shrink();
}

template<scalar S>
inline void split_in_blocks(DiagInfo<S> &diag, const SubspaceStructure &substruct, const bool discard = true) {
  for(auto &[I, eig]: diag)
    split_in_blocks_Eigen(eig, substruct.at(I), discard);
}

template<scalar S, typename Matrix = Matrix_traits<S>>
inline void shrink(Eigen<S> &e, const SubspaceDimensions &) {
  e.vectors.shrink();
}

template<scalar S>
inline void shrink(DiagInfo<S> &diag, const SubspaceStructure &substruct) {
  for(auto &[I, eig]: diag)
    shrink(eig, substruct.at(I));
}

template<scalar S>
inline void h5save_blocks_Eigen(H5Easy::File &fd, const std::string &name, const Eigen<S> &eig, const SubspaceDimensions &sub)
{
  // For certain symmetry types (SPSU2) the same ancestor subspace may contribute more than one
  // time, thus we add a suffix with a counter (nr).
  const auto range = range0(sub.combs());
  for (const auto & [nr, i]: range | ranges::views::enumerate) {
    if (!sub.exists(i)) continue;
    const auto &block = eig.U.get(i);
    my_assert(size1(block) == eig.getnrstored());
    my_assert(size2(block) == sub.rmax(i));
    h5_dump_matrix(fd, name + "/" + sub.ancestor(i).name() + "|nr=" + std::to_string(nr), block);
  }
}

template<scalar S>
inline void h5save_blocks(H5Easy::File &fd, const std::string &name, const DiagInfo<S> &diag, const SubspaceStructure &substruct) {
  for(const auto &[I, eig]: diag)
    h5save_blocks_Eigen(fd, name + I.name(), eig, substruct.at(I));
}

// Rows of a recalculation table that share the ket block ip are summed before the product with U(ip)^dagger. This
// saves matrix products: in the doublet and triplet tables there are on average 1.7 rows per ip, in the <||f||>
// tables 1.2. The rows are processed in table order within each ip group, independently of the grouping.
template<typename Row> auto group_rows_by_ip(std::vector<Row> rows) {
  std::stable_sort(rows.begin(), rows.end(), [](const auto &a, const auto &b) { return a.ip < b.ip; });
  return rows;
}

// Recalculates the irreducible matrix elements <I1|| f || Ip>. Called from recalc_irreduc() in nrg-recalc-* files.
template<scalar S> template<typename T>
auto Symmetry<S>::recalc_f(const DiagInfo<S> &diag,
                           const Invar &I1, // bra
                           const Invar &Ip, // ket
                           const T &table) const
{
  nrglog('f', "*** recalc_f() ** f: I1=(" << I1 << ") Ip=(" << Ip << ")");
  if (!recalc_f_coupled(I1, Ip, this->Invar_f)) return Matrix(0,0); // exception for QST and SPSU2T
  const auto &diagI1 = diag.at(I1);
  const auto &diagIp = diag.at(Ip);
  const auto dim1 = diagI1.getnrstored(); // # of states in I1, i.e. the dimension of the <||f||> matrix.
  const auto dimp = diagIp.getnrstored();
  nrglog('f', "dim1=" << dim1 << " dimp=" << dimp);
  my_assert(dim1 < 1000000 && dimp < 1000000); // bug trap
  auto f = zero_matrix<S>(dim1, dimp);
  const auto backend = mult_backend(P.mult);
  if (backend == MultBackend::cuda) {
    CudaRecalcAccumulator<S, Matrix> cuda_f(f);
    if (cuda_f.active()) {
      for (const auto &[i1, ip, factor]: table) {
        nrglog('f', "** i1=" << i1 << " ip=" << ip << " factor=" << factor);
        const auto &U1 = diagI1.U(i1);
        const auto &Up = diagIp.U(ip);
        nrglog('f', "* norm1=" << frobenius_norm(U1) << " normp=" << frobenius_norm(Up));
        if (finite_size(U1) && finite_size(Up)) nrglog('f', "* U1(0,0)=" << U1(0,0) << " Up(0,0)=" << Up(0,0));
        product_CUDA_accumulate<S>(cuda_f, factor, U1, Up);
      }
      cuda_f.copy_back();
      if (P.logletter('F')) dump_matrix(f);
      return f;
    }
    for (const auto &[i1, ip, factor]: table)
      product<S>(backend, f, factor, diagI1.U(i1), diagIp.U(ip));
    if (P.logletter('F')) dump_matrix(f);
    return f;
  }
  // <I1||f||Ip> gets contributions from various |QSr> states. These are given by i1, ip in the Recalc_f type tables.
  // f += sum_ip [sum_{i1} factor U(i1)] U(ip)^dagger, see group_rows_by_ip().
  const auto rows = group_rows_by_ip(std::vector<Recalc_f<S>>(std::begin(table), std::end(table)));
  Matrix X;
  for (auto it = rows.begin(); it != rows.end();) {
    const auto ip = it->ip;
    const auto &Up = diagIp.U(ip);
    bool have_X = false;
    for (; it != rows.end() && it->ip == ip; ++it) {
      const auto &[i1, ip_, factor] = *it;
      nrglog('f', "** i1=" << i1 << " ip=" << ip_ << " factor=" << factor);
      const auto &U1 = diagI1.U(i1);
      nrglog('f', "* norm1=" << frobenius_norm(U1) << " normp=" << frobenius_norm(Up));
      if (!finite_size(U1) || !finite_size(Up)) continue;
      nrglog('f', "* U1(0,0)=" << U1(0,0) << " Up(0,0)=" << Up(0,0));
      my_assert(size1(U1) == dim1 && size2(U1) == size2(Up));
      my_assert(my_isfinite(factor));
      if (!have_X) {
        X = factor * U1;
        have_X = true;
      } else {
        X += factor * U1;
      }
    }
    if (have_X) gemm_nh<S>(S(1), X, Up, S(1), f);
  }
  if (P.logletter('F')) dump_matrix(f);
  return f;
}

// Recalculate the (irreducible) matrix elements of various operators. This is the most important routine in this
// program, so it is heavily instrumentalized for debugging purposes. It is called from recalc_doublet(),
// recalc_singlet(), and other routines. The inner-most for() loops can be found here, so this is the right spot that
// one should try to hand optimize.
template<scalar S> template<typename T>
std::optional<Matrix_traits<S>> Symmetry<S>::recalc_general(const DiagInfo<S> &diag,
                                                            const SubspaceStructure &substruct,
                                                            const MatrixElements<S> &cold,
                                                            const Invar &I1,             // target subspace (bra)
                                                            const Invar &Ip,             // target subspace (ket)
                                                            const T &table,
                                                            const Invar &Iop) const      // quantum numbers of the operator
{
  if (P.logletter('r')) std::cout << "*** recalc_general: " << nrgdump3(I1, Ip, Iop) << std::endl;
  if (!triangle_inequality(I1, Ip, Iop)) return {};
  const auto &diagI1 = diag.at(I1);
  const auto &diagIp = diag.at(Ip);
  const auto &substructI1 = substruct.at(I1);
  const auto &substructIp = substruct.at(Ip);
  const auto dim1 = diagI1.getnrstored();
  const auto dimp = diagIp.getnrstored();
  auto cn = zero_matrix<S>(dim1, dimp);
  if (dim1 == 0 || dimp == 0) return cn; // return empty matrix
  // The table rows that contribute: cn += factor U(i1) <IN1||O||INp> U(ip)^dagger.
  struct Row {
    size_t i1;
    size_t ip;
    const Matrix *O;
    t_coef factor;
  };
  std::vector<Row> rows;
  rows.reserve(std::size(table));
  for (const auto &[i1, ip, IN1, INp, factor]: table) {
    my_assert(1 <= i1 && i1 <= nr_combs() && 1 <= ip && ip <= nr_combs());
    if (P.logletter('r')) std::cout << nrgdump7(i1, ip, IN1, ancestor(I1,i1-1), INp, ancestor(Ip,ip-1), factor) << std::endl;
    if (!Invar_allowed(IN1) || !Invar_allowed(INp)) continue;
    my_assert(IN1 == ancestor(I1, i1-1));
    my_assert(INp == ancestor(Ip, ip-1));
    const auto rmax1 = substructI1.rmax(i1-1);
    const auto rmaxp = substructIp.rmax(ip-1);
    if (rmax1 == 0 || rmaxp == 0) continue;
    const auto cold_it = cold.find({IN1, INp});
    if (cold_it == cold.end()) continue;
    my_assert(isfinite(factor));
    rows.push_back({i1, ip, &cold_it->second, factor});
  }
  const auto backend = mult_backend(P.mult);
  if (backend == MultBackend::cuda) {
    CudaRecalcAccumulator<S, Matrix> cuda_cn(cn);
    if (cuda_cn.active()) {
      for (const auto &r : rows) transform_CUDA_accumulate<S>(cuda_cn, r.factor, diagI1.U(r.i1), *r.O, diagIp.U(r.ip));
      cuda_cn.copy_back();
    } else {
      for (const auto &r : rows) transform<S>(backend, cn, r.factor, diagI1.U(r.i1), *r.O, diagIp.U(r.ip));
    }
    if (P.logletter('R')) dump_matrix(cn);
    return cn;
  }
  // cn += sum_ip [sum_{i1} factor U(i1) O] U(ip)^dagger, see group_rows_by_ip().
  Matrix X;
  const auto grouped = group_rows_by_ip(std::move(rows));
  for (auto it = grouped.begin(); it != grouped.end();) {
    const auto ip = it->ip;
    const auto &Up = diagIp.U(ip);
    bool have_X = false;
    for (; it != grouped.end() && it->ip == ip; ++it) {
      const auto &U1 = diagI1.U(it->i1);
      if (!finite_size(U1) || !finite_size(Up)) continue;
      my_assert(size2(U1) == size1(*it->O) && size2(*it->O) == size2(Up));
      if (!have_X) {
        X.resize(size1(U1), size2(*it->O));
        gemm_nn<S>(S(it->factor), U1, *it->O, S(0), X);
        have_X = true;
      } else {
        gemm_nn<S>(S(it->factor), U1, *it->O, S(1), X);
      }
    }
    if (have_X) gemm_nh<S>(S(1), X, Up, S(1), cn);
  }
  if (P.logletter('R')) dump_matrix(cn);
  return cn;
}

// This routine is used for recalculation of global operators in nrg-recalc-*.cc
template<scalar S>
void Symmetry<S>::recalc1_global(const DiagInfo<S> &diag, // XXX: pass Eigen instead
                                 const Invar &I,
                                 Matrix &m, // modified, not produced!
                                 const size_t i1,
                                 const size_t ip,
                                 const t_coef value) const
{
  const auto &diagI = diag.at(I);
  product<S>(mult_backend(P.mult), m, value, diagI.U(i1), diagI.U(ip));
}

} // namespace

#endif
