// matrix.cc - Symmetry dependent code for Hamiltonian matrix generation
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _matrix_hpp_
#define _matrix_hpp_

#include <stdexcept>

#include "invar.hpp"
#include "eigen.hpp"
#include "symmetry.hpp"
#include "step.hpp"
#include "operators.hpp"
#include "subspaces.hpp"

#define FMT_HEADER_ONLY
#include <fmt/format.h>

namespace NRG {

// +++ Construct an offdiagonal part of the Hamiltonian. +++

// We test if the block (i,j) exists at all. If not, factor is not evaluated. This prevents divisions by zero.
#define offdiag_function(step, i, j, ch, fnr, factor, h, qq, In, opch) \
  { if (qq.offdiag_contributes(i, j)) this->offdiag_function_impl(step, i, j, ch, fnr, factor, h, qq, In, opch); }

// i,j - indexes of the out-of-diagonal matrix block that we are constructing
// ch,fnr - channel and extra index to locate the correct block of the <||f||> matrix (irreducible matrix elements)
// h - matrix being built
// qq - matrix dimensions data (to determine the position in the matrix)
// In - In[i] and In[j] are the invariant subspaces of required <||f||> matrix elements
// factor  - the coefficient which multiplies the irreducible matrix elements. This coefficient takes into account
//           the multiplicities.
// NOTE: the offdiagonal part depends on xi(N), while zeta(N) affect the diagonal part of the Hamiltonian matrix!
template<scalar S>
void Symmetry<S>::offdiag_function_impl(const Step &step, const size_t i, const size_t j,
                                             const size_t ch,      // channel number
                                             const size_t fnr,     // extra index for <||f||>, usually 0
                                             const t_coef factor,  // may be complex (in principle)
                                             Matrix &h, const SubspaceDimensions &qq, const InvarVec &In, const Opch<S> &opch) const
{
  my_assert(1 <= i && i <= qq.combs() && 1 <= j && j <= qq.combs());
  if (!my_isfinite(factor))
    throw std::runtime_error(fmt::format("offdiag_function(): factor not finite {} {} {} {}", i, j, ch, fnr));
  if (const auto f = opch[ch][fnr].find({In[i-1], In[j-1]}); f != opch[ch][fnr].cend()) {   // < In[i] r | f^\dag | In[j] r' >
    const Matrix & mat = f->second;
    my_assert(qq.rmax(i-1) == size1(mat) && qq.rmax(j-1) == size2(mat));
    const auto factor_scaled = factor / step.scale();
    // We are building the upper triangular part of the Hermitian Hamiltonian. Thus usually i < j. If not, we must
    // conjugate transpose the contribution!
    const bool conj_transpose = i > j;
    if (conj_transpose) {
      auto hsub = submatrix(h, qq.part_mma(j), qq.part_mma(i));
      hsub += conj_me(factor_scaled) * herm(mat);
    } else {
      auto hsub = submatrix(h, qq.part_mma(i), qq.part_mma(j));
      hsub += factor_scaled * mat;
    }
  } else
    throw std::runtime_error(fmt::format("offdiag_function(): matrix not found {} {} {} {}", i, j, ch, fnr));
}

// +++ Shift the diagonal matrix elements by the number of electrons multiplied by the required constant(s) zeta. +++
//
// 'number' is the number of electrons for channel 'ch' in invariant subspaces indexed by 'i'. Note that 'number' is
// a floating point number: this is required, for example, in the LR symmetric basis sets.
//
// NOTE: for problems where a given invariant subspace does not correspond to a fixed number of added electrons, a
// generalized routine should be used. 
template<scalar S>
void Symmetry<S>::diag_function_impl(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                                          Matrix &h, const SubspaceDimensions &qq, const double f) const
{
  my_assert(1 <= i && i <= qq.combs());
  // For convenience we subtract the average site occupancy. XXX: how does this affect the total energy??
  const auto avgoccup = (double)P.spin / 2; // multiplicity divided by 2
  // Energy shift of the diagonal matrix elements in the NRG Hamiltonian.
  for (const auto j: qq.view_mma(i)) h(j, j) += sc_zeta * (number - f*avgoccup) / step.scale();
}

template<scalar S>
void Symmetry<S>::diag_function(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                                     Matrix &h, const SubspaceDimensions &qq) const
{
  diag_function_impl(step, i, number, sc_zeta, h, qq, 1);
}

template<scalar S>
void Symmetry<S>::diag_function_half(const Step &step, const size_t i, const double number, const t_coef sc_zeta,
                                          Matrix &h, const SubspaceDimensions &qq) const
{
  diag_function_impl(step, i, number, sc_zeta, h, qq, 0.5);
}

// +++ Shift the offdiagonal matrix elements by factor. +++

template<scalar S>
void Symmetry<S>::diag_offdiag_function(const Step &step, const size_t i, const size_t j, const t_coef factor,
                                             Matrix &h, const SubspaceDimensions &qq) const
{
  my_assert(1 <= i && i <= qq.combs() && 1 <= j && j <= qq.combs());
  if (i > j) return; // only upper triangular part
  const auto begin1 = qq.offset(i-1);
  const auto size1  = qq.rmax(i-1);
  const auto begin2 = qq.offset(j-1);
  const auto size2  = qq.rmax(j-1);
  const auto contributes = (size1 > 0) && (size2 > 0);
  if (!contributes) return;
  my_assert(size1 == size2);
  const t_coef factor_scaled = factor / step.scale();
  for (const auto l: range0(size1)) h(begin1 + l, begin2 + l) += factor_scaled;
}

} // namespace NRG

#endif // _matrix_cc_
