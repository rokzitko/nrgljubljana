// matrix.cc - Symmetry dependent code for Hamiltonian matrix generation
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _matrix_cc_
#define _matrix_cc_

// +++ Construct an offdiagonal part of the Hamiltonian. +++

inline bool Symmetry::offdiag_contributes(const size_t i, const size_t j, const Rmaxvals &qq) const { // i,j are 1-based
  my_assert(1 <= i && i <= qq.combs() && 1 <= j && j <= qq.combs());
  my_assert(i != j);
  return qq.exists(i-1) && qq.exists(j-1);
}

// We test if the block (i,j) exists at all. If not, factor is not evaluated. This prevents divisions by zero.
#define offdiag_function(step, i, j, ch, fnr, factor, h, qq, In, opch) \
  { if (offdiag_contributes(i, j, qq)) offdiag_function_impl(step, i, j, ch, fnr, factor, h, qq, In, opch); }

// i,j - indexes of the out-of-diagonal matrix block that we are constructing
// ch,fnr - channel and extra index to locate the correct block of the <||f||> matrix (irreducible matrix elements)
// h - matrix being built
// qq - matrix dimensions data (to determine the position in the matrix)
// In - In[i] and In[j] are the invariant subspaces of required <||f||> matrix elements
// factor  - the coefficient which multiplies the irreducible matrix elements. This coefficient takes into account
//           the multiplicities.
// NOTE: the offdiagonal part depends on xi(N), while zeta(N) affect the diagonal part of the Hamiltonian matrix! 
void Symmetry::offdiag_function_impl(const Step &step, const size_t i, const size_t j,
                                     const size_t ch,      // channel number
                                     const size_t fnr,     // extra index for <||f||>, usually 0
                                     const t_matel factor, // may be complex (in principle)
                                     Matrix &h, const Rmaxvals &qq, const InvarVec &In, const Opch &opch) const 
{
  my_assert(1 <= i && i <= qq.combs() && 1 <= j && j <= qq.combs());
  if (!my_isfinite(factor))
    throw std::runtime_error(fmt::format("offdiag_function(): factor not finite {} {} {} {}", i, j, ch, fnr));
  const auto begin1 = qq.offset(i-1);
  const auto begin2 = qq.offset(j-1);
  const auto size1  = qq.rmax(i-1);
  const auto size2  = qq.rmax(j-1);
  // already checked in offdiag_contributes(), but we do it again out of paranoia...
  my_assert(size1 && size2);
  // < In[i] r | f^\dag | In[j] r' >
  if (const auto f = opch[ch][fnr].find({In[i-1], In[j-1]}); f != opch[ch][fnr].cend()) {
    const auto &mat = f->second;
    my_assert(size1 == mat.size1() && size2 == mat.size2());
    const auto factor_scaled = factor / step.scale();
    // We are building the upper triangular part of the Hermitian Hamiltonian. Thus usually i < j. If not, we must
    // conjugate transpose the contribution!
    const bool conj_transpose = i > j;
    if (conj_transpose) {
      ublas::matrix_range<Matrix> hsub(h, ublas::range(begin2, begin2 + size2), ublas::range(begin1, begin1 + size1));
      noalias(hsub) += CONJ_ME(factor_scaled) * herm(mat);
    } else {
      ublas::matrix_range<Matrix> hsub(h, ublas::range(begin1, begin1 + size1), ublas::range(begin2, begin2 + size2));
      noalias(hsub) += factor_scaled * mat;
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
void Symmetry::diag_function_impl(const Step &step, const size_t i, const size_t ch, const double number, const t_coef sc_zeta,
                                  Matrix &h, const Rmaxvals &qq, const double f) const 
{
  my_assert(1 <= i && i <= qq.combs());
  const auto begin1 = qq.offset(i-1);
  const auto size1  = qq.rmax(i-1);
  // For convenience we subtract the average site occupancy.
  const auto avgoccup = (double)P.spin / 2; // multiplicity divided by 2
  // Energy shift of the diagonal matrix elements in the NRG Hamiltonian. WARNING: for N=0, we are not adding the
  //  first site of the Wilson chain (indexed as 0), but the second one (indexed as 1). Therefore the appropriate
  //  zeta is not zeta(0), but zeta(1). zeta(0) is the shift applied to the f[0] orbital in initial.m !!!
  const auto shift = sc_zeta * (number - f*avgoccup) / step.scale();
  for (const auto j: boost::irange(begin1, begin1+size1)) h(j, j) += shift;
}

void Symmetry::diag_function(const Step &step, const size_t i, const size_t ch, const double number, const t_matel sc_zeta, 
                             Matrix &h, const Rmaxvals &qq) const 
{
  diag_function_impl(step, i, ch, number, sc_zeta, h, qq, 1);
}

void Symmetry::diag_function_half(const Step &step, const size_t i, const size_t ch, const double number, const t_matel sc_zeta, 
                                  Matrix &h, const Rmaxvals &qq) const 
{
  diag_function_impl(step, i, ch, number, sc_zeta, h, qq, 0.5);
}

// +++ Shift the offdiagonal matrix elements by factor. +++

void Symmetry::diag_offdiag_function(const Step &step, const size_t i, const size_t j, const size_t chin, const t_matel factor, 
                                     Matrix &h, const Rmaxvals &qq) const 
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
  const auto factor_scaled = factor / step.scale();
  for (const auto l: range0(size1)) h(begin1 + l, begin2 + l) += factor_scaled;
}

#endif // _matrix_cc_
