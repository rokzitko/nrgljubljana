// matrix.cc - Symmetry dependent code for Hamiltonian matrix generation
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _matrix_cc_
#define _matrix_cc_

// TO DO: remove P., move all the required information to the
// symmetry class where it belongs.

/* +++ Construct an offdiagonal part of the Hamiltonian. +++

 i,j - indexes of matrix blocks that constitue the total Hamiltonian matrix
       in the given invariant subspace. Not all possible blocks do contribute.
 ch - contribution of which channel is being considered? 0,1,..
 factor0 - the coefficient which multiplies the irreducible matrix
           elements. This coefficient takes into account the multiplicities.
           For simple symmetries, this factor is +-1, where the sign is
           related to the fermionic anticommutations that are required to
           shift around the creation/annihilation operators that appear in the
           hopping operator.

 NOTE: the offdiagonal part depends on xi(N), while zeta(N) affect the
 diagonal part of the Hamiltonian matrix! 
*/

bool offdiag_contributes(size_t i, size_t j, size_t ch, const Rmaxvals &qq) {
  P.allowed_ijch(i, j, ch);
  my_assert(i != j);
  const size_t size1 = qq.rmax(i);
  const size_t size2 = qq.rmax(j);
  // Here we test if the block (i,j) exists at all. If it doesn't, we
  // immediately return, otherwise the following assertion will trigger a
  // false error message and abort the execution.
  const bool contributes = (size1 > 0) && (size2 > 0);
  return contributes;
}

// i,j - indexes of the out-of-diagonal matrix block that we are
//       constructing
// ch,fnr - channel and extra index to locate the correct block of
//          the <||f||> matrix (irreducible matrix elements)
// factor - appropriate factor (includes the Clebsch-Gordan coefficients
//          as well as the xi hopping parameter)
// h - matrix being built
// qq - matrix dimensions data (to determine the position in the matrix)
// In - In[i] and In[j] are the invariant subspaces of required <||f||>
//      matrix elements
void offdiag_function(size_t i, size_t j,
                      size_t ch,      // channel number
                      size_t fnr,     // extra index for <||f||>, usually 0
                      t_matel factor, // may be complex (in principle)
                      Matrix &h, const Rmaxvals &qq, const InvarVec &In, const Opch &opch) {
  if (!offdiag_contributes(i, j, ch, qq))
    return;
  if (!my_isfinite(factor)) {
    cout << "offdiag_function() critical error: factor is not finite." << endl;
    cout << "i=" << i << " j=" << j << " ch=" << ch << " fnr=" << fnr << " factor=" << factor << endl;
    exit(1);
  }
  const size_t begin1 = qq.offset(i);
  const size_t begin2 = qq.offset(j);
  const size_t size1  = qq.rmax(i);
  const size_t size2  = qq.rmax(j);
  // already checked in offdiag_contributes(), but we do it again out of paranoia...
  my_assert(size1 && size2);
  // < In[i] r | f^\dag | In[j] r' >
  const Twoinvar II {In[i], In[j]};
  const t_matel factor_scaled = factor / step.scale();
  const auto f = opch[ch][fnr].find(II);
  if (f == opch[ch][fnr].cend()) {
    cout << "offdiag_function() critical error: subspace does not exist." << endl;
    cout << "i=" << i << " j=" << j << " ch=" << ch << " fnr=" << fnr << " factor_scaled=" << factor_scaled << endl;
    cout << "source subspaces II=" << II << endl;
    exit(1);
  }
  const auto &mat = f->second;
  my_assert(size1 == mat.size1() && size2 == mat.size2());
  nrglog('i', "offdiag i=" << i << " j=" << j << " factor_scaled=" << factor_scaled << " II=" << II);
  // We are building the upper triangular part of the Hermitian Hamiltonian.
  // Thus usually i < j. If not, we must conjugate transpose the contribution!
  const bool conj_transpose = i > j;
  if (conj_transpose) {
    matrix_range<Matrix> hsub(h, range(begin2, begin2 + size2), range(begin1, begin1 + size1));
    noalias(hsub) += CONJ_ME(factor_scaled) * herm(mat);
  } else {
    matrix_range<Matrix> hsub(h, range(begin1, begin1 + size1), range(begin2, begin2 + size2));
    noalias(hsub) += factor_scaled * mat;
  }
}

/* +++ Shift the diagonal matrix elements by the number of electrons
 multiplied by the required constant(s) zeta. +++

 'number' is the number of electrons for channel 'ch' in invariant
 subspaces indexed by 'i'. Note that 'number' is a floating point
 number: this is required, for example, in the LR symmetric basis
 sets.

 NOTE: for problems where a given invariant subspace does not
 correspond to a fixed number of added electrons, a generalized
 routine should be used. 
*/
void diag_function(size_t i, size_t ch, double number, t_coef sc_zeta, Matrix &h, const Rmaxvals &qq) {
  P.allowed_block_index(i);
  my_assert(number >= 0.0 && number <= 14.0);
  const size_t begin1 = qq.offset(i);
  const size_t size1  = qq.rmax(i);
  // For convenience we subtract the average site occupancy.
  const double avgoccup = ((double)P.spin) / 2; // multiplicity divided by 2
  /* Energy shift of the diagonal matrix elements in the NRG Hamiltonian.
   WARNING: for N=0, we are not adding the first site of the Wilson chain
   (indexed as 0), but the second one (indexed as 1). Therefore the
   appropriate zeta is not zeta(0), but zeta(1). zeta(0) is the shift
   applied to the f[0] orbital in initial.m !!! */
  const t_coef shift = sc_zeta * (number - avgoccup) / step.scale();
  nrglog('i', "diag i=" << i << " shift=" << shift);
  for (size_t j = begin1; j < begin1 + size1; j++) h(j, j) += shift;
}

// Compare with diag_function()
void diag_function_half(size_t i, size_t ch, double number, t_matel sc_zeta, Matrix &h, const Rmaxvals &qq) {
  P.allowed_block_index(i);
  my_assert(0.0 <= number && number <= P.spin);
  const size_t begin1 = qq.offset(i);
  const size_t size1  = qq.rmax(i);
  // For convenience we subtract the average site occupancy.
  const double avgoccup = ((double)P.spin) / 2; // multiplicity divided by 2
  // avgoccup is divided by a further factor of 2 compared
  // to diag_function() above!
  const t_matel shift = sc_zeta * (number - avgoccup / 2) / step.scale();
  nrglog('i', "diag_half i=" << i << " shift=" << shift);
  for (size_t j = begin1; j < begin1 + size1; j++) h(j, j) += shift;
}

// Compare with diag_function() above.
void spinz_function(size_t i, size_t j, size_t ch, t_matel spinz, Matrix &h, const Rmaxvals &qq) {
  P.allowed_ijch(i, j, ch);
  my_assert(i == j);

  // compare with the ISOSPINX macro
  const t_matel shift = spinz * double(P.globalB) / step.scale();

  const size_t begin1 = qq.offset(i);
  const size_t size1  = qq.rmax(i);
  nrglog('i', "spinz i=" << i << " shift=" << shift);
  for (size_t k = begin1; k < begin1 + size1; k++) h(k, k) += shift;
}

void spinx_function(size_t i, size_t j, size_t ch, t_matel spinx, Matrix &h, const Rmaxvals &qq) {
  P.allowed_ijch(i, j, ch);
  const t_matel shift = spinx * double(P.globalBx) / step.scale();
  if (i > j) return; // only upper triangular part
  size_t begin1    = qq.offset(i);
  size_t size1     = qq.rmax(i);
  size_t begin2    = qq.offset(j);
  size_t size2     = qq.rmax(j);
  bool contributes = (size1 > 0) && (size2 > 0);
  if (!contributes) return;
  my_assert(size1 == size2);
  nrglog('i', "spinx i=" << i << " shift=" << shift);
  for (size_t l = 0; l < size1; l++) h(begin1 + l, begin2 + l) += shift;
}

// +++ Shift the offdiagonal matrix elements by factor. +++

void diag_offdiag_function(size_t i, size_t j, size_t chin, t_matel factor, Matrix &h, const Rmaxvals &qq) {
  P.allowed_ijch(i, j, chin);
  if (i > j) return; // only upper triangular part
  size_t begin1    = qq.offset(i);
  size_t size1     = qq.rmax(i);
  size_t begin2    = qq.offset(j);
  size_t size2     = qq.rmax(j);
  bool contributes = (size1 > 0) && (size2 > 0);
  if (!contributes) return;
  my_assert(size1 == size2);
  const t_matel factor_scaled = factor / step.scale();
  nrglog('i', "diag_offdiag i=" << i << " j=" << j << " factor_scaled=" << factor_scaled);
  for (size_t l = 0; l < size1; l++) h(begin1 + l, begin2 + l) += factor_scaled;
}

#endif // _matrix_cc_
