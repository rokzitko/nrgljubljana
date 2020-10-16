#ifndef _recalc_cc_
#define _recalc_cc_

// We split the matrices of eigenvectors in blocks according to the partition into "ancestor subspaces". At the price
// of some copying, this increases memory localisation of data and thus improves numerical performence of gemm calls
// in the recalculation of matrix elements. Note that the original (matrix) data is discarded after the splitting had
// completed!
void split_in_blocks_Eigen(const Invar &I, Eigen &e, const QSrmax &qsrmax) {
  const auto combs = qsrmax.at(I).combs();
  e.blocks.resize(combs);
  const auto nr = e.getnrstored(); // nr. of eigenpairs
  my_assert(nr > 0);
  my_assert(nr <= e.getdim()); // rmax = length of eigenvectors
  for (const auto block: range0(combs)) {
    const auto rmax   = qsrmax.at(I).rmax(block + 1); // offset 1
    const auto offset = qsrmax.at(I).offset(block + 1); // RRR
    my_assert(e.matrix.size1() >= nr);
    my_assert(e.matrix.size2() >= offset + rmax);
    ublas::matrix_range<Matrix> Up(e.matrix, ublas::range(0, nr), ublas::range(offset, offset + rmax));
    e.blocks[block] = Matrix(Up);
    my_assert(e.blocks[block].size1() == nr);
    my_assert(e.blocks[block].size2() == rmax);
  }
  e.matrix = Matrix(0, e.getdim()); // We don't need the matrix anymore, but we keep the information about the dimensionality!!
}

void split_in_blocks(DiagInfo &diag, const QSrmax &qsrmax) {
  for(auto &[I, eig]: diag)
    split_in_blocks_Eigen(I, eig, qsrmax);
}

// Recalculates irreducible matrix elements <I1|| f || Ip>. Called from recalc_irreduc() in nrg-recalc-* files.
Matrix Symmetry::recalc_f(const DiagInfo &diag, 
                          const QSrmax &qsrmax, 
                          const Invar &I1,
                          const Invar &Ip, 
                          const struct Recalc_f table[], 
                          const size_t jmax)
{
  nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
  if (!recalc_f_coupled(I1, Ip, Invar_f)) {
    nrglog('f', "Does not fulfill the triangle inequalities.");
    return Matrix(0,0);
  }
  const Eigen &diagI1 = diag.at(I1);
  const Eigen &diagIp = diag.at(Ip);
  // Number of states in Ip and in I1, i.e. the dimension of the <||f||> matrix of irreducible matrix elements.
  const auto dim1 = diagI1.getnrstored();
  const auto dimp = diagIp.getnrstored();
  nrglog('f', "dim1=" << dim1 << " dimp=" << dimp);
  const Twoinvar II = {I1, Ip};
  Matrix f = Matrix(dim1, dimp, 0);
  if (dim1 && dimp) {
    // <I1||f||Ip> gets contributions from various |QSr> states. These are given by i1, ip in the Recalc_f type tables.
    for (const auto j: range0(jmax)) {
      // rmax1, rmaxp are the dimensions of the invariant subspaces
      const auto rmax1 = qsrmax.at(I1).rmax(table[j].i1);
      const auto rmaxp = qsrmax.at(Ip).rmax(table[j].ip);
      if (!(rmax1 > 0 && rmaxp > 0)) continue;
      if (P.logletter('f'))
        nrgdump6(j, table[j].i1, table[j].ip, table[j].factor, rmax1, rmaxp);
      my_assert(my_isfinite(table[j].factor) && rmax1 == rmaxp);
      const Matrix &U1 = diagI1.blocks[table[j].i1 - 1]; // offset 1.. argh!
      const Matrix &Up = diagIp.blocks[table[j].ip - 1];
      my_assert(U1.size1() == dim1 && Up.size1() == dimp && U1.size2() == Up.size2());
      my_assert(rmax1 == U1.size2() && rmaxp == Up.size2());
      atlas::gemm(CblasNoTrans, CblasConjTrans, t_factor(table[j].factor), U1, Up, t_factor(1.0), f);
    } // loop over j
  }
  if (P.logletter('F')) dump_matrix(f);
  return f;
}

// Recalculate the (irreducible) matrix elements of various operators. This is the most important routine in this
// program, so it is heavily instrumentalized for debugging purposes. It is called from recalc_doublet(),
// recalc_singlet(), and other routines. The inner-most for() loops can be found here, so this is the right spot that
// one should try to hand optimize.

Matrix Symmetry::recalc_general(const DiagInfo &diag, 
                                const QSrmax &qsrmax,        // information about the matrix structure
                                const MatrixElements &cold,
                                const Invar &I1,             // target subspace (bra)
                                const Invar &Ip,             // target subspace (ket)
                                const struct Recalc table[],
                                const size_t jmax,           // length of table
                                const Invar &Iop) const      // quantum numbers of the operator
{
  if (P.logletter('r')) {
    cout << "recalc_general: ";
    nrgdump3(I1, Ip, Iop) << endl;
  }
  const Eigen &diagI1 = diag.at(I1);
  const Eigen &diagIp = diag.at(Ip);
  const auto dim1 = diagI1.getnrstored();
  const auto dimp = diagIp.getnrstored();
  const Twoinvar II = {I1, Ip};
  Matrix cn = Matrix(dim1, dimp, 0);
  if (dim1 == 0 || dimp == 0) return cn; // empty matrix
  for (const auto j: range0(jmax)) { // loop over combinations of i/ip
    if (P.logletter('r')) {
      nrgdump3(j, I1, Ip);
      nrgdump2(table[j].i1, table[j].ip);
      nrgdump2(table[j].IN1, table[j].INp);
      nrgdump(table[j].factor) << endl;
    }
    if (!Invar_allowed(table[j].IN1) || !Invar_allowed(table[j].INp)) continue;
    const auto rmax1 = qsrmax.at(I1).rmax(table[j].i1);
    const auto rmaxp = qsrmax.at(Ip).rmax(table[j].ip);
    // Proceed if this combination of i1/ip contributes.
    if (rmax1 == 0 || rmaxp == 0) continue;
    const auto IN1 = ancestor(I1, table[j].i1-1); // RRR
    const auto INp = ancestor(Ip, table[j].ip-1); // RRR
    my_assert(IN1 == table[j].IN1 && INp == table[j].INp);
    const Twoinvar ININ = {table[j].IN1, table[j].INp};
    const auto cnt    = cold.count(ININ); // Number of (IN1,INp) subspaces.
    my_assert(cnt == 0 || cnt == 1);        // Anything other than 0 or 1 is a bug!
    // There are exceptions when a subspace might not contribute. Two are known: 1. for triplet operators, two
    // singlet states give no contribution [SU(2) symmetry]; 2. some subspaces might not exist at low iteration
    // steps. If triangle inequality is not satisfied and there are indeed no states for the given subspace pair,
    // this is OK and we just skip this case.
    const bool triangle = triangle_inequality(table[j].IN1, Iop, table[j].INp);
    if (!triangle && cnt == 0) continue;
    // All exceptions should be handled by now. If cnt != 1, this is a bug, probably related to symmetry properties
    // and the coefficient tables for NRG transformations of matrices.
    my_assert(cnt == 1);
    my_assert(my_isfinite(table[j].factor));
    if (P.logletter('r')) cout << "Contributes: rmax1=" << rmax1 << " rmaxp=" << rmaxp << endl;
    // RECALL: rmax1 - dimension of the subspace of invariant subspace I1 spanned by the states originating from the
    //  combination |I1>_i1, where i1=1...P.combs. It is clearly equal to the dimension of the invariant subspace IN1
    //  from the previous (N-th) iteration.
    const Matrix &m = cold.at(ININ);     // m: irreducible elements at previous stage
    my_assert(rmax1 == m.size1() && rmaxp == m.size2());
    const Matrix &U1 = diagI1.blocks[table[j].i1 - 1]; // offset 1.. argh!
    const Matrix &Up = diagIp.blocks[table[j].ip - 1];
    my_assert(U1.size1() == dim1 && U1.size2() == rmax1 && Up.size1() == dimp && Up.size2() == rmaxp);
    // &&&& Performace hot-spot. Ensure that you're using highly optimised BLAS library.
    Matrix temp(rmax1, dimp);
    atlas::gemm(CblasNoTrans, CblasConjTrans, t_factor(1.0), m, Up, t_factor(0.0), temp);
    atlas::gemm(CblasNoTrans, CblasNoTrans, t_factor(table[j].factor), U1, temp, t_factor(1.0), cn);
  } // over table
  if (P.logletter('R'))
    std::cout << std::endl << "Matrix dump, I1=" << I1 << " Ip=" << Ip << ":" << std::endl << cn << std::endl;
  return cn;
}

// This routine is used for recalculation of global operators in
// nrg-recalc-*.cc
void Symmetry::recalc1_global(const DiagInfo &diag,
                              const QSrmax &qsrmax,
                              const Invar &I, 
                              Matrix &m, // XXX: return this one
                              const size_t i1, 
                              const size_t ip, 
                              const t_factor value) const {
  const Eigen &diagI = diag.at(I);
  const auto dim = diagI.getnrstored();
  if (dim == 0) return;
  const auto rmax1 = qsrmax.at(I).rmax(i1);
  const auto rmaxp = qsrmax.at(I).rmax(ip);
  my_assert(rmax1 == rmaxp);
  if (rmax1 == 0 || rmaxp == 0) return;
  const Matrix &U1 = diagI.blocks[i1-1];
  const Matrix &Up = diagI.blocks[ip-1];
  my_assert(U1.size1() == dim && U1.size2() == rmax1);
  my_assert(Up.size1() == dim && Up.size2() == rmaxp);
  // m = m + value * U1 * Up^trans
  atlas::gemm(CblasNoTrans, CblasConjTrans, t_factor(value), U1, Up, t_factor(1.0), m);
}

#endif // _recalc_cc_
