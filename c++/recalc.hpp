#ifndef _recalc_hpp_
#define _recalc_hpp_

#include "debug.hpp"
#include "traits.hpp"
#include "invar.hpp"
#include "eigen.hpp"
#include "subspaces.hpp"
#include "operators.hpp"
#include "symmetry.hpp"
#include "numerics.hpp"

namespace NRG {

// We split the matrices of eigenvectors in blocks according to the partition into "ancestor subspaces". At the price
// of some copying, this increases memory localisation of data and thus improves numerical performence of gemm calls
// in the recalculation of matrix elements. Note that the original (matrix) data is discarded after the splitting had
// completed!
template<scalar S, typename Matrix = Matrix_traits<S>>
inline void split_in_blocks_Eigen(Eigen<S> &e, const SubspaceDimensions &sub) {
  const auto combs = sub.combs();
  e.blocks.resize(combs);
  const auto nr = e.getnrstored(); // nr. of eigenpairs
  my_assert(0 < nr && nr <= e.getdim()); // dim = length of eigenvectors
  for (const auto block: range0(combs)) {
    const auto Up = submatrix(e.matrix, {0, nr}, sub.part(block));
    e.blocks[block] = Up;
  }
  e.matrix = Matrix(0, e.getdim()); // We don't need the matrix anymore, but we keep the information about the dimensionality!!
}

template<scalar S>
inline void split_in_blocks(DiagInfo<S> &diag, const SubspaceStructure &substruct) {
  for(auto &[I, eig]: diag)
    split_in_blocks_Eigen(eig, substruct.at(I));
}

template<scalar S>
inline void h5save_blocks_Eigen(H5Easy::File &fd, const std::string &name, const Eigen<S> &eig, const SubspaceDimensions &sub)
{
  for (const auto i: range0(sub.combs()))
    h5_dump_matrix(fd, name + "/" + sub.ancestor(i).name(), eig.blocks[i]);
}

template<scalar S>
inline void h5save_blocks(H5Easy::File &fd, const std::string &name, const DiagInfo<S> &diag, const SubspaceStructure &substruct) {
  for(const auto &[I, eig]: diag)
    h5save_blocks_Eigen(fd, name + I.name(), eig, substruct.at(I));
}

// Recalculates the irreducible matrix elements <I1|| f || Ip>. Called from recalc_irreduc() in nrg-recalc-* files.
template<scalar S> template<typename T>
auto Symmetry<S>::recalc_f(const DiagInfo<S> &diag,
                           const SubspaceStructure &substruct,
                           const Invar &I1, // bra
                           const Invar &Ip, // ket
                           const T &table) const
{
  nrglog('f', "recalc_f() ** f: I1=(" << I1 << ") Ip=(" << Ip << ")");
  if (!recalc_f_coupled(I1, Ip, this->Invar_f)) return Matrix(0,0); // exception for QST and SPSU2T
  const auto & [diagI1, diagIp] = diag.subs(I1, Ip);
  const auto & [dim1, dimp]     = diag.dims(I1, Ip);   // # of states in Ip and in I1, i.e. the dimension of the <||f||> matrix.
  nrglog('f', "dim1=" << dim1 << " dimp=" << dimp);
  Matrix f = Matrix(dim1, dimp, 0);
  // <I1||f||Ip> gets contributions from various |QSr> states. These are given by i1, ip in the Recalc_f type tables.
  for (const auto &[i1, ip, factor]: table)
    product<S>(f, factor, diagI1.Ublock(i1), diagIp.Ublock(ip));
  if (P.logletter('F')) dump_matrix(f);
  return f;
}

// Recalculate the (irreducible) matrix elements of various operators. This is the most important routine in this
// program, so it is heavily instrumentalized for debugging purposes. It is called from recalc_doublet(),
// recalc_singlet(), and other routines. The inner-most for() loops can be found here, so this is the right spot that
// one should try to hand optimize.
template<scalar S> template<typename T>
auto Symmetry<S>::recalc_general(const DiagInfo<S> &diag,
                                 const SubspaceStructure &substruct, // XXX: drop
                                 const MatrixElements<S> &cold,
                                 const Invar &I1,             // target subspace (bra)
                                 const Invar &Ip,             // target subspace (ket)
                                 const T &table,
                                 const Invar &Iop) const      // quantum numbers of the operator
{
  if (P.logletter('r')) std::cout << "recalc_general: " << nrgdump3(I1, Ip, Iop) << std::endl;
  const auto & [diagI1, diagIp] = diag.subs(I1, Ip);
  const auto & [dim1, dimp]     = diag.dims(I1, Ip);
  const Twoinvar II = {I1, Ip};
  auto cn = Matrix(dim1, dimp, 0);
  if (dim1 == 0 || dimp == 0) return cn; // return empty matrix
  for (const auto &[i1, ip, IN1, INp, factor]: table) {
    my_assert(1 <= i1 && i1 <= nr_combs() && 1 <= ip && ip <= nr_combs());
    if (P.logletter('r')) std::cout << nrgdump5(i1, ip, IN1, INp, factor) << std::endl;
    if (!Invar_allowed(IN1) || !Invar_allowed(INp)) continue;
    my_assert(IN1 == ancestor(I1, i1-1) && INp == ancestor(Ip, ip-1));
    const Twoinvar ININ = {IN1, INp};
    if (cold.count(ININ) == 0) continue;
    transform<S>(cn, factor, diagI1.Ublock(i1), cold.at(ININ), diagIp.Ublock(ip));
  } // over table
  if (P.logletter('R')) dump_matrix(cn);
  return cn;
}
   
// This routine is used for recalculation of global operators in nrg-recalc-*.cc
template<scalar S>
void Symmetry<S>::recalc1_global(const DiagInfo<S> &diag, // XXX: pass Eigen instead
                                 const SubspaceStructure &substruct, // XXX: drop
                                 const Invar &I,
                                 Matrix &m, // modified, not produced!
                                 const size_t i1,
                                 const size_t ip,
                                 const t_coef value) const
{
  const auto &diagI = diag.at(I);
  product<S>(m, value, diagI.Ublock(i1), diagI.Ublock(ip));
}

} // namespace

#endif
