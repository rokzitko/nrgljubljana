// dmnrg.h - Density-matrix NRG
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _dmnrg_h_
#define _dmnrg_h_

// Choose one!
#define LINEBYLINE
//#define WHOLEMATRIX

#ifdef WHOLEMATRIX
void saveMatrix(boost::archive::binary_oarchive &oa, const Matrix &m) { oa << m; }
void loadMatrix(boost::archive::binary_iarchive &ia, Matrix &m) { ia >> m; }
void saveEigen(boost::archive::binary_oarchive &oa, const Eigen &m) { oa << m; } // everything!
void loadEigen(boost::archive::binary_iarchive &ia, Eigen &m) { ia >> m; } // everything!
#endif

// This approach is required for large problem sizes where an
// individual matrix could exceed 2GB. This has been an issue with
// serialization of matrices for MPI due to limited maximal MPI
// message size.
#ifdef LINEBYLINE
void saveMatrix(boost::archive::binary_oarchive &oa, const Matrix &m) {
  const size_t size1 = m.size1();
  const size_t size2 = m.size2();
  oa << size1 << size2;
  for (size_t i = 0; i < size1; i++) {
    ublas::vector<t_matel> vec = ublas::matrix_row<const Matrix>(m, i);
    oa << vec;
  }
}

void loadMatrix(boost::archive::binary_iarchive &ia, Matrix &m) {
  size_t size1, size2;
  ia >> size1 >> size2;
  m = Matrix(size1, size2);
  for (size_t i = 0; i < size1; i++) {
    ublas::vector<t_matel> vec;
    ia >> vec;
    ublas::matrix_row<Matrix>(m, i) = vec;
  }
}

void saveEigen(boost::archive::binary_oarchive &oa, const Eigen &m) {
  // RawEigen
  oa << m.value_orig;
  saveMatrix(oa, m.matrix);
  // Eigen
  oa << m.value_zero << m.nrpost;
  oa << m.absenergy << m.absenergyG << m.absenergyN;
}

void loadEigen(boost::archive::binary_iarchive &ia, Eigen &m) {
  // RawEigen
  ia >> m.value_orig;
  loadMatrix(ia, m.matrix);
  // Eigen
  ia >> m.value_zero >> m.nrpost;
  ia >> m.absenergy >> m.absenergyG >> m.absenergyN;
}
#endif

void saveRho(size_t N, const string &prefix, const DensMatElements &rho, const Params &P) {
  nrglog('H', "Storing density matrices [N=" << N << "]... ");
  my_assert(P.Ninit <= N && N <= P.Nmax - 1);
  const string fn = workdir.rhofn(prefix, N);
  std::ofstream MATRIXF(fn, ios::binary | ios::out);
  if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
  boost::archive::binary_oarchive oa(MATRIXF);
  size_t nr = rho.size();
  oa << nr;
  size_t cnt   = 0;
  size_t total = 0;
  for (const auto &[I, mat] : rho) {
    oa << I;
    saveMatrix(oa, mat);
    if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn));  // Check each time
    cnt++;
    total += mat.size1();
  }
  my_assert(cnt == nr);
  nrglog('H', "[total=" << total << " subspaces=" << nr << "]");
  MATRIXF.close();
}

DensMatElements loadRho(size_t N, const string &prefix, const Params &P) {
  nrglog('H', "Loading density matrices [N=" << N << "]...");
  my_assert(P.Ninit <= N && N <= P.Nmax - 1);
  DensMatElements rho;
  const string fn = workdir.rhofn(prefix, N);
  std::ifstream MATRIXF(fn, ios::binary | ios::in);
  if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
  boost::archive::binary_iarchive ia(MATRIXF);
  size_t nr;
  ia >> nr;
  size_t total = 0;
  for (size_t cnt = 0; cnt < nr; cnt++) {
    Invar inv;
    ia >> inv;
    loadMatrix(ia, rho[inv]);
    if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));  // Check each time
    total += rho[inv].size1();
  }
  nrglog('H', "[total=" << total << " subspaces=" << nr << "]");
  MATRIXF.close();
  if (P.removefiles)
    if (remove(fn)) throw std::runtime_error(fmt::format("Error removing {}", fn));
  return rho;
}

// save_transformations() stores all required information (energies, transformation matrices, subspace labels,
// dimensions of 'alpha' subspaces) that is needed to calculate reduced density matrix in the DM-NRG technique.
// Matrices are stored on disk as binary files.

// This function is called after the diagonalisation but prior to truncation (at iteration N) with DiagInfo diag,
// i.e. with *all calculated* (eigenvalue, eigenvector) pairs. NOTE: if diag=dsyevr, we effectively perform a
// truncation at the moment of the partial diagonalization!!

void save_transformations(size_t N, const DiagInfo &diag, const Params &P) {
  // P.Ninit-1 corresponds to the zero-th step, when diag contains the
  // eigenvalues from the initial diagonalization and the unitary matrices
  // are all identity matrices.
  if (!P.ZBW) {
    my_assert(N + 1 >= P.Ninit && N + 1 <= P.Nmax);
  } else
    my_assert(N == P.Ninit);
  nrglog('H', "Storing transformation matrices (N=" << N << ")...");
  const string fn = workdir.unitaryfn(N);
  ofstream MATRIXF(fn, ios::binary | ios::out);
  if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for writing.", fn));
  boost::archive::binary_oarchive oa(MATRIXF);
  size_t nr = diag.size();
  oa << nr;
  size_t cnt   = 0;
  size_t total = 0;
  for(const auto &[I, eig]: diag) {
    oa << I;
    saveEigen(oa, eig);
    if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error writing {}", fn)); // Check after each write.
    cnt++;
    total += eig.getnr();
  }
  my_assert(cnt == nr);
  nrglog('H', "[total=" << total << " subspaces=" << cnt << "]");
  MATRIXF.close();
}

void remove_transformation_files(size_t N, const Params &P) {
  remove(workdir.unitaryfn(N));
}

DiagInfo load_transformations(size_t N, const Params &P, bool remove_files = false) {
  DiagInfo diag;
  if (!P.ZBW) {
    my_assert(N + 1 >= P.Ninit && N + 1 <= P.Nmax);
  } else
    my_assert(N == P.Ninit);
  nrglog('H', "Loading transformation matrices (N=" << N << ")...");
  const string fn = workdir.unitaryfn(N);
  std::ifstream MATRIXF(fn, ios::binary | ios::in);
  if (!MATRIXF) throw std::runtime_error(fmt::format("Can't open file {} for reading", fn));
  boost::archive::binary_iarchive ia(MATRIXF);
  size_t nr; // Number of subspaces
  ia >> nr;
  size_t total = 0;
  for (size_t cnt = 0; cnt < nr; cnt++) {
    Invar inv;
    ia >> inv;
    loadEigen(ia, diag[inv]);
    if (MATRIXF.bad()) throw std::runtime_error(fmt::format("Error reading {}", fn));
    total += diag[inv].getnr();
  }
  nrglog('H', "[total=" << total << " subspaces=" << nr << "]");
  MATRIXF.close();
  if (remove_files) remove_transformation_files(N, P);
  return diag;
}

// Calculation of the contribution from subspace I1 of rhoN (density
// matrix at iteration N) to rhoNEW (density matrix at iteration N-1)
void cdmI(const size_t i,        // Subspace index (alpha=1,...,P.combs)
          const Invar &I1,       // Quantum numbers corresponding to subspace i
          const Matrix &rhoN,    // rho^N
          const Eigen &diagI1,   // contains U_{I1}
          Matrix &rhoNEW,        // rho^{N-1}
          const size_t N,
          const t_factor factor, // multiplicative factor that accounts for multiplicity
          const AllSteps &dm,
          const Params &P) 
{
  nrglog('D', "cdmI i=" << i << " I1=" << I1 << " factor=" << factor);
  // Range of indexes r and r' in matrix C^{QS,N}_{r,r'}, cf. Eq. (3.55)
  // in my dissertation.
  const size_t dim = rhoNEW.size2();
  my_assert(rhoNEW.size1() == rhoNEW.size2()); // quadratic matrix
  // number of states taken into account in the density-matrix at
  // *current* (Nth) stage (in subspace I1)
  const size_t nromega = rhoN.size2();
  my_assert(rhoN.size1() == rhoN.size2()); // quadratic matrix
  // continue only if connection exists
  if (nromega == 0 || dim == 0) return;
  // rmax (info[I1].rmax[i]) is the range of r in U^N_I1(omega|ri), only
  // those states that we actually kept..
  const size_t rmax = dm[N].at(I1).rmax.rmax(i);
  // rmax can be zero in the case a subspace has been completely truncated
  if (rmax == 0) return;
  // Otherwise, rmax must equal dim
  my_assert(rmax == dim);
  // Check range of omega: do the dimensions of C^N_I1(omega omega') and
  // U^N_I1(omega|r1) match?
  // diagI1.getnr() is the total number of states
  // obtained after diagonalisation. It can exceed nromega due to
  // truncation.
  // We do this test at this point, to ensure rmax!=0 and dim!=0
  // and nromega!=0, otherwise there is no contribution anyway.
  const size_t I1nr = diagI1.getnr();
  my_assert(nromega <= I1nr);
  // offset gives the offset that is added to r1,rp to find the
  // elements ri in U^N_I1(omega|ri)
  const size_t offset = dm[N].at(I1).rmax.offset(i);
  const size_t dim1   = diagI1.matrix.size1();
  const size_t dim2   = diagI1.matrix.size2();
  my_assert(nromega <= dim1 && offset + dim <= dim2);
  const ublas::matrix_range<const Matrix> U(diagI1.matrix, ublas::range(0, nromega), ublas::range(offset, offset + dim));
  Matrix T(dim, nromega);
  // T <- U^dag rhoN
  atlas::gemm(CblasConjTrans, CblasNoTrans, t_factor(1.0), U, rhoN, t_factor(0.0), T);
  // rhoNEW <- rhoNEW + factor T U
  // Note that we are *adding* to rhoNEW
  atlas::gemm(CblasNoTrans, CblasNoTrans, t_factor(factor), T, U, t_factor(1.0), rhoNEW);
}

// Calculation of the shell-N REDUCED DENSITY MATRICES:
// Calculate rho at previous iteration (N-1, rhoPrev) from rho at
// the current iteration (N, rho)

void calc_densitymatrix_iterN(const DiagInfo &diag,
                              const DensMatElements &rho, // input
                              DensMatElements &rhoPrev,   // output
                              size_t N, const AllSteps &dm, const Params &P) {
  nrglog('D', "calc_densitymatrix_iterN N=" << N);
  for (const auto &[I, dimsub] : dm[N - 1]) { // loop over all subspaces at *previous* iteration
    const InvarVec subs = Sym->dmnrg_subspaces(I);
    size_t dim          = dimsub.kept;
    rhoPrev[I]          = Matrix(dim, dim);
    if (!dim) continue;
    rhoPrev[I].clear();
    for (size_t i = 1; i <= P.combs; i++) {
      Invar sub = subs[i];
      const auto x = rho.find(sub);
      const auto y = diag.find(sub);
      if (x != rho.end() && y != diag.end()) {
        const t_factor coef = double(Sym->mult(sub)) / double(Sym->mult(I));
        cdmI(i, sub, x->second, y->second, rhoPrev[I], N, coef, dm, P);
      }
    }
  } // loop over invariant spaces
}

inline bool file_exists(const std::string &fn)
{
  ofstream F(fn, ios::binary | ios::out);
  return bool(F);
}

// Returns true if all the required density matrices are already
// saved on the disk.
bool already_computed(const std::string &prefix, const Params &P) {
  for (auto N = P.Nmax - 1; N > P.Ninit; N--) {
    const std::string fn = workdir.rhofn(prefix, N - 1); // note the minus 1
    if (!file_exists(fn)) {
      cout << fn << " not found. Computing." << endl;
      return false;
    }
  }
  return true;
}

/* calc_densitymatrix() is called prior to starting the NRG procedure for
 the second time. Here we calculate the shell-N density matrices for all
 iteration steps. */

void calc_densitymatrix(DensMatElements &rho, const AllSteps &dm, const Params &P, 
                        const std::string filename = FN_RHO) {
  if (P.resume && already_computed(filename, P)) {
    cout << "Not necessary: already computed!" << endl;
    return;
  }
  check_trace_rho(rho); // Must be 1.
  if (P.ZBW) return;
  TIME("DM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    cout << "[DM] " << N << endl;
    DiagInfo diag_loaded = load_transformations(N, P);
    DensMatElements rhoPrev;
    calc_densitymatrix_iterN(diag_loaded, rho, rhoPrev, N, dm, P);
    check_trace_rho(rhoPrev); // Make sure rho is normalized to 1.
    saveRho(N - 1, filename, rhoPrev, P);
    rho.swap(rhoPrev);
  }
}

// ****************** Calculation of the FULL REDUCED DENSITY MATIRICES

// Calculate rho(N), the shell-N density matrix, computed using
// the discarded states at shell N.
// Must be called AFTER calc_ZnD().
// Called fron nrg.cc immediately after the first NRG run (with N=Nmax-1),
// and also from calc_fulldensity_iterN (with lower N).
// A. Weichselbaum, J. von Delft, Phys. Rev. Lett. 99, 076402 (2007)
// T. A. Costi, V. Zlatic, Phys. Rev. B 81, 235127 (2010)
// H. Zhang, X. C. Xie, Q. Sun, Phys. Rev. B 82, 075111 (2010)
DensMatElements init_rho_FDM(size_t N, const AllSteps &dm, const Stats &stats, const Params &P) { // XXX: dm
  DensMatElements rhoFDM;
  double tr = 0.0;
  for (const auto &[I, ds] : dm[N]) {
    rhoFDM[I]     = Matrix(ds.max(), ds.max());
    rhoFDM[I].clear();
    Matrix &rhoI = rhoFDM[I];
    for (size_t i = ds.min(); i < ds.max(); i++) {
      const double betaE = ds.eig.absenergyN[i] / P.T;
      const double ratio = stats.wn[N] / stats.ZnDNd[N];
      double val2        = exp(-betaE) * ratio;
      val2               = std::isfinite(val2) ? val2 : 0.0;
      rhoI(i, i)         = val2;
      tr += Sym->mult(I) * val2;
    }
  }
  // Trace should be equal to the total weight of the shell-N contribution to the FDM.
  const double diff = (tr - stats.wn[N]) / stats.wn[N]; // relative error
  nrglog('w', "tr=" << tr << " diff=" << diff);
  if (std::isfinite(diff) && !num_equal(diff, 0.0, 1e-8))
    my_assert(stats.wn[N] < 1e-12);    // ..OK if small enough overall.
  return rhoFDM;
}

void calc_fulldensitymatrix_iterN(const Step &step, // only required for step::last()
                                  const DiagInfo &diag,
                                  const DensMatElements &rhoFDM, // input
                                  DensMatElements &rhoFDMPrev,   // output
                                  size_t N, const AllSteps &dm, const Stats &stats,
                                  const Params &P) {
  nrglog('D', "calc_fulldensitymatrix_iterN N=" << N);
  DensMatElements rhoDD;
  if (!step.last(N))
    rhoDD = init_rho_FDM(N, dm, stats, P);
  for (const auto &[I, ds] : dm[N - 1]) { // loop over all subspaces at *previous* iteration
    const InvarVec subs = Sym->dmnrg_subspaces(I);
    size_t dim          = ds.kept;
    rhoFDMPrev[I]       = Matrix(dim, dim);
    if (!dim) continue;
    rhoFDMPrev[I].clear();
    for (size_t i = 1; i <= P.combs; i++) {
      const auto sub = subs[i];
      // DM construction for non-Abelian symmetries: must include
      // the ratio of multiplicities as a coefficient.
      const t_factor coef = double(Sym->mult(sub)) / double(Sym->mult(I));
      // Contribution from the KK sector.
      const auto x1 = rhoFDM.find(sub);
      const auto y = diag.find(sub);
      if (x1 != rhoFDM.end() && y != diag.end())
        cdmI(i, sub, x1->second, y->second, rhoFDMPrev[I], N, coef, dm, P);
      // Contribution from the DD sector. rhoDD -> rhoFDMPrev
      if (!step.last(N)) {
        const auto x2 = rhoDD.find(sub);
        if (x2 !=rhoDD.end() && y != diag.end())
          cdmI(i, sub, x2->second, y->second, rhoFDMPrev[I], N, coef, dm, P);
      }
      // (Exception: for the N-1 iteration, the rhoPrev is already initialized with the DD sector of the last iteration.) }
    } // over combinations
  } // over subspaces
}

// Sum of statistical weights from site N to the end of the Wilson chain.
double sum_wn(size_t N, const Stats &stats, const Params &P) {
  double sum = 0.0;
  for (size_t n = P.Nmax - 1; n >= N; n--) sum += stats.wn[n];
  return sum;
}

void calc_fulldensitymatrix(const Step &step, DensMatElements &rhoFDM, const AllSteps &dm, const Stats &stats, const Params &P,
                            const std::string filename = FN_RHOFDM) {
  if (P.resume && already_computed(filename, P)) {
    cout << "Not necessary: already computed!" << endl;
    return;
  }
  if (P.ZBW) return;
  TIME("FDM");
  for (size_t N = P.Nmax - 1; N > P.Ninit; N--) {
    cout << "[FDM] " << N << endl;
    DiagInfo diag_loaded = load_transformations(N, P);
    DensMatElements rhoFDMPrev;
    calc_fulldensitymatrix_iterN(step, diag_loaded, rhoFDM, rhoFDMPrev, N, dm, stats, P);
    double tr       = trace(rhoFDMPrev);
    double expected = sum_wn(N, stats, P);
    double diff     = (tr - expected) / expected;
    nrglog('w', "tr[rhoFDM(" << N << ")]=" << tr << " sum(wn)=" << expected << " diff=" << diff);
    my_assert(num_equal(diff, 0.0));
    saveRho(N - 1, filename, rhoFDMPrev, P);
    rhoFDM.swap(rhoFDMPrev);
  }
}

#endif // _dmnrg_h_
