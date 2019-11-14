// dmnrg.h - Density-matrix NRG
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _dmnrg_h_
#define _dmnrg_h_

// Routines for saving the density matrix to disk
string saverhofn(const string &prefix, int N) { return P::workdir + "/" + prefix + tostring(N); }

// Choose one!
#define LINEBYLINE
//#define WHOLEMATRIX

#ifdef WHOLEMATRIX
void saveMatrix(boost::archive::binary_oarchive &oa, const Matrix &m) { oa << m; }

void loadMatrix(boost::archive::binary_iarchive &ia, Matrix &m) { ia >> m; }

void saveEigen(boost::archive::binary_oarchive &oa, const Eigen &m) { oa << m; }

void loadEigen(boost::archive::binary_iarchive &ia, Eigen &m) { ia >> m; }
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
    matrix_row<const Matrix> mr = matrix_row<const Matrix>(m, i);
    ublas::vector<t_matel> vec  = ublas::vector<t_matel>(size2);
    vec                         = mr;
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
    matrix_row<Matrix>(m, i) = vec;
  }
}

void saveEigen(boost::archive::binary_oarchive &oa, const Eigen &m) {
  oa << m.nr << m.rmax << m.nrpost;
  oa << m.value << m.shift << m.absenergy;
  saveMatrix(oa, m.matrix0);
}

void loadEigen(boost::archive::binary_iarchive &ia, Eigen &m) {
  ia >> m.nr >> m.rmax >> m.nrpost;
  ia >> m.value >> m.shift >> m.absenergy;
  loadMatrix(ia, m.matrix0);
}
#endif

void saveRho_split(size_t N, const string &prefix, const DensMatElements &rho) {
  const string fn = saverhofn(prefix, N);
  std::ofstream MATRIXF(fn.c_str(), ios::binary | ios::out);
  if (!MATRIXF) my_error("Can't open file %s for writing.", fn.c_str());
  boost::archive::binary_oarchive oa(MATRIXF);
  size_t nr = rho.size();
  oa << nr;
  size_t cnt   = 0;
  size_t total = 0;
  for (const auto &i : rho) {
    const Invar inv = i.first;
    oa << inv;
    ostringstream fni;
    fni << fn << "-" << cnt;
    std::ofstream ONEMATRIX(fni.str().c_str(), ios::binary | ios::out);
    if (!ONEMATRIX) my_error("Can't open file %s for writing.", fni.str().c_str());
    boost::archive::binary_oarchive oaone(ONEMATRIX);
    saveMatrix(oaone, i.second);
    // Check each time..
    if (ONEMATRIX.bad()) my_error("Error writing %s", fni.str().c_str());
    ONEMATRIX.close();
    cnt++;
    total += i.second.size1();
  }
  my_assert(cnt == nr);
  nrglog('H', "[total=" << total << " nr subspaces=" << nr << "]");
  // Check once at the end.
  if (MATRIXF.bad()) my_error("Error writing %s", fn.c_str());
  MATRIXF.close();
}

void saveRho(size_t N, const string &prefix, const DensMatElements &rho) {
  my_assert(P::Ninit <= N && N <= P::Nmax - 1);
  nrglog('H', "Storing density matrices... ");
  saveRho_split(N, prefix, rho);
}

int remove(string filename) { return remove(filename.c_str()); }

void loadRho_split(size_t N, const string &prefix, DensMatElements &rho) {
  const string fn = saverhofn(prefix, N);
  std::ifstream MATRIXF(fn.c_str(), ios::binary | ios::in);
  if (!MATRIXF) my_error("Can't open file %s for reading", fn.c_str());
  boost::archive::binary_iarchive ia(MATRIXF);
  size_t nr;
  ia >> nr;
  size_t total = 0;
  for (size_t cnt = 0; cnt < nr; cnt++) {
    Invar inv;
    ia >> inv;
    ostringstream fni;
    fni << fn << "-" << cnt;
    std::ifstream ONEMATRIX(fni.str().c_str(), ios::binary | ios::in);
    if (!ONEMATRIX) my_error("Can't open file %s for reading.", fni.str().c_str());
    boost::archive::binary_iarchive iaone(ONEMATRIX);
    loadMatrix(iaone, rho[inv]);
    my_assert(rho[inv].size1() <= diag[inv].value.size());
    if (ONEMATRIX.bad()) my_error("Error reading %s", fni.str().c_str());
    ONEMATRIX.close();
    if (P::removefiles)
      if (remove(fni.str())) my_error("Error removing %s", fni.str().c_str());
    total += diag[inv].value.size();
  }
  nrglog('H', "[total=" << total << " nr subspaces=" << nr << "]");
  if (MATRIXF.bad()) my_error("Error reading %s", fn.c_str());
  MATRIXF.close();
  if (P::removefiles)
    if (remove(fn)) my_error("Error removing %s", fn.c_str());
}

void loadRho(size_t N, const string &prefix, DensMatElements &rho) {
  my_assert(P::Ninit <= N && N <= P::Nmax - 1);
  rho.clear();
  nrglog('H', "Loading density matrices...");
  loadRho_split(N, prefix, rho);
}

string unitaryfn(size_t N) { return P::workdir + "/" + FN_UNITARY + tostring(N); }

/* storetransformations() stores all required information (energies,
 transformation matrices, subspace labels, dimensions of 'alpha' subspaces)
 that is needed to calculate reduced density matrix in the DM-NRG
 technique. Matrices are stored on disk as binary files.

 This function is called after the diagonalisation but prior to truncation
 (at iteration N) with DiagInfo diag, i.e. with *all calculated*
 (eigenvalue, eigenvector) pairs. NOTE: if diag=dsyevr, we effectively
 perform a truncation at the moment of the partial diagonalization!! */

void store_transformations_split(size_t N, const DiagInfo &diag) {
  const string fn = unitaryfn(N);
  ofstream MATRIXF(fn.c_str(), ios::binary | ios::out);
  if (!MATRIXF) my_error("Can't open file %s for writing.", fn.c_str());
  boost::archive::binary_oarchive oa(MATRIXF);
  size_t nr = diag.size();
  oa << nr;
  size_t cnt   = 0;
  size_t total = 0;
  LOOP_const(diag, i) {
    const Invar I = INVAR(i);
    oa << I;
    ostringstream fni;
    fni << fn << "-" << cnt;
    ofstream ONEMATRIX(fni.str().c_str(), ios::binary | ios::out);
    if (!ONEMATRIX) my_error("Can't open file %s for writing.", fni.str().c_str());
    boost::archive::binary_oarchive oaone(ONEMATRIX);
    saveEigen(oaone, i.second);
    // Here check after each write.
    if (ONEMATRIX.bad()) my_error("Error writing %s", fni.str().c_str());
    ONEMATRIX.close();
    cnt++;
    total += i.second.value.size();
  }
  my_assert(cnt == nr);
  nrglog('H', "[total=" << total << " nr subspaces=" << cnt << "]");
  // Here make a single check at the end, not for each write..
  if (MATRIXF.bad()) my_error("Error writing %s", fn.c_str());
  MATRIXF.close();
}

void store_transformations(size_t N, const DiagInfo &diag) {
  // P::Ninit-1 corresponds to the zero-th step, when diag contains the
  // eigenvalues from the initial diagonalization and the unitary matrices
  // are all identity matrices.
  if (!P::ZBW) {
    my_assert(N + 1 >= P::Ninit && N + 1 <= P::Nmax);
  } else
    my_assert(N == P::Ninit);
  nrglog('H', "Storing transformation matrices...");
  store_transformations_split(N, diag);
}

void load_transformations_split(size_t N, DiagInfo &diag) {
  const string fn = unitaryfn(N);
  std::ifstream MATRIXF(fn.c_str(), ios::binary | ios::in);
  if (!MATRIXF) my_error("Can't open file %s for reading", fn.c_str());
  boost::archive::binary_iarchive ia(MATRIXF);
  size_t nr; // Number of subspaces
  ia >> nr;
  size_t total = 0;
  for (size_t cnt = 0; cnt < nr; cnt++) {
    Invar inv;
    ia >> inv;
    ostringstream fni;
    fni << fn << "-" << cnt;
    std::ifstream ONEMATRIX(fni.str().c_str(), ios::binary | ios::in);
    if (!ONEMATRIX) my_error("Can't open file %s for reading.", fni.str().c_str());
    my_assert(diag.count(inv) == 0); // added 18.8.2015
    boost::archive::binary_iarchive iaone(ONEMATRIX);
    loadEigen(iaone, diag[inv]);
    diag[inv].perform_checks(); // basic checks, added 9.11.2015
    if (ONEMATRIX.bad()) my_error("Error reading %s", fni.str().c_str());
    ONEMATRIX.close();
    total += diag[inv].value.size();
  }
  nrglog('H', "[total=" << total << " nr subspaces=" << nr << "]");
  if (MATRIXF.bad()) my_error("Error reading %s", fn.c_str());
  MATRIXF.close();
}

void load_transformations(size_t N, DiagInfo &diag) {
  if (!P::ZBW) {
    my_assert(N + 1 >= P::Ninit && N + 1 <= P::Nmax);
  } else
    my_assert(N == P::Ninit);
  diag.clear(); // blank state, no subspaces present
  nrglog('H', "Loading transformation matrices...");
  load_transformations_split(N, diag);
}

void remove_transformation_files_split(size_t N) {
  const string fn = unitaryfn(N);
  if (remove(fn)) my_error("Error removing %s", fn.c_str());
  for (size_t cnt = 0; cnt < diag.size(); cnt++) {
    ostringstream fni;
    fni << fn << "-" << cnt;
    if (remove(fni.str())) my_error("Error removing %s", fni.str().c_str());
  }
}

void remove_transformation_files(size_t N) {
  if (P::removefiles) remove_transformation_files_split(N);
}

// Calculation of the contribution from subspace I1 of rhoN (density
// matrix at iteration N) to rhoNEW (density matrix at iteration N-1)
void cdmI(size_t i,            // Subspace index (alpha=1,...,P::combs)
          const Invar &I1,     // Quantum numbers corresponding to subspace i
          const Matrix &rhoN,  // rho^N
          const Eigen &diagI1, // contains U_{I1}
          Matrix &rhoNEW,      // rho^{N-1}
          size_t N,
          t_factor factor) // multiplicative factor that accounts for multiplicity
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
  const size_t rmax = dm[N][I1].rmax.rmax(i);
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
  if (!(nromega <= I1nr)) {
    cout << "nromega=" << nromega << " I1nr=" << I1nr << " dim=" << dim << endl;
    my_error("Range mismatch in cdmI().");
  }
  // offset gives the offset that is added to r1,rp to find the
  // elements ri in U^N_I1(omega|ri)
  const size_t offset = dm[N][I1].rmax.offset(i);
  const size_t d1     = diagI1.matrix0.size1();
  const size_t d2     = diagI1.matrix0.size2();
  if (!(nromega <= d1 && offset + dim <= d2)) {
    cout << "nromega=" << nromega << " d1=" << d1 << endl;
    cout << "offset=" << offset << " dim=" << dim << " offset+dim=" << offset + dim << " d2=" << d2 << endl;
    my_error("Matrix dimensions errorin cdmI().");
  }
  const matrix_range<const Matrix> U(diagI1.matrix0, range(0, nromega), range(offset, offset + dim));
  Matrix T(dim, nromega);
  // T <- U^dag rhoN
  atlas::gemm(CblasConjTrans, CblasNoTrans, t_factor(1.0), U, rhoN, t_factor(0.0), T);
  // rhoNEW <- rhoNEW + factor T U
  // Note that we are *adding* to rhoNEW
  atlas::gemm(CblasNoTrans, CblasNoTrans, t_factor(factor), T, U, t_factor(1.0), rhoNEW);
}

InvarVec dmnrg_subspaces(const Invar &I) {
  InvarVec In = input_subspaces();
  for (size_t i = 1; i <= P::combs; i++) {
    In[i].inverse();
    In[i].combine(I);
  }
  return In;
}

// Calculation of the shell-N REDUCED DENSITY MATRICES:
// Calculate rho at previous iteration (N-1, rhoPrev) from rho at
// the current iteration (N, rho)

void calc_densitymatrix_iterN(DiagInfo &diag,
                              DensMatElements &rho,     // input
                              DensMatElements &rhoPrev, // output
                              size_t N) {
  nrglog('D', "calc_densitymatrix_iterN N=" << N);
  // loop over all subspaces at *previous* iteration
  for (const auto &ii : dm[N - 1]) {
    const Invar I     = ii.first;
    const InvarVec In = dmnrg_subspaces(I);
    size_t dim        = ii.second.kept;
    nrglog('D', "dm I=" << I << " dim=" << dim);
    rhoPrev[I] = Matrix(dim, dim);
    rhoPrev[I].clear();
    if (!dim) continue;
    for (size_t i = 1; i <= P::combs; i++) {
      const t_factor coef = double(mult(In[i])) / double(mult(I));
      nrglog('D', "i=" << i << " In[i]=" << In[i] << " coef=" << coef);
      const size_t cnt1 = rho.count(In[i]);
      const size_t cnt2 = diag.count(In[i]);
      if (cnt1 && cnt2) cdmI(i, In[i], rho[In[i]], diag[In[i]], rhoPrev[I], N, coef);
    }
  } // loop over invariant spaces
}

// Returns true if all the required density matrices are already
// saved on the disk.
bool already_computed(string prefix) {
  for (size_t N = P::Nmax - 1; N > P::Ninit; N--) {
    const string fn = saverhofn(prefix, N - 1); // note the minus 1
    size_t nrsubs   = dm[N].size();
    for (size_t i = 0; i <= nrsubs; i++) {
      ostringstream fni;
      fni << fn << "-" << i;
      ofstream F(fni.str().c_str(), ios::binary | ios::out);
      if (!F) {
        cout << fn << " not found. Computing." << endl;
        return false;
      }
    }
  }
  return true;
}

/* calc_densitymatrix() is called prior to starting the NRG procedure for
 the second time. Here we calculate the shell-N density matrices for all
 iteration steps. */

void calc_densitymatrix(DensMatElements &rho) {
  nrglog('@', "@ calc_densitymatrix");
  if (P::resume && already_computed(FN_RHO)) {
    cout << "Not necessary: already computed!" << endl;
    return;
  }
  check_trace_rho(rho); // Must be 1.
  if (P::ZBW) return;
  TIME("DM");
  for (size_t N = P::Nmax - 1; N > P::Ninit; N--) {
    cout << "[DM] " << N << endl;
    DiagInfo diag;
    load_transformations(N, diag);
    DensMatElements rhoPrev;
    calc_densitymatrix_iterN(diag, rho, rhoPrev, N);
    check_trace_rho(rhoPrev); // Make sure rho is normalized to 1.
    saveRho(N - 1, FN_RHO, rhoPrev);
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
void init_rho_FDM(DensMatElements &rhoFDM, size_t N) {
  nrglog('@', "@ init_rho_FDM(" << N << ")");
  double ZZ = STAT::ZnD[N];
  rhoFDM.clear();
  long double tr1 = 0.0;
  long double tr2 = 0.0;
  for (const auto &j : dm[N]) {
    const Invar I    = j.first;
    const size_t min = (LAST_ITERATION(N) ? 0 : j.second.kept);
    const size_t max = j.second.total;
    rhoFDM[I]        = Matrix(max, max);
    rhoFDM[I].clear();
    Matrix &rhoI = rhoFDM[I];
    for (size_t i = min; i < max; i++) {
      const long double Eabs = j.second.absenergy[i] - STAT::GSenergy;
      my_assert(Eabs >= 0.0);
      const long double betaE = Eabs / P::T;
      long double val1        = expl(-betaE) / ZZ;
      val1                    = std::isfinite(val1) ? val1 : 0.0;
      tr1 += mult(I) * val1;
      const long double ratio = STAT::wn[N] / ZZ;
      long double val2        = expl(-betaE) * ratio;
      val2                    = std::isfinite(val2) ? val2 : 0.0;
      rhoI(i, i)              = double(val2);
      tr2 += mult(I) * val2;
    }
  }
  nrglog('w', "tr1=" << tr1 << " tr2=" << tr2);
  // Sanity checks. The trace tr2 should be equal to the total
  // weight of the shell-N contribution to the full density matrix.
  if (std::isfinite(tr1) && !num_equal(tr1, 1.0)) my_warning("tr1=%24.16Lf", tr1); // This is not fatal! (but weird)
  const double diff = (tr2 - STAT::wn[N]) / STAT::wn[N];                           // relative error
  if (std::isfinite(diff) && !num_equal(diff, 0.0, 1e-8)) {
    my_warning("diff=%24.16lf", diff); // This is more serious..
    my_assert(STAT::wn[N] < 1e-12);    // ..but OK if small enough overall.
  }
}

void calc_fulldensitymatrix_iterN(DiagInfo &diag,
                                  DensMatElements &rhoFDM,     // input
                                  DensMatElements &rhoFDMPrev, // output
                                  size_t N) {
  nrglog('D', "calc_fulldensitymatrix_iterN N=" << N);
  DensMatElements rhoDD;
  if (!LAST_ITERATION(N)) init_rho_FDM(rhoDD, N);
  // loop over all subspaces at *previous* iteration
  for (const auto &ii : dm[N - 1]) {
    const Invar I     = ii.first;
    const InvarVec In = dmnrg_subspaces(I);
    size_t dim        = ii.second.kept;
    nrglog('D', "fdm I=" << I << " dim=" << dim);
    rhoFDMPrev[I] = Matrix(dim, dim);
    rhoFDMPrev[I].clear();
    if (!dim) continue;
    for (size_t i = 1; i <= P::combs; i++) {
      // DM construction for non-Abelian symmetries: must include
      // the ratio of multiplicities as a coefficient.
      const t_factor coef = double(mult(In[i])) / double(mult(I));
      // Contribution from the KK sector. (Exception: for the N-1 iteration,
      // the rhoPrev is already initialized with the DD sector of the last
      // iteration.)
      nrglog('D', "KK sector");
      cdmI(i, In[i], rhoFDM[In[i]], diag[In[i]], rhoFDMPrev[I], N, coef);
      // Contribution from the DD sector. rhoDD -> rhoFDMPrev
      if (!LAST_ITERATION(N)) {
        nrglog('D', "DD sector");
        cdmI(i, In[i], rhoDD[In[i]], diag[In[i]], rhoFDMPrev[I], N, coef);
      }
    }
  } // loop over invariant spaces
}

// Sum of statistical weights from site N to the end of the Wilson chain.
double sum_wn(size_t N) {
  double sum = 0.0;
  for (size_t n = P::Nmax - 1; n >= N; n--) sum += STAT::wn[n];
  return sum;
}

void calc_fulldensitymatrix(DensMatElements &rhoFDM) {
  nrglog('@', "@ calc_densitymatrix");
  if (P::resume && already_computed(FN_RHOFDM)) {
    cout << "Not necessary: already computed!" << endl;
    return;
  }
  if (P::ZBW) return;
  TIME("FDM");
  for (size_t N = P::Nmax - 1; N > P::Ninit; N--) {
    cout << "[FDM] " << N << endl;
    DiagInfo diag;
    load_transformations(N, diag);
    DensMatElements rhoFDMPrev;
    calc_fulldensitymatrix_iterN(diag, rhoFDM, rhoFDMPrev, N);
    double tr       = trace(rhoFDMPrev);
    double expected = sum_wn(N);
    double diff     = (tr - expected) / expected;
    nrglog('w', "tr[rhoFDM(" << N << ")]=" << tr << " sum(wn)=" << expected << " diff=" << diff);
    my_assert(num_equal(diff, 0.0));
    saveRho(N - 1, FN_RHOFDM, rhoFDMPrev);
    rhoFDM.swap(rhoFDMPrev);
  }
}

#endif // _dmnrg_h_
