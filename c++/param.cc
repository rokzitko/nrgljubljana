#include <utility>



// param.cc - Parameter parsing
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _param_cc_
#define _param_cc_

// Abstraction type for (constant) values that may be initialized
// only once and then only read. [Non-strict version: allow rewriting
// by the same value using force_set().]
template <class T> class readonly {
  private:
  T data;
  bool init;

  public:
  readonly() { init = false; };
  inline operator const T &() const {
    my_assert(init);
    return data;
  }
  void operator=(T newvalue) {
    if (!init) {
      data = newvalue;
      init = true;
    } else // the value must not change!
      my_assert(data == newvalue);
  }
  // Don't abuse this one!
  void force_set(T newvalue) {
    data = newvalue;
    init = true;
  }
};

// Base class for parameter containers.
class parambase {
  protected:
  string _keyword;
  string _desc;
  string _value;

  public:
  parambase(string keyword, string desc, string defaultv) : _keyword(std::move(keyword)), _desc(std::move(desc)), _value(std::move(defaultv)){};
  virtual ~parambase()= default;
  virtual void setvalue_str(string newvalue) = 0;
  virtual void dump()                        = 0;
  string getkeyword() const { return _keyword; }
  string getdesc() const { return _desc; }
};

// Container for parameters
list<parambase *> allparams;

// Templated specialized classes for various storage types (int, double,
// string, bool)

template <class T> class param : public parambase {
  private:
  T data;
  bool defaultval = true;

  public:
  // Constructor: keyword is a CASE SENSITIVE name of the parameter,
  // desc is at this time used as in-line documentation and defaultv
  // is a string containing a default value which is immediately
  // parsed.
  param(const string &keyword, const string &desc, const string &defaultv) : parambase(keyword, desc, defaultv) {
    data = fromstring<T>(_value);
    for (auto &i : allparams) {
      if (i->getkeyword() == keyword) my_error("Internal error: keyword conflict.");
    }
    allparams.push_back((parambase *)this);
  }
  void dump() override { cout << _keyword << "=" << data << (!defaultval ? " *" : "") << endl; }
  // This line enables to access parameters using an object as a rvalue
  inline operator const T &() const { return data; }
  inline T value() const { return data; }
  void setvalue_str(string newvalue) override {
    _value     = newvalue;
    data       = fromstring<T>(newvalue);
    defaultval = false;
  }
  void setvalue(T newdata) { data = newdata; }
};

enum dr_value { undefined, diagdsyev, diagdsyevr, diagzheev, diagzheevr };

// CONVENTION: parameters that are user configurable are declared as
// param<T>, other parameters (set at runtime) are defined as basic
// types T.

// Interfacing:
// S = solver, C = constructor, L = low-level, * = hide (deprecated & experimental)
// //! lines define parameters that are only used in the high-level interface

namespace P {
  // *************************************************************
  // Parameters controlling discretization scheme and Wilson chain

  param<double> Lambda("Lambda", "Logarithmic discretization parameter", "2.0"); // S

  // Discretization scheme: Y)oshida-Whitaker-Oliveira, C)ampo-Oliveira, Z)itko-Pruschke
  param<string> discretization("discretization", "Discretization scheme", "Z"); // N

  // Twist parameter for the mesh. See Yoshida, Whitaker, Oliveira PRB 41 9403 1990
  param<double> z("z", "Parameter z in the logarithmic discretization", "1.0"); // N

  //! param<int> Nz ("Nz", "Number of discretization meshes", "1"); // S
  //! param<double> xmax ("xmax", "Largest x in discretization ODE solver", "20"); // S

  // Support of the hybridisation function,
  // [-bandrescale:bandrescale]. The automatic rescaling is
  // implemented as follows: 1. The spectral function is rescaled on
  // the x-axis by 1/bandrescale and on the y-axis by bandrescale. 2.
  // The Wilson chain coefficients are scaled by bandrescale. 3. In
  // the data file, all eigenvalues are reduced by the factor of
  // bandrescale. 4. In the NRG code, the current scale is multiplied
  // by bandrescale in function SCALE(N).
  param<double> bandrescale("bandrescale", "Band rescaling factor", "1.0"); // C

  // If polarized=true, the Wilson chain depends on the spin.
  param<bool> polarized("polarized", "Spin-polarized Wilson chain", "false"); // N

  // If pol2x2=true for symtype=U1, the hybridisation function (and
  // thus Wilson chain) depends on two spin indexes, i.e., it is
  // described by a full 2x2 matrix in the spin space.
  param<bool> pol2x2("pol2x2", "2x2 spin structure in Wilson chain", "false"); // N

  // Support for channel-mixing terms in the Wilson chain
  // (implemented for QS and QSZ symmetry types)
  param<bool> rungs("rungs", "Channel-mixing terms in Wilson chain", "false"); // N

  // Length of the Wilson chain. The number of iterations ise
  // determined either by Nmax or, if Tmin>0, by the lowest
  // temperature(=energy) scale still considered, Tmin. In other
  // words, if Tmin>0, the parameter Nmax is recomputed so that T at
  // the Nmax-th iteration is equal or higher than the minimal
  // temperature Tmin. NOTE: STAT::N runs from 0 to Nmax-1
  // (including), the coefficient indexes run from 0 to Nmax (due to
  // the inclusion of the zero-th site of the Wilson chain)!
  size_t Nmax = 0; // 0 means non-initialized.

  //! param<double> Tmin ("Tmin", "Lowest scale on the Wilson chain", "1e-4"); // S

  // If tri=cpp, we do the tridiagonalisation in the C++ part of the
  // code. In other cases,. we make use of external tools or Mathematica.
  param<string> tri("tri", "Tridiagonalisation approach", "old");               // N
  param<size_t> preccpp("preccpp", "Precision for tridiagonalisation", "2000"); // N

  // ************************
  // NRG iteration parameters

  param<string> diag("diag", "Eigensolver routine (dsyev|dsyevr|zheev|zheevr|default)", "default"); // N

  // For partial diagonalisation routines (dsyevr, zheevr), diagratio controls the fraction
  // of eigenspectrum that we compute.
  param<double> diagratio("diagratio", "Ratio of eigenstates computed in partial diagonalisation", "1.0"); // N
  param<size_t> dsyevrlimit("dsyevrlimit", "Minimal matrix size for dsyevr", "100");                       // N
  param<size_t> zheevrlimit("zheevrlimit", "Minimal matrix size for zheevr", "100");                       // N

  // If an insufficient number of states is computed during an
  // iteration with partial diagonalisations, diagratio can be
  // increased by a factor of restartfactor and calculation restarted
  // if restart=true.
  param<bool> restart("restart", "Restart calculation to achieve truncation goal?", "true"); // N
  param<double> restartfactor("restartfactor", "Rescale factor for restart=true", "2.0");    // N

  // Truncation parameters. If keepenergy>0.0, then the cut-off
  // energy truncation scheme will be used. Parameter 'keep' is then
  // used to cap the maximum number of states. If 'keepmin' is set to
  // a positive value, it will ensure that a chosen minimal number of
  // states is always kept.
  param<size_t> keep("keep", "Maximum number of states to keep at each step", "100");     // S
  param<double> keepenergy("keepenergy", "Cut-off energy for truncation", "-1.0");        // S
  param<size_t> keepmin("keepmin", "Minimum number of states to keep at each step", "0"); // S

  // Safeguard feature: keep more states in case of a near degeneracy
  // near the truncation cutoff. This prevents to truncate in the
  // middle of a (near-)degenerate cluster of states, leading to
  // unphysical effects. Disabled if set to 0, otherwise keeps states
  // that differ at most by P::safeguard.
  param<double> safeguard("safeguard", "Additional states to keep in case of a near degeneracy", "1e-5"); // N
  param<size_t> safeguardmax("safeguardmax", "Maximal number of additional states", "200");               // N

  // Fix artificial splitting of states due to floating point
  // round-off. Note that the use of this feature with large value of
  // epsilon may be dangerous. (The default conservative value of
  // 1e-15 cures differences in the 2-3 least significant bits.)
  param<double> fixeps("fixeps", "Threshold value for eigenvalue splitting corrections", "1e-15"); // N

  // ******************************************************
  // Physical temperature for finite-temperature quantities
  param<double> T("T", "Temperature, k_B T/D,", "0.001"); // S

  // \bar{\beta}, defines the effective temperature for computing the
  // thermodynamic quantities at iteration N. See Krishna-Murthy,
  // page 1009. The default value is 1.0, which is somewhat large,
  // but it is a safer default compared to 0.46 that requires large
  // truncation cutoff energy to obtain converged results.
  param<double> betabar("betabar", "Parameter \bar{\beta} for thermodynamics", "1.0"); // N

  // *************************************************************
  // Parameters controlling the problem and quantities of interest

  //! param<string> problem ("problem", "Model considered (templated)", "SIAM"); // C

  // List of operators being considered
  param<string> ops("ops", "Operators to be calculated", ""); // S

  // Dynamical quantities of interest
  param<string> specs("specs", "Spectral functions (singlet ops) to compute", "");           // S
  param<string> specd("specd", "Spectral functions (doublet ops) to compute", "");           // S
  param<string> spect("spect", "Spectral functions (triplet ops) to compute", "");           // S
  param<string> specq("specq", "Spectral functions (quadruplet ops) to compute", "");        // S
  param<string> specot("specot", "Spectral functions (orbital triplet ops) to compute", ""); // S

  // Calculation of the temperature-dependent conductance G(T) &
  // first and second moment of A(w)(-df/dw), which are related to
  // the thermopower and the heat conductance.
  param<string> specgt("specgt", "Conductance curves to compute", "");  // S
  param<string> speci1t("speci1t", "I_1 curves to compute", "");        // S
  param<string> speci2t("speci2t", "I_2 curves to compute", "");        // S
  param<double> gtp("gtp", "Parameter p for G(T) calculations", "0.7"); // N

  // Calculation of the temperature-depenedent susceptibility chi(T)
  param<string> specchit("specchit", "Susceptibilities to compute", "");      // S
  param<double> chitp("chitp", "Parameter p for chi(T) calculations", "1.0"); // N

  // If chitp_ratio>0, chitp=chitp_ratio/betabar.
  param<double> chitp_ratio("chitp_ratio", "Determine p from betabar", "-999"); // *

  // 3-point vertex functions
  param<string> specv3("specv3", "3-leg vertex functions to compute?", "");               // S
  param<bool> v3mm("v3mm", "Compute 3-leg vertex on matsubara/matsubara axis?", "false"); // S

  // **************
  // NRG algorithms

  // Calculate finite-temperature spectral functions using Costi,
  // Hewson, Zlastic approach.
  param<bool> finite("finite", "Perform Costi-Hewson-Zlatic finite-T calculation", "false"); // N

  // Perform DMNRG spectral function calculation. W. Hofstetter,
  // Phys. Rev. Lett. 85, 1508 (2000).
  param<bool> dmnrg("dmnrg", "Perform DMNRG (density-matrix NRG) calculation", "false"); // S

  // Perform complete-Fock-space calculation. R. Peters, T. Pruschke,
  // F. B. Anders, Phys. Rev. B 74, 245113 (2006)
  param<bool> cfs("cfs", "Perform CFS (complete Fock space) calculation", "false"); // S

  // Support for calculating greater and lesser correlation functions
  param<bool> cfsgt("cfsgt", "CFS greater correlation function", "false"); // N
  param<bool> cfsls("cfsls", "CFS lesser correlation function", "false");  // N

  // Perform full-density-matrix NRG calculation. Weichselbaum, J. von Delft, PRL 99 076402 (2007).
  param<bool> fdm("fdm", "Perform FDM (full-density-matrix) calculation", "false"); // S

  // Support for calculating greater and lesser correlation functions
  param<bool> fdmgt("fdmgt", "FDM greater correlation function?", "false"); // N
  param<bool> fdmls("fdmls", "FDM lesser correlation function?", "false");  // N

  // Calculate the expectation value <O>(T) using the FDM algorithm.
  param<bool> fdmexpv("fdmexpv", "Calculate expectation values using FDM", "false"); // S
  param<size_t> fdmexpvn("fdmexpvn", "Iteration where we evaluate expv", "0");       // N

  // Dynamical quantity calculations on the Mastubara axis
  param<bool> finitemats("finitemats", "T>0 calculation on Matsubara axis", "false"); // N
  param<bool> dmnrgmats("dmnrgmats", "DMNRG calculation on Matsubara axis", "false"); // S
  param<bool> fdmmats("fdmmats", "FDM calculation on Matsubara axis", "false");       // S
  param<size_t> mats("mats", "Number of Matsubara points to collect", "100");         // S

  // If dm is set to true, density matrices are computed.
  // Automatically enabled when needed (DMNRG, CFS, FDM).
  param<bool> dm("dm", "Compute density matrixes?", "false"); // N

  // **************
  // Frequency mesh
  param<double> broaden_max("broaden_max", "Broadening mesh maximum frequency", "10");             // N
  param<double> broaden_min("broaden_min", "Broadening mesh minimum frequency", "-99.");           // N
  param<double> broaden_min_ratio("broaden_min_ratio", "Auto-tune broaden_min parameter", "3.0");  // N
  param<double> broaden_ratio("broaden_ratio", "Common ration of the geometric sequence", "1.05"); // N

  //! param<double> mesh_max ("mesh_max", "Mesh maximum frequency", "10"); // C
  //! param<double> mesh_min ("mesh_min", "Mesh minimum frequency", "1e-4"); // C
  //! param<double> mesh_ratio ("mesh_ratio", "Common ratio of the geometric sequence", "1.05"); // C

  // Broadening parameters. cf. cond-mat/0607497
  param<double> alpha("alpha", "Width of logarithmic gaussian", "0.3");           // S
  param<double> omega0("omega0", "Smallest energy scale in the problem", "-1.0"); // N
  param<double> omega0_ratio("omega0_ratio", "omega0 = omega0_ratio x T", "1.0"); // N
  //! param<double> gamma ("gamma", "Parameter for Gaussian convolution step", "0.2"); // S

  // ******************************************************
  // Parameters for fine-grained control of NRG calculation

  // Number of concurrent threads for matrix diagonalisation
  param<int> diagth("diagth", "Diagonalisation threads", "1"); // N

  // Interleaved diagonalization
  param<bool> substeps("substeps", "Interleaved diagonalization", "false"); // N

  // For "strategy=all", we recompute for the next iteration the
  // matrix elements corresponding to all eigenpairs computed in the
  // diagonalization step. For "strategy=kept", we only recompute
  // those matrix elements that are actually required in the next
  // iteration, i.e. those that correspond to the eigenpairs that we
  // keep after truncation. Note 1: using "strategy=kept" also
  // affects the calculation of thermodynamic and dynamic quantities
  // (in non-CFS algorithms), since the sum over all combinations of
  // states is correspondingly truncated. However, since due to the
  // gain in efficiency with "strategy=kept" we can increase the
  // number of states taken into account, this can be easily
  // compensated since the contributions are weighted with
  // exponentially decreasing Boltzmann factors. Note 2: for CFS/FDM
  // calculations, where all eigenpairs need to be recalculated, the
  // strategy is automatically switched to "all"!
  param<string> strategy("strategy", "Recalculation strategy", "kept"); // N

  // It is possible to include more than the zero-th site in the
  // initial Wilson chain. This is controlled by parameter Ninit,
  // which is equal to the last f operator index still retained.
  // Ninit < Nmax.
  param<size_t> Ninit("Ninit", "Initial Wilson chain ops", "0"); // N

  // Output real and imaginary parts of calculated correlators (specs).
  param<bool> reim("reim", "Output imaginary parts of correlators?", "false"); // N

  // Number of eigenstates to save in "annotated.dat" per iteration.
  // Use dumpabs=true together with dumpscaled=false to trace total
  // energies of excitations.
  param<size_t> dumpannotated("dumpannotated", "Number of eigenvalues to dump", "0");                      // N
  param<bool> dumpabs("dumpabs", "Dump in terms of absolute energies", "false");                           // N
  param<bool> dumpscaled("dumpscaled", "Dump using omega_N energy units", "true");                         // N
  param<size_t> dumpprecision("dumpprecision", "Dump with # digits of precision", "8");                    // N
  param<bool> dumpgroups("dumpgroups", "Dump by grouping degenerate states", "true");                      // N
  param<double> grouptol("grouptol", "Energy tolerance for considering two states as degenerate", "1e-6"); // N

  // Dump diagonal matrix elements of singlet operators. This is
  // particularly useful in the presence of the gap, where the matrix
  // elements are different in the long-chain limit. "dumpdiagonal"
  // matrix elements are saved for each singlet operator. The output
  // goes to the standard output.
  param<size_t> dumpdiagonal("dumpdiagonal", "Dump diagonal matrix elements", "0"); // N

  // ********************
  // Binning & broadening

  param<bool> savebins("savebins", "Save binned (unbroadened) data", "false"); // N
  param<bool> broaden("broaden", "Enable broadening of spectra", "true");      // N

  // Overrides for the binning interval limits
  param<double> emin("emin", "Lower binning limit", "-1.0"); // N
  param<double> emax("emax", "Upper binning limit", "-1.0"); // N

  // Number of bins per energy decade for accumulating spectral data
  param<size_t> bins("bins", "bins/decade for spectral data", "1000"); // N

  // If P::accumulation != 0, the accumulation point of the
  // logarithmic mesh will be shifted away from omega=0. This may be
  // used in superconducting cases to get good spectral resolution
  // near the gap edges (cf. T. Hecht et al, JPCM).
  param<double> accumulation("accumulation", "Shift of the accumulation points for binning", "0.0"); // N

  // If P::linstep != 0 (and P::accumulation != 0), a linear mesh
  // will be used in the interval [-accumulation:accumulation]. This
  // may be used to properly resolve the sub-gap peaks (Shiba bound
  // states or resonances).
  param<double> linstep("linstep", "Bin width for linear mesh", "0"); // N

  // *************
  // Peak trimming

  // DISCARD_TRIM is relative value wrt "interval width" between
  // consecutive "bins". A "bin" is trimmed away if its weight is
  // less than (energy_bin(this)-energy_bin(next))*DISCARD_TRIM. This
  // is consistent with the fact that spectral density is spectral
  // weight per interval.
  param<double> discard_trim("discard_trim", "Peak clipping at the end of the run", "1e-16"); // N

  // DISCARD_IMMEDIATELY is relative value wrt peak energy: we clip
  // if peak's weight is lower than energy*DISCARD_IMMEDIATELY.
  param<double> discard_immediately("discard_immediately", "Peak clipping on the fly", "1e-16"); // N

  // Optimization in 3-pt vertex calculations: drop small terms.
  param<double> v3mmcutoff("v3mmcutoff", "Cutoff for small terms", "1e-16"); // *

  // ********
  // Patching

  // Parameter that controls the N/N+n patching schemes in non-FDM algorithms
  param<double> goodE("goodE", "Energy window parameter for patching", "2.0"); // N

  // Do N/N+1 patching, rather than N/N+2 which is the default
  param<bool> NN1("NN1", "Do N/N+1 patching?", "false"); // N

  // In N/N+2 patching we can use the information from even
  // iterations (NN2even=true, which is the default) or from odd
  // iterations (NN2even=false).
  param<bool> NN2even("NN2even", "Use even iterations in N/N+2 patching", "true"); // N

  // NN2avg=true is equivalent to averaging results from a
  // NN2even=true and a NN2even=false calculations.
  param<bool> NN2avg("NN2avg", "Average over even and odd N/N+2 spectra", "false"); // N

  // If NNtanh > 0, tanh window function is used instead of linear
  // (i.e. a smooth heap instead of a triangle).
  param<double> NNtanh("NNtanh", "a in tanh[a(x-0.5)] window function", "0.0"); // N

  // **************************
  // Formatting of output files
  param<size_t> width_td("width_td", "Widht of columns in 'td'", "16");               // N
  param<size_t> width_custom("width_custom", "Width of columns in 'custom'", "16");   // N
  param<size_t> prec_td("prec_td", "Precision of columns in 'td'", "10");             // N
  param<size_t> prec_custom("prec_custom", "Precision of columns in 'custom'", "10"); // N
  param<size_t> prec_xy("prec_xy", "Precision of spectral function output", "10");    // N

  // Checkpoint-restart functionality. Attempt to restart calculation
  // by reading the "unitary*" files from a previous calculation.
  // Automatically determines the number of files.
  param<bool> resume("resume", "Attempt restart?", "false"); // N
  int laststored;                                            // int, because -1 indicates that no stored data was found

  /* Fine-grained control over data logging with the following tokens:
   i - iteration (subspaces list)
   s - dump ancestor subspaces during Hamiltonian matrix building
   m - dump Hamiltonian matrix of each subspace [very verbose!]
   H - details about storing unitary matrices to files
   f - follow recalc_f() [low-level]
   F - matrix elements in recalc_f() [very verbose!]
   0 - functions calls in recalculations, etc. [high-level]
   r - follow recalc_general() [low-level]
   R - matrix elements in recalc_general() [very verbose!]
   g - follow calc_generic() [low-level]
   c - details about spectral function calculation [high-level]
   A - eigensolver diagnostics (routine used, matrix size)
   t - timing for eigensolver routines
   e - dump eigenvalues in function diagonalize_h()
   % - min and max energy (see also tdeA)
   d - eigenvalue computation [low-level]
   T - spectral peak trimming statistics
   * - merging details
   w - calculation of weights w_n
   X - show energy clusters in find_clusters()
   $ - Hamiltonian matrix sorting diagnostics
   M - MPI parallelization details
   ! - debug internal variables
   D - DMNRG calculation details
   @ - follow the program flow
   Useful combinations:
    fr - debug recalculation of irreducible matrix elements <||f||>
    ies - debug matrix construction and diagonalization
  */
  param<string> log("log", "list of tokens to define what to log", ""); // N
  param<bool> logall("logall", "Log everything", "false");              // N

  // ********************************************
  // Parameters where overrides are rarely needed

  param<bool> done("done", "Create DONE file?", "true");                         // N
  param<bool> calc0("calc0", "Perform calculations at 0-th iteration?", "true"); // N

  // If dmnrg=true, cfs=false, fdm=false, setting lastall=true will
  // enforce "keeping" all states in the last step of the NRG
  // iteration. For cfs=true and/or fdm=true, lastall is
  // automatically enable (unless this is overriden by setting
  // lastalloverride=true).
  param<bool> lastall("lastall", "Keep all states in the last iteratio for DMNRG", "false"); // N

  // If lastalloverride=true, then "lastall" is not set to true
  // automatically for CFS and FDM approaches.
  param<bool> lastalloverride("lastalloverride", "Override automatic lastall setting", "false"); // N

  // ***********************************
  // Deprecated (candidates for removal)

  // Used only for symtype=QSZ and U1. IMPORTANT: this is not the correct way to
  // add a global magnetic field for the said symmetry types.
  param<double> globalB("globalB", "Global magnetic field in z direction", "0.0"); // *

  // Used only for symtype=U1. IMPORTANT: See the comment for globalB.
  param<double> globalBx("globalBx", "Global magnetic field in x direction", "0.0"); // *

  // Setting 'tdht' to an integer value >0 will perform the
  // calculations of thermodynamic quantities (output file td) at
  // temperatures above the bandwidth D.
  param<int> tdht("tdht", "Calculate TD properties for T>D", "0"); // *

  // ******************************************************************
  // Parameters mostly relevant for low-level debugging and RG analysis

  // Store information about subspaces (total, calculated and kept states).
  param<bool> dumpsubspaces("dumpsubspaces", "Save detailed subspace info", "false"); // N

  // Enable dumps of the matrix elements of <f> channel by channel,
  // subspace pair by subspace pair.
  param<bool> dump_f("dump_f", "Dump <f> matrix elements", "false"); // N

  param<bool> dumpenergies("dumpenergies", "Dump (all) energies to file?", "false");  // N
  param<size_t> logenumber("logenumber", "# of eigenvalues to show for log=e", "10"); // N

  // stopafter=nrg, stops calculation after the first sweep
  // stopafter=rho, stops calculation after computing the density matrix
  param<string> stopafter("stopafter", "Stop calculation at some point?", ""); // N

  // For testing or for partial NRG calculations. If set to non-zero
  // value, the calculation is stopped at chosen step.
  param<int> forcestop("forcestop", "Stop iteration?", "-1"); // N

  // If set to false, the unitary transformation matrix and density
  // matrix files are kept after the calculation.
  param<bool> removefiles("removefiles", "Remove temporary data files?", "true"); // N

  // Output imaginary parts of expectation values. Default is OFF!
  param<bool> noimag("noimag", "Do not output imaginary parts of expvs", "true"); // N

  param<bool> checksumrules("checksumrules", "Check operator sumrules", "false"); // N

  // It is possible that diagonalisation fails for some reason. If
  // checkvectors=true, additional tests are performed: finiteness
  // of all values, normalization and orthogonality tests.
  param<bool> checkdiag("checkdiag", "Test diag results", "false"); // N

  // Check if the density matrix for CFS has trace 1.
  param<bool> checkrho("checkrho", "Test tr(rho)=1", "false"); // N

  // *******************************************
  // Internal parameters, not under user control

  readonly<size_t> channels;     // Number of channels
  readonly<size_t> coeffactor;   // coefchannels = coeffactor * channels (typically coeffactor=1)
  readonly<size_t> coefchannels; // Number of coefficient sets (typically coefchannels=channels)
  readonly<size_t> perchannel;   // f-matrices per channel (typically 1)
  readonly<size_t> combs;        // dimension of new shell Hilbert space, 4 for single-channel, 16 for two-channel, etc.
  bool ZBW    = false;           // Zero-bandwidth calculation if Nmax=Ninit.
  size_t Nlen = 0;               // Nlen=Nmax for regular calculations. Nlen=1 for ZBW.
  // Nlen is the length of wn, wnfactor, ZnD and dm vectors.
  double Tpi; // T*pi
  // Spin expressed in terms of the spin multiplicity, 2S+1. For SL &
  // SL 3 symmetry types, P::spin is 1. Default value of 2 is valid
  // for all other symmetry types.
  size_t spin;
  // Directory where we keep temporary files during the computation.
  // It should point to a storage device with ample space.
  string workdir = ".";
  dr_value diagroutine;
} // namespace P

// Is i an allowed block index?
void allowed_block_index(size_t i) { my_assert(1 <= i && i <= P::combs); }

// Is ch an allowed channel index (0..P::channels-1)?
void allowed_channel(size_t ch) { my_assert(ch < P::channels); }

// Is ch an allowed index for coefficient table? There can be more
// coefficient tables than actual channels (for example in the case of
// spin-polarized conduction bands).
void allowed_coefchannel(size_t ch) {
  my_assert(P::coefchannels >= P::channels);
  my_assert(ch < P::coefchannels);
}

// Returns true if any of the CFS or related spectral function
// calculations are requested.
bool cfs_flags() { return (P::cfs || P::fdm); }

// Calculate some constant parameters (invariants), etc.
// Called after parameters have been parsed and VALIDATED.
void calculate_invariants() {
  // Take the first character (for backward compatibility)
  P::discretization.setvalue(string(P::discretization, 0, 1));
  if (P::chitp_ratio > 0.0) P::chitp.setvalue(P::chitp_ratio / P::betabar);
  P::Tpi = P::T * M_PI;
}

#endif // _param_cc_
