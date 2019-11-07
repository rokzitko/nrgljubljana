// param.cc - Parameter parsing
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _param_cc_
#define _param_cc_

// Abstraction type for (constant) values that may be initialized only once
// and then only read. [Non-strict version: allow rewriting by the same
// value.]
template<class T> class readonly 
{
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
      } else  // the value must not change!
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
  parambase(string keyword, string desc, string defaultv) :
    _keyword(keyword), _desc(desc), _value(defaultv) { };
  virtual ~parambase() { };
  virtual void setvalue_str(string newvalue) = 0;
  virtual void dump() = 0;
  string getkeyword() const { return _keyword; }
  string getdesc() const { return _desc; }
};

// Container for parameters
list<parambase*> allparams;

// Templated specialized classes for various storage types (int, double,
// string, bool)

template <class T>
class param : public parambase
{
private:
  T data;
  bool defaultval = true;
public:
  // Constructor: keyword is a CASE SENSITIVE name of the parameter,
  // desc is at this time used as in-line documentation and defaultv
  // is a string containing a default value which is immediately
  // parsed.
  param(const string &keyword, const string &desc, const string &defaultv) 
     : parambase(keyword, desc, defaultv) {
	data = fromstring<T>(_value);
	for(auto &i : allparams) {
	   if (i->getkeyword() == keyword)
	      my_error("Internal error: keyword conflict.");
	}	
	allparams.push_back((parambase*)this);
     }
  virtual void dump() { cout << _keyword << "=" << data << (!defaultval ? " *" : "") << endl; }
  // This line enables to access parameters using an object as a rvalue
  inline operator const T &() const { return data; }
  inline T value(void) const { return data; }
  virtual void setvalue_str(string newvalue) {
    _value = newvalue;
    data = fromstring<T>(newvalue);
    defaultval = false;
  }
  void setvalue(T newdata) { data = newdata; }
};

// CONVENTION: parameters that are user configurable are declared as
// param<T>, other parameters (set at runtime) are defined as basic
// types T.

enum dr_value { undefined, diagdsyev, diagdsyevr, diagzheev, diagzheevr };

namespace P {
  readonly<size_t> channels; // Number of channels
  readonly<size_t> coeffactor; // coefchannels = coeffactor * channels (typically = 1)
  readonly<size_t> coefchannels; // Number of coefficient sets (typically = channels)
  readonly<size_t> perchannel; // f-matrices per channel
  readonly<size_t> combs; // 4 for single-channel, 16 for two-channel, etc.
   
  /* Support for interleaved diagonalization */
  param<bool> substeps ("substeps", "Partial diagonalizations", "false");

  /* The number of iterations is determined either by Nmax or, if Tmin>0,
   by the lowest temperature(=energy) scale still considered, Tmin. In
   other words, if Tmin>0, the parameter Nmax is recomputed so that T at
   the Nmax-th iteration is equal or higher than the minimal temperature
   Tmin. See also parameter T below, which corresponds to the physical
   temperature for spectral function calculations, and parameter
   Tmin_ratio. NOTE: STAT::N runs from 0 to Nmax-1 (including). Note
   that coefficients run from 0 to Nmax (due to the zero-th site
   of the Wilson chain)! */
   
  size_t Nmax = 0; // Not readonly! 0 means non-initialized.
  bool ZBW = false; // Zero-bandwidth calculation if Nmax=Ninit.
  size_t Nlen = 0; // Nlen=Nmax for regular calculations. Nlen=1 for ZBW.
   // Nlen is length of wn, wnfactor, ZnD and dm.

  /* It is possible to include more than the zero-th site in the
     initial Wilson chain. This is controlled by parameter Ninit, which is
     equal to the last f operator index still retained. Ninit < Nmax. */
  param<size_t> Ninit ("Ninit", "Initial Wilson chain ops", "0");
   
  param<double> Lambda ("Lambda", "Logarithmic discretization parameter", "1.0");

  /* For "strategy=all", we recompute for the next iteration the matrix
   elements corresponding to all eigenpairs computed in the diagonalization
   step. For "strategy=kept", we only recompute those matrix elements that
   are actually required in the next iteration, i.e. those that correspond
   to the eigenpairs that we keep after truncation. Note 1: using
   "strategy=kept" also affects the calculation of thermodynamic and
   dynamic quantities, since the sum over all combinations of states is
   correspondingly truncated. However, since due to gain in efficiency with
   "strategy=kept" we can increase the number of states taken into account,
   this can be easily compensated since the contributions are weighted with
   exponentially decreasing Boltzmann factors. Note 2: when we do
   CFS calculations ("cfs=true"), where all eigenpairs need to be
   recalculated, the strategy is automatically switched to "all"! */
  param<string> strategy ("strategy", 
			  "Which vectors to use in recalculations?", 
			  "kept"); // Prior to 11.3.2010, the default was "all"
  param<bool> lastall ("lastall",
		       "Keep all states in last iteration",
		       "false");
   
   /* If lastalloverride=true, then "lastall" is not set to true automatically
    *      for CFS and FDM approaches. */
   param<bool> lastalloverride ("lastalloverride",
				"", "false");
   
  // cf. Yoshida, Whitaker, Oliveira PRB 41 9403 1990
  param<double> z ("z", "Parameter z in the logarithmic discretization", "1.0");

  // Valid discretization schemes: Y)oshida-Whitaker-Oliveira, C)ampo-Oliveira, Z)itko
  param<string> discretization ("discretization", 
				"Discretization approach used", "Y");

  // See initial.m for more information! Possibilities are "flat", "cosine"
  // and "dmft". The latter corresponds to arbitrary density of states
  // which requires special handling. NOTE: P::band is only used in the C++
  // part of the code for generating the header in the output files.
  param<string> band ("band", "Density of states in the band", "flat");
   
  // This options allows to reuse the existing discretization code
  // which assumes the conduction band to be contained within the
  // interval [-1:1], extending effectively this interval to
  // [-bandrescale:bandrescale]. This is implemented as follows:
  // 1. The band DOS is rescaled on the x-axis by 1/bandrescale and
  //    on the y-axis by bandrescale.
  // 2. The Wilson chain coefficients are scaled by bandrescale.
  // 3. In the data file, all eigenvalues reduced by the factor
  //    of bandrescale.
  // 4. In the NRG code, the current scale is multiplied
  //    by bandrescale, cf. function SCALE(N).
  // All input and output quantities can thus be in natural units.
  param<double> bandrescale ("bandrescale", "Automatic rescaling of the energies", "1.0");

  // \bar{\beta}, defines the effective temperature for computing
  // thermodynamic quantities at iteration N. See Krishna-Murthy p. 1009.
  // The default value is 1.0, which is somewhat large, but it is a safer
  // default compared to 0.46 that requires large truncation cutoff energy
  // to obtain converged results. The difference in extracted T_K is less
  // than 1% when comparing results obtained with betabar=1.0 and 0.75 with
  // keepenergy=10 and Lambda=2.
  param<double> betabar ("betabar", "\bar{\beta}", 
			 "1.0");

  /* The physical temperature for finite-T NRG calculations. See the
   discussion for parameter T_ratio, where a different default value is set
   unless this action is explicitly overriden by user. */
  param<double> T ("T", "Temperature, k_B T/D,", "-999.0");
  double Tpi; // T*pi

  /* Truncation parameters. If keepenergy>0.0, then the cut-off energy
   truncation scheme will be used. Parameter 'keep' is then used to cap the
   maximum number of states. If 'keepmin' is set to a positive value, it
   will ensure that a chosen minimal number of states is kept, in spite of
   the energy cut-off. */
   // TO DO: synonym "keepmax"
  param<size_t> keep ("keep", 
		      "Maximum number of states to keep at each stage", 
		      "100");
  param<double> keepenergy ("keepenergy", 
			    "The cut-off energy for truncation.", 
			    "-1.0");
  param<size_t> keepmin ("keepmin", 
			 "Minimum number of states to keep at each stage", 
			 "0");
   
  /* "NRG Ljubljana" allows fine-grained control over data logging. The
   following tokens are recognized:
  
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
  param<string> log ("log", "String with tokens that define what to log", "");
  param<bool> logall ("logall", "Log EVERYTHING.", "false");

  param<size_t> logenumber ("logenumber", "Nr. of eigenvalues to dump for 'e'", "10");

  // Broadening parameters. cf. cond-mat/0607497
  param<double> alpha ("alpha", "Width of logarithmic gaussian", "0.3");
  param<double> omega0 ("omega0", "Smallest energy scale in the problem", "-1.0");
  param<double> omega0_ratio ("omega0_ratio", "omega0 = omega0_ratio x T", "1.0");
   
  // P::bins is the number of bins per energy decade.
  param<size_t> bins ("bins", "Perform binning of spectral data: bins/decade", "100");
   
  // Overrides for the binning interval limits
  param<double> emin ("emin", "Lower binning limit", "-1.0");
  param<double> emax ("emax", "Upper binning limit", "-1.0");
   
  // If P::accumulation != 0, the accumulation point of the 
  // logarithmic mesh will be shifted away from omega=0. This may be
  // used in superconducting cases to get good spectral resolution
  // near the gap edges (cf. T. Hecht et al, JPCM).
  param<double> accumulation ("accumulation", "Shift of the accumulation points for binning", "0.0");
   
  // If P::linstep != 0 (and P::accumulation != 0), a linear mesh
  // will be used in the interval [-accumulation:accumulation]. This
  // may be used to properly resolve the sub-gap peaks (Shiba bound
  // states or resonances).
  param<double> linstep ("linstep", "Bin width for linear mesh", "0");
   
  // Peak trimming. One should better use conservative default values
  // here! The other two parameters control how we clip spectral
  // peaks with very small weight in order to gain on efficiency.
   
  // DISCARD_TRIM is relative value wrt "interval width" between
  // consecutive "bins". A "bin" is trimmed away if its weight is less
  // than (energy_bin(this)-energy_bin(next))*DISCARD_TRIM. This is
  // consistent with the fact that spectral density is spectral weight per
  // interval.
  param<double> DISCARD_TRIM ("discard_trim", 
			      "Peak clipping at the end of the run", "1e-16");

  // DISCARD_IMMEDIATELY is relative value wrt peak energy: we clip if
  // peak's weight is lower than energy*DISCARD_IMMEDIATELY.
  param<double> DISCARD_IMMEDIATELY ("discard_immediately", 
				     "Peak clipping on the fly", "1e-16");

  // Do N/N+1 patching (rather than N/N+2, which is the default)?
  param<bool> NN1 ("NN1", "Do N/N+1 patching?", "false");

  // In N/N+2 patching we can use the information from even iterations
  // (NN2even=true, which is the default) or from odd iterations
  // (NN2even=false).
  param<bool> NN2even ("NN2even", "Use even iterations in N/N+2 patching",
		       "true");
   
  // NN2avg=true is equivalent to averaging results from a NN2even=true
  // and a NN2even=false calculations.
  param<bool> NN2avg ("NN2avg", 
		      "Effectively averages even and odd N/N+2 spectra",
		      "false");
   
  // If NNtanh > 0, tanh window function is used instead of linear (i.e. a
  // smooth heap instead of a triangle).
  param<double> NNtanh ("NNtanh", "a in tanh[a(x-0.5)] window fnc", "0.0");
   
  param<double> goodE ("goodE", 
		       "Energy at which we calculate spectral density", "2.0");

  param<bool> savebins ("savebins", "Save binned (unbroadened) data", "false");
   
  param<bool> broaden ("broaden", "Enable broadening of spectra", "true");
  param<double> broaden_max ("broaden_max", 
			     "Parameters for producing smooth spectral density", "2.5");
  param<double> broaden_min ("broaden_min", "See broaden() for details", "-99.");
  param<double> broaden_min_ratio ("broaden_min_ratio", "See broaden() for details", "3.0");
  param<double> broaden_ratio ("broaden_ratio", "Broadening ratio", "1.05");
   
  param<bool> dumpenergies ("dumpenergies", "Dump (all) energies to file?", "false");
   
  /* If dumpannotated != 0, then a file "annotated.dat" is generated
   which contains 'dumpannotated' eigenvalue I lines for each
   iteration. Iterations are separated by blank lines. */
  param<size_t> dumpannotated ("dumpannotated", 
			      "Number of annotated eigenvalues to dump", "0");
  param<bool> dumpscaled ("dumpscaled",
			  "Use omega_N energy units in annotated.dat?", "true");

  // Use dumpabs=true together with dumpscaled=false to trace total
  // energies of excitations.
  param<bool> dumpabs ("dumpabs",
		       "Use total energies in annotated.dat?", "false");
  param<size_t> dumpprecision ("dumpprecision",
			    "Number of digits of precision in energy dumps", "8");
  param<bool> dumpgroups ("dumpgroups",
			  "Group degenerate states",
			  "true");
  param<double> grouptol ("grouptol",
			  "Energy tolerance for considering two states degenerate",
			  "1e-6");

  // Dump diagonal matrix elements of singlet operators. This is
  // particularly useful in the presence of the gap, where the matrix
  // elements are different in the long-chain limit.
  // "dumpdiagonal" matrix elements are saved for each singlet operator.
  // The output goes to the log file (i.e., standard output).
  param<size_t> dumpdiagonal ("dumpdiagonal", 
			      "Dump diagonal matrix elements of singlet ops.", 
			      "0");
   
  param<string> ops ("ops", "String with operators to be calculated", "");

  param<string> specs ("specs", "Spectral functions (singlet ops) to compute?", "");
  param<string> specd ("specd", "Spectral functions (doublet ops) to compute?", "");
  param<string> spect ("spect", "Spectral functions (triplet ops) to compute?", "");
  param<string> specq ("specq", "Spectral functions (quadruplet ops) to compute?", "");
  param<string> specot ("specot", "Spectral functions (orbital triplet ops) to compute?", "");
  
  // Output real and imaginary parts of calculated correlators (specs).
  // Be careful about the definitions, though!
  param<bool> reim ("reim", "Output imaginary parts of correlators?", "false");
   
  // Calculation of the temperature-dependent conductance G(T)
  param<string> specgt ("specgt", "Conductance curves to compute?", "");
  param<double> gtp ("gtp", "Parameter p for G(T) calculations", "0.7");

  // It is possible to compute the first and second moment of
  // A(w)(-df/dw), which are related to the thermopower and the heat
  // conductance. (These two use the same parameter 'gtp' as the
  // calculations of the temperature-dependent conductance G(T).)
  param<string> speci1t ("speci1t", "I_1 curves to compute?", "");
  param<string> speci2t ("speci2t", "I_2 curves to compute?", "");

  // Calculation of the temperature-depenedent susceptibility chi(T)
  param<string> specchit ("specchit", "Susceptibilities to compute?", "");
  param<double> chitp ("chitp", "Parameter p for chi(T) calculations", "1.0");

  // If chitp_ratio>0, chitp=chitp_ratio/betabar.
  param<double> chitp_ratio ("chitp_ratio", "Determine p from betabar", "-999");
   
  // 3-point vertex functions
  param<string> specv3 ("specv3", "3-leg vertex functions to compute?", "");
  param<bool> v3mm ("v3mm", "Compute 3-leg vertex on matsubara/matsubara axis?", "false");
   
  // If dm is set to true, density matrices are computed.
  param<bool> dm ("dm", "Compute density matrixes?", "false");
   
  // If finite=true, we calculate finite-temperature spectral
  // functions using Costi, Hewson, Zlastic approach.
  param<bool> finite ("finite", "Do a finite-temperature calculation?", "false");

  // Perform DM-NRG spectral function calculation. Cf. W. Hofstetter, Phys.
  // Rev. Lett. 85, 1508 (2000).
  param<bool> dmnrg ("dmnrg", "Do Hoffstetter's DM-NRG calculation?", "false");

  // If cfs is set to true, spectral functions are computed using the
  // complete basis set of the Wilson chain. See R. Peters, T. Pruschke, F.
  // B. Anders, Phys. Rev. B 74, 245113 (2006) 
  param<bool> cfs ("cfs", "Complete Fock space calculation?", "false");
   
  // Support for calculating greater and lesser correlation functions
  param<bool> cfsgt ("cfsgt", "Greater correlation function?", "false");
  param<bool> cfsls ("cfsls", "Lesser correlation function?", "false");

  // See A. Weichselbaum, J. von Delft, PRL 99 076402 (2007).
  param<bool> fdm ("fdm", "Full-density-matrix calculation?", "false");

  // Support for calculating greater and lesser correlation functions
  param<bool> fdmgt ("fdmgt", "Greater correlation function?", "false");
  param<bool> fdmls ("fdmls", "Lesser correlation function?", "false");
   
  // Calculate the expectation value <O>(T) using the FDM algorithm.
  param<bool> fdmexpv ("fdmexpv", "Expectation value using FDM", "false");
  param<size_t> fdmexpvn ("fdmexpvn", "At which iteration do we evaluate expv", "0");

  // Dynamical quantity calculations on the Mastubara axis
  param<bool> finitemats ("finitemats", "T>0 calculation on Matsubara axis", "false");
  param<bool> dmnrgmats ("dmnrgmats", "DMNRG calculation on Matsubara axis", "false");
  param<bool> fdmmats ("fdmmats", "FDM calculation on Matsubara axis", "false");
  param<size_t> mats ("mats", "Number of Matsubara points to collect", "100");

  // Formatting of output files
  param<size_t> width_td     ("width_td", "Widht of columns in 'td'", "16");
  param<size_t> width_custom ("width_custom", "Width of columns in 'custom'", "16");
  param<size_t> prec_td      ("prec_td", "Precision of columns in 'td'", "10");
  param<size_t> prec_custom  ("prec_custom", "Precision of columns in 'custom'", "10");

  param<size_t> prec_xy      ("prec_xy", "Precision of spectral function output", "10");

  // Safeguard feature: keep more states in case of a near degeneracy
  // Disabled if set to 0, otherwise keeps states that differ at most by P::safeguard
  // Prior to 3.11.2010, the default was 0, now it is finite, but small (1e-5).
  param<double> safeguard ("safeguard",
			   "Keep extra states in case of a near degeneracy", "1e-5");
  param<size_t> safeguardmax ("safeguardmax",
			      "Maximal number of additional states", "200");

  // Fix artificial splitting of states due to floating point round-off.
  // Note that the use of this feature with large value of epsilon may be
  // dangerous. (Default value of 1e-15 cures differences in the 2-3 least
  // significant bits.)
  param<double> fixeps ("fixeps", "Epsilon of splitting correction", "1e-15");
   
  param<string> diag ("diag", "Eigensolver routine", "default");
  dr_value diagroutine;

  param<double> diagratio ("diagratio", 
			   "Ratio of eigenvalues in partial diagonalisation", 
			   "1.0");
  param<size_t> dsyevrlimit ("dsyevrlimit", 
			     "Minimal size of matrix for dsyevr", 
			     "100");
  param<size_t> zheevrlimit ("zheevrlimit", 
			     "Minimal size of matrix for zheevr",
			     "100");
   
  /* If an insufficient number of states is computed during an iteration
   with partial diagonalisations, diagratio can be increased by a factor of
   restartfactor and calculation restarted if restart==true. */
  param<bool> restart ("restart", "Restart calculations?", "true");
  param<double> restartfactor ("restartfactor", "Multiply by?", "2.0");

   // It is possible that diagonalisation fails for some reason. If
   // checkvectors=true, additional tests are performed: finiteness
   // of all values, normalization and orthogonality tests.
  param<bool> checkdiag ("checkdiag", "Test diag results", "true");
   
   // Check if the density matrix for CFS has trace 1.
  param<bool> checkrho ("checkrho", "Test tr(rho)=1", "false");
   
  /* Usually all calculations (expectation values, spectral densities)
   are computed *after* each NRG iteration. This leaves out the very start
   of the procedure, before the first diagonalisation is performed. By
   setting calc0 to true, we compute all quantities also at the beginning.
  */
  param<bool> calc0 ("calc0", "Perform calculations at 0-th iteration?", "true");
   
  // Setting 'tdht' to an integer value >0 will perform the
  // calculations of thermodynamic quantities (output file td) at
  // temperatures above the bandwidth D. 
  param<int> tdht ("tdht", "Calculate TD properties for T>D", "0");

  /* Spin expressed in terms of the spin multiplicity, 2S+1. For SL code,
   P::spin is 1. Default value of 2 is valid for all other symmetry types.
  */
  param<size_t> spin ("spin", "Conduction-band electron spin (multiplicity)", "2");

  // old, orth or cpp. If tri=cpp, we do the tridiagonalisation in the C++
  // part of the code, see tridiag().
  param<string> tri ("tri", "Tridiagonalisation approach", "old");
  param<size_t> preccpp ("preccpp", "Arbitrary precision", "2000");

  param<bool> checksumrules ("checksumrules", "Check operator sumrules", "false");

  // Used only for symtype=QSZ and U1. IMPORTANT: this is not the correct way to
  // add a global magnetic field for the said symmetry types.
  param<double> globalB ("globalB", "Global magnetic field in z direction", "0.0");

  // Used only for symtype=U1. IMPORTANT: See the comment for globalB.
  param<double> globalBx ("globalBx", "Global magnetic field in x direction", "0.0");

   // If polarized=true, the DOS depends on the spin; this can be
   // used to describe ferromagnetic leads in a proper way. This
   // functionality can be used in the case of symtype=QSZ.
  param<bool> polarized ("polarized", "Spin-polarized conduction band", "false");
   
  // If pol2x2=true for symtype=U1, the DOS depends on two spin indexes,
  // i.e., it is described by a full 2x2 matrix in the spin space.
  param<bool> pol2x2 ("pol2x2", "2x2 structure in spin space", "false"); 

  // The following one is parsed from the command line! This is the
  // directory where we keep temporary files during the computation, so
  // workdir should point to some local storage device with ample space.
  string workdir = ".";
   
  // If NRG Ljubljana is compiled with MPI support, the default is
  // to use the MPI parallelization for spreading the diagonalizations
  // of the Hamiltonian matrices to all the cooperating processes.
#ifdef HAVE_BOOST_MPI
  param<bool> mpi ("mpi", "MPI parallelization", "true");
#else
  param<bool> mpi ("mpi", "MPI parallelization", "false");
#endif
         
  // Output imaginary parts of expectation values. Default is OFF! 
  param<bool> noimag ("noimag", "Do not output imaginary parts of expvs", "true");

  // Store information about subspaces (total, calculated and kept states).
  // For troubleshooting.
  param<bool> dumpsubspaces ("dumpsubspaces", "Save subspace info", "false");

  // Create file 'DONE' when finished. Useful for running on clusters.
  param<bool> done ("done", "Create DONE file?", "true");

  // Attempt to restart calculation by reading in the "unitary*"
  // files from a previous calculation. Automatically determines the
  // number of files.
  param<bool> resume ("resume", "Attempt restart?", "false");
  int laststored; // int, because -1 indicates that no stored data was found
   
  // For testing or for partial NRG calculations. If set to non-zero
  // value, the calculation is stopped at chosen step.
  param<int> forcestop ("forcestop", "Stop iteration?", "-1");

  // stopafter=nrg, stops calculation after the first sweep
  // stopafter=rho, stops calculation after computing the density matrix
  param<string> stopafter ("stopafter", 
			   "Stop calculation at some point?", "");

  // If set to false, the unitary transformation matrix and density
  // matrix files are kept after the calculation.
  param<bool> removefiles ("removefiles",
			   "Remove temporary data files?", "true");
   
  // If true, dumps the matrix elements of <f> channel by channel,
  // subspace pair by subspace pair.
  param<bool> dump_f ("dump_f", "Dump <f> matrix elements", "false");
   
  // support for channel-mixing terms in the Wilson chain
  param<bool> rungs ("rungs", "Channel-mixing terms in chain", "false");
   
  // Optimization in 3-pt vertex calculations: drop small terms.
  param<double> v3mmcutoff ("v3mmcutoff", "Cutoff for small terms", "1e-16");
   
  // Number of concurrent threads for matrix diagonalisation
  param<int> diagth ("diagth", "Diagonalisation threads", "1");
}

// Is i an allowed block index?
void allowed_block_index(size_t i)
{
   my_assert(1 <= i && i <= P::combs);
}

// Is ch an allowed channel index (0..P::channels-1)?
void allowed_channel(size_t ch)
{
   my_assert(ch < P::channels);
}

// Is ch an allowed index for coefficient table? There can be more
// coefficient tables than actual channels (for example in the case of
// spin-polarized conduction bands).
void allowed_coefchannel(size_t ch)
{
   my_assert(P::coefchannels >= P::channels);
   my_assert(ch < P::coefchannels);
}

// Returns true if any of the CFS or related spectral function
// calculations are requested.
bool cfs_flags(void)
{
   return (P::cfs || P::fdm);
}

// Calculate some constant parameters (invariants), etc.
// Called after parameters have been parsed and VALIDATED.
void calculate_invariants()
{
  // Take the first character (for backward compatibility)
  P::discretization.setvalue(string(P::discretization, 0, 1));
  if (P::chitp_ratio > 0.0)
     P::chitp.setvalue(P::chitp_ratio/P::betabar);
   P::Tpi = P::T * M_PI;
}

#endif // _param_cc_
