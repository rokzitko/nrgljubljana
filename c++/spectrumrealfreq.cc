// Real-frequency spectral function
class SpectrumRealFreq : public Spectrum {
 private:
   Bins fspos, fsneg; // Full spectral information
   void mergeNN2half(Bins &fullspec, const Bins &cs, const Step &step);
   void mergeNN2(spCS_t, const Step &);
   void mergeCFS(spCS_t);
   void weight_report(const double imag_tolerance = 1e-10);
   void trim();
   void savebins();
   void continuous();
 public:
   SpectrumRealFreq(const string &opname, const string &filename, shared_ptr<Algo> algotype, const Params &P) :
     Spectrum(opname, filename, algotype, P), fspos(P), fsneg(P) {};
   void merge(spCS_t, const Step &step) override;
   ~SpectrumRealFreq() override;
};

SpectrumRealFreq::~SpectrumRealFreq() {
  cout << "Spectrum: " << opname << " " << algotype->name() << " ->"; // appended in savebins() & continuous()
  trim();
  savebins();
  continuous();
  weight_report();
}

// Merge the spectrum for a finite Wilson chain into the "true" NRG spectrum. For complete Fock space NRG
// calculation, the merging tricks are not necessary: we just collect all the delta peaks from all iterations.
void SpectrumRealFreq::merge(spCS_t cs, const Step &step) {
  if (algotype->merge() == "NN2") return mergeNN2(cs, step);
  if (algotype->merge() == "CFS") return mergeCFS(cs);
  my_assert_not_reached();
}

// Spectrum merging for complete Fock space calculation.
void SpectrumRealFreq::mergeCFS(spCS_t cs) {
  auto csb = dynamic_pointer_cast<ChainSpectrumBinning>(cs);
  fspos.merge(csb->spos);
  fsneg.merge(csb->sneg);
}

// Window functions: f(0)=0, f(1)=1.
inline double fnc_linear(const double x) { return x; }

inline double fnc_tanh_0(const double x, const double NNtanh) { return tanh(NNtanh * (x - 0.5)); }

inline double fnc_tanh(const double x, const double NNtanh) {
  const double f0 = fnc_tanh_0(0, NNtanh);
  const double fx = fnc_tanh_0(x, NNtanh);
  const double f1 = fnc_tanh_0(1, NNtanh);
  return (fx - f0) / (f1 - f0);
}

inline double fnc(double x, const Params &P) {
  return P.NNtanh > 0 ? fnc_tanh(x, P.NNtanh) : fnc_linear(x);
}

inline double windowfunction(const double E, const double Emin, const double Ex, const double Emax, const Step &step, const Params &P) {
  // Exception 1
  if (E <= Ex && step.last()) return 1.0;
  // Exception 2
  if (E >= Ex && step.first()) return 1.0;
  // Exception 3
  if (P.ZBW) return 1.0;
  if (E <= Emin || E >= Emax) return 0.0;
  if (Emin < E && E <= Ex) return fnc((E - Emin) / (Ex - Emin), P);
  if (Ex < E && E < Emax) return 1.0 - fnc((E - Ex) / (Emax - Ex), P);
  my_assert_not_reached();
}

// Here we perform the actual merging of data using the N/N+2 scheme.
// Note that we use a windowfunction (see above) to accomplish the
// smooth combining of data.
void SpectrumRealFreq::mergeNN2half(Bins &fullspec, const Bins &cs, const Step &step) {
  double Emin = step.scale() * P.getEmin(); // p
  double Ex   = step.scale() * P.getEx();   // p Lambda
  double Emax = step.scale() * P.getEmax(); // p Lambda^2
  if (P.ZBW) {                              // override for zero bandwidth calculation
    Emin = 0;
    Emax = std::numeric_limits<double>::max(); // infinity
  }
  const auto len = fullspec.bins.size();
  my_assert(len == cs.bins.size()); // We require equivalent bin sets!!
  for (const auto i: range0(len)) {
    const double energy   = cs.bins[i].first;
    const t_weight weight = cs.bins[i].second;
    if (Emin < energy && energy < Emax && weight != 0.0) {
      const double factor = (P.NN2avg ? 0.5 : 1.0);
      fullspec.bins[i].second += factor * weight * windowfunction(energy, Emin, Ex, Emax, step, P);
    }
  }
}

// Spectrum merging using the N/N+n patching. One has to be careful about
// the specifics, in particular the values of Emin, Ex and Emax. The
// current choice seems to be working quite all right.
// See R. Bulla, T. A. Costi, D. Vollhardt, Phys. Rev. B 64, 045103 (2001).

void SpectrumRealFreq::mergeNN2(spCS_t cs, const Step &step) {
  auto csb = dynamic_pointer_cast<ChainSpectrumBinning>(cs);
  if (!step.N_for_merging()) return;
  mergeNN2half(fspos, csb->spos, step);
  mergeNN2half(fsneg, csb->sneg, step);
}

void SpectrumRealFreq::weight_report(const double imag_tolerance) {
  auto fmt = [imag_tolerance](t_weight x) -> string { return abs(x.imag()) < imag_tolerance ? to_string(x.real()) : to_string(x); };
  const t_weight twneg = fsneg.total_weight();
  const t_weight twpos = fspos.total_weight();
  cout << endl << "[" << opname << "]"
       << " pos=" << fmt(twpos) << " neg=" << fmt(twneg) << " sum= " << fmt(twpos + twneg) << endl;
  for (int m = 1; m <= 4; m++) {
    const auto mom = moment(fsneg.bins, fspos.bins, m);   // Spectral moments from delta-peaks
    cout << fmt::format("mu{}={} ", m, fmt(mom));
  }
  std::cout << std::endl;
  const t_weight f = fd_fermi(fsneg.bins, fspos.bins, P.T);
  const t_weight b = fd_bose(fsneg.bins, fspos.bins, P.T);
  cout << "f=" << fmt(f) << " b=" << fmt(b) << endl;
}

// Save binary raw (binned) spectral function. If using complex numbers and P.reim==true, we save triplets
// (energy,real part,imag part).
void SpectrumRealFreq::savebins() {
  if (!P.savebins) return;
  const std::string fn = filename + ".bin";
  cout << " " << fn;
  ofstream Fbins = safe_open(fn, true); // true=binary!
  for (const auto &[e, w] : fspos.bins) {
    const auto [wr, wi] = reim(w);
    my_assert(e > 0.0);
    Fbins.write((char *)&e, sizeof(double));
    Fbins.write((char *)&wr, sizeof(double));
    if (P.reim)
      Fbins.write((char *)&wi, sizeof(double));
  }
  for (const auto &[abse, w] : fsneg.bins) {
    const auto [wr, wi] = reim(w);
    const double e = -abse;
    my_assert(e < 0.0); // attention!
    Fbins.write((char *)&e, sizeof(double));
    Fbins.write((char *)&wr, sizeof(double));
    if (P.reim) 
      Fbins.write((char *)&wi, sizeof(double));
  }
}

void SpectrumRealFreq::trim() {
  fspos.trim();
  fsneg.trim();
}

// Energy mesh for spectral functions
std::vector<double> make_mesh(const Params &P) {
  const double broaden_min = P.get_broaden_min();
  std::vector<double> vecE; // Energies on the mesh
  for (double E = P.broaden_max; E > broaden_min; E /= P.broaden_ratio) vecE.push_back(E);
  return vecE;
}

void SpectrumRealFreq::continuous() {
  if (!P.broaden) return;
  const double alpha  = P.alpha;
  const double omega0 = P.omega0 < 0.0 ? P.omega0_ratio * P.T : P.omega0;
  Spikes densitypos, densityneg;
  const std::vector<double> vecE = make_mesh(P); // Energies on the mesh
  for (const double E : vecE) {
    weight_bucket valpos, valneg;
    for (const auto &[e, w] : fspos.bins) {
      my_assert(e > 0.0);
      valpos += w * BR_NEW(E, e, alpha, omega0);
      valneg += w * BR_NEW(-E, e, alpha, omega0);
    }
    for (const auto &[e, w] : fsneg.bins) {
      my_assert(e > 0.0); // attention!
      valneg += w * BR_NEW(-E, -e, alpha, omega0);
      valpos += w * BR_NEW(E, -e, alpha, omega0);
    }
    densitypos.emplace_back(E, valpos);
    densityneg.emplace_back(-E, valneg);
  }
  ranges::sort(densityneg, sortfirst());
  ranges::sort(densitypos, sortfirst());
  const std::string fn = filename + ".dat";
  cout << " " << fn;
  ofstream Fdensity = safe_open(fn);
  densityneg.save(Fdensity, P.prec_xy, P.reim);
  densitypos.save(Fdensity, P.prec_xy, P.reim);
}
