#ifndef _spectrum_hpp_
#define _spectrum_hpp_

#include <algorithm>
#include "traits.hpp"
#include "params.hpp"
#include "bins.hpp"
#include "matsubara.hpp"

namespace NRG {

template<scalar S, typename t_weight = weight_traits<S>>
class ChainBinning {
 private:
   const Params &P;
   Bins<S> spos, sneg;
 public:
   explicit ChainBinning(const Params &P) : P(P), spos(P), sneg(P) {}
   void add(const double energy, const t_weight weight) {
     if (energy >= 0.0)
       spos.add(energy, weight);
     else
       sneg.add(-energy, weight);
   }
   void add(const std::pair<double, t_weight> &p, const t_weight factor) {
      const auto & [energy, weight] = p;
      this->add(energy, factor * weight);
   }
   auto total_weight() const { return spos.total_weight() + sneg.total_weight(); }
   template<scalar U> friend class SpectrumRealFreq;
};

template<scalar S, typename t_weight = weight_traits<S>>
class ChainMatsubara {
 private:
   const Params &P;
   Matsubara<S> m;
 public:
   explicit ChainMatsubara(const Params &P, const gf_type gt) : P(P), m(P.mats, gt, P.T){};
   void add(const size_t n, const t_weight w) { m.add(n, w); }
   template<scalar U> friend class GFMatsubara;
};

template<scalar S, typename t_weight = weight_traits<S>>
class ChainTempDependence {
 private:
   const Params &P;
   Temp<S> v;
 public:
   explicit ChainTempDependence(const Params &P) : P(P), v(P) {}
   void add(const double T, const t_weight value) { v.add_value(T, value); }
   template<scalar U> friend class TempDependence;
};

// Real-frequency spectral function
template<scalar S>
class SpectrumRealFreq {
 private:
   const std::string name, algoname, filename; // e.g. "A_d-A_d", "FT", "spec_A_d-A_d_dens_FT.dat"
   const Params &P;
   Bins<S> fspos, fsneg; // Full spectral information, separately for positive and negative frequencies
   void mergeNN2half(Bins<S> &fullspec, const Bins<S> &cs, const Step &step);
   void weight_report(const double imag_tolerance = 1e-10);
   void trim() {
     fspos.trim();
     fsneg.trim();
   }
   void savebins();
   void continuous();
 public:
   SpectrumRealFreq(const std::string &name, const std::string &algoname, const std::string &filename, const Params &P) :
     name(name), algoname(algoname), filename(filename), P(P), fspos(P), fsneg(P) {}
   void mergeCFS(const ChainBinning<S> &cs) {
     fspos.merge(cs.spos); // Collect delta peaks
     fsneg.merge(cs.sneg);
   }
   void mergeNN2(const ChainBinning<S> &cs, const Step &step) {
     if (!step.N_for_merging()) return;
     mergeNN2half(fspos, cs.spos, step); // Spectrum merging using the N/N+n patching.
     mergeNN2half(fsneg, cs.sneg, step);
   }
   void save() {
     fmt::print(fmt::emphasis::bold, "Spectrum: {} {} -> ", name, algoname); // savebins() & continuous() append the filenames
     trim();
     if (P.savebins) savebins();
     if (P.broaden) continuous();
     weight_report();
   }
};

inline double windowfunction(const double E, const double Emin, const double Ex, const double Emax, const Step &step, const Params &P) {
  if (E <= Ex && step.last()) return 1.0;  // Exception 1
  if (E >= Ex && step.first()) return 1.0; // Exception 2
  if (P.ZBW()) return 1.0;                   // Exception 3
  if (E <= Emin || E >= Emax) return 0.0;  // Optimization
  // Window functions: f(0)=0, f(1)=1.
  const auto fnc_linear = [](const auto x) { return x ; };
  const auto fnc_tanh_0 = [](const double x, const double NNtanh) { return tanh(NNtanh * (x - 0.5)); };
  const auto fnc_tanh = [fnc_tanh_0](const double x, const double NNtanh) {
    const auto f0 = fnc_tanh_0(0, NNtanh);
    const auto fx = fnc_tanh_0(x, NNtanh);
    const auto f1 = fnc_tanh_0(1, NNtanh);
    return (fx - f0) / (f1 - f0);
  };
  const auto fnc = [&P, fnc_linear, fnc_tanh](const auto x) { return P.NNtanh > 0 ? fnc_tanh(x, P.NNtanh) : fnc_linear(x); };
  if (Emin < E && E <= Ex) return fnc((E - Emin) / (Ex - Emin));
  if (Ex < E && E < Emax) return 1.0 - fnc((E - Ex) / (Emax - Ex));
  my_assert_not_reached();
}

// Here we perform the actual merging of data using the N/N+2 scheme. Note that we use a windowfunction (see above)
// to accomplish the smooth combining of data.
// See R. Bulla, T. A. Costi, D. Vollhardt, Phys. Rev. B 64, 045103 (2001)
template<scalar S>
void SpectrumRealFreq<S>::mergeNN2half(Bins<S> &fullspec, const Bins<S> &cs, const Step &step) {
  const auto Emin = step.scale() * P.getEmin(); // p
  const auto Ex   = step.scale() * P.getEx();   // p Lambda
  const auto Emax = step.scale() * P.getEmax(); // p Lambda^2
  const auto len = fullspec.bins.size();
  my_assert(len == cs.bins.size()); // We require equivalent bin sets!!
  for (const auto i: range0(len)) {
    const auto [energy, weight] = cs.bins[i];
    if (Emin < energy && energy < Emax && weight != 0.0) {
      const auto factor = P.NN2avg ? 0.5 : 1.0;
      fullspec.bins[i].second += factor * weight * windowfunction(energy, Emin, Ex, Emax, step, P);
    }
  }
}

template<scalar S>
void SpectrumRealFreq<S>::weight_report(const double imag_tolerance) {
  const auto fmt = [imag_tolerance](const auto x) -> std::string { return abs(x.imag()) < imag_tolerance ? to_string(x.real()) : to_string(x); };
  const auto twneg = fsneg.total_weight();
  const auto twpos = fspos.total_weight();
  std::cout << std::endl << "pos=" << fmt(twpos) << " neg=" << fmt(twneg) << " sum= " << fmt(twpos + twneg) << std::endl;
  for (int m = 1; m <= 4; m++) {
    const auto mom = moment(fsneg.bins, fspos.bins, m);   // Spectral moments from delta-peaks
    std::cout << fmt::format("mu{}={} ", m, fmt(mom));
  }
  std::cout << std::endl;
  const auto f = fd_fermi(fsneg.bins, fspos.bins, P.T);
  const auto b = fd_bose (fsneg.bins, fspos.bins, P.T);
  std::cout << "f=" << fmt(f) << " b=" << fmt(b) << std::endl;
}

// Save binary raw (binned) spectral function. If using complex numbers and P.reim==true, we save triplets
// (energy,real part,imag part).
template<scalar S>
void SpectrumRealFreq<S>::savebins() {
  const auto fn = filename + ".bin";
  std::cout << " " << fn;
  std::ofstream Fbins = safe_open(fn, true); // true=binary!
  for (const auto &[e, w] : fspos.bins) {
    const auto [wr, wi] = reim(w);
    my_assert(e > 0.0);
    Fbins.write((char *)&e,  sizeof(double));
    Fbins.write((char *)&wr, sizeof(double));
    if (P.reim)
      Fbins.write((char *)&wi, sizeof(double));
  }
  for (const auto &[abse, w] : fsneg.bins) {
    const auto [wr, wi] = reim(w);
    const auto e = -abse;
    my_assert(e < 0.0); // attention!
    Fbins.write((char *)&e,  sizeof(double));
    Fbins.write((char *)&wr, sizeof(double));
    if (P.reim) 
      Fbins.write((char *)&wi, sizeof(double));
  }
}

// Energy mesh for spectral functions
inline std::vector<double> make_mesh(const Params &P) {
  const double broaden_min = P.get_broaden_min();
  std::vector<double> vecE; // Energies on the mesh
  for (double E = P.broaden_max; E > broaden_min; E /= P.broaden_ratio) vecE.push_back(E);
  return vecE;
}

template<scalar S>
void SpectrumRealFreq<S>::continuous() {
  using t_weight = weight_traits<S>;
  const double alpha  = P.alpha;
  const double omega0 = P.omega0 < 0.0 ? P.omega0_ratio * P.T : P.omega0;
  Spikes<S> densitypos, densityneg;
  const auto vecE = make_mesh(P); // Energies on the mesh
  for (const auto E : vecE) {
    t_weight valpos{}, valneg{};
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
  const auto fn = filename + ".dat";
  std::cout << " " << fn;
  std::ofstream Fdensity = safe_open(fn);
  densityneg.save(Fdensity, P.prec_xy, P.reim);
  densitypos.save(Fdensity, P.prec_xy, P.reim);
}

template<scalar S>
class GFMatsubara {
 private:
   const std::string name, algoname, filename;
   const Params &P;
   Matsubara<S> results;
 public:
   GFMatsubara(const std::string &name, const std::string &algoname, const std::string &filename, gf_type gt, const Params &P) : 
     name(name), algoname(algoname), filename(filename), P(P), results(P.mats, gt, P.T) {}
   void merge(const ChainMatsubara<S> &cm) {
     results.merge(cm.m);
   }
   void save() {
     fmt::print(fmt::emphasis::bold, "GF Matsubara: {} {} -> {}\n", name, algoname, filename);
     results.save(safe_open(filename + ".dat"), P.prec_xy);
   }
};

template<scalar S>
class TempDependence {
 private:
   const std::string name, algoname, filename;
   const Params &P;
   Spikes<S> results;
 public:
   TempDependence<S>(const std::string &name, const std::string &algoname, const std::string &filename, const Params &P) : 
     name(name), algoname(algoname), filename(filename),  P(P) {}
   void merge(const ChainTempDependence<S> &ctd) {
     std::copy(ctd.v.cbegin(), ctd.v.cend(), std::back_inserter(results));
   }
   void save() {
     fmt::print(fmt::emphasis::bold, "Temperature dependence: {} {} -> {}\n", name, algoname, filename);
     ranges::sort(results, sortfirst());
     results.save(safe_open(filename + ".dat"), P.prec_xy, P.reim);
   }
};

} // namespace

#endif
