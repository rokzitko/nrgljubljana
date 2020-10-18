// spectral.h - Code for handling spectral function data
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _spectral_h_
#define _spectral_h_

// Container for holding spectral information represented by delta peaks. "Weight" is of type t_weight (complex).
template<typename S>
using t_delta_peak_tmpl = std::pair<double, typename traits<S>::t_weight>;
using t_delta_peak = t_delta_peak_tmpl<scalar>;

inline void outputxy(std::ostream &F, const double x, const std::complex<double> z, const bool imagpart, const double clip_tol_imag = 1e-10) {
  const auto [r, i] = reim(z);
  F << x << " " << r;
  if (imagpart) F << " " << (abs(i)>abs(r)*clip_tol_imag ? i : 0);
  F << std::endl;
}

template<typename S>
class Spikes_tmpl : public std::vector<t_delta_peak_tmpl<S>> {
 public:
   template<typename T>
     void save(T&& F, const int prec, const bool imagpart) {
       F << setprecision(prec);
       for (const auto &[e, w] : *this) outputxy(F, e, w, imagpart);
     }
   t_weight sum_weights() const {
     return ranges::accumulate(*this, t_weight{}, [](const auto &sum, const auto &p) { return sum+p.second; });
   }
};
using Spikes = Spikes_tmpl<scalar>;

#ifndef M_SQRTPI
#define M_SQRTPI 1.7724538509055160273
#endif

// cf. A. Weichselbaum and J. von Delft, cond-mat/0607497

// Modified log-Gaussian broadening kernel. For gamma=alpha/4, the
// kernel is symmetric in both arguments.
inline double BR_L(double e, double ept, double alpha, double omega0) {
  if ((e < 0.0 && ept > 0.0) || (e > 0.0 && ept < 0.0)) return 0.0;
  if (ept == 0.0) return 0.0;
  const double gamma = alpha/4.0;
  return exp(-sqr(log(e / ept) / alpha - gamma)) / (alpha * abs(e) * M_SQRTPI);
}

// Normalized to 1, width omega0. The kernel is symmetric in both
// arguments.
inline double BR_G(double e, double ept, double omega0) { return exp(-sqr((e - ept) / omega0)) / (omega0 * M_SQRTPI); }

// Note: 'ept' is the energy of the delta peak in the raw spectrum,
// 'e' is the energy of the data point in the broadened spectrum.
inline double BR_NEW(double e, double ept, double alpha, double omega0) {
  double part_l = BR_L(e, ept, alpha, omega0);
  // Most of the time we only need to compute part_l (loggaussian)
  if (abs(e) > omega0) return part_l;
  // Note: this is DIFFERENT from the broadening kernel proposed by
  // by Weichselbaum et al. This BR_h is a function of 'e', not
  // of 'ept'! This breaks the normalization at finite temperatures.
  // On the other hand, it gives nicer spectra when used in conjunction
  // with the self-energy trick.
  double BR_h = exp(-sqr(log(abs(e) / omega0) / alpha));
  my_assert(BR_h >= 0.0 && BR_h <= 1.0);
  return part_l * BR_h + BR_G(e, ept, omega0) * (1.0 - BR_h);
}

// Calculate "moment"-th spectral moment.
template<typename S>
CONSTFNC auto moment(const Spikes_tmpl<S> &s_neg, const Spikes_tmpl<S> &s_pos, const int moment) {
  using t_weight = typename traits<S>::t_weight;
  auto sumA = ranges::accumulate(s_pos, t_weight{}, [moment](auto s, const auto &x){ const auto &[e,w] = x; return s+w*pow(e,moment); });
  auto sumB = ranges::accumulate(s_neg, t_weight{}, [moment](auto s, const auto &x){ const auto &[e,w] = x; return s+w*pow(-e,moment); });
  return sumA+sumB;
}

CONSTFNC double fermi_fnc(const double omega, const double T) { 
  return 1 / (1 + exp(-omega / T)); 
}

CONSTFNC double bose_fnc(const double omega, const double T) {
  const auto d = 1.0 - exp(-omega / T);
  return d != 0.0 ? 1.0/d : std::numeric_limits<double>::quiet_NaN();
}

template<typename F, typename S>
auto sum(const Spikes_tmpl<S> &s, const bool invert, F && f) {
  using t_weight = typename traits<S>::t_weight;
  return ranges::accumulate(s, t_weight{}, [&f,invert](auto s, const auto &x){ const auto &[e,w] = x; return s+w*f(invert ? -e : e); });
}

// Integrated spectral function with a kernel as in FDT for fermions
template<typename S>
CONSTFNC auto fd_fermi(const Spikes_tmpl<S> &s_neg, const Spikes_tmpl<S> &s_pos, double const T) {
  auto fnc = [T](const auto x) { return fermi_fnc(x, T); };
  return sum(s_neg, true, fnc) + sum(s_pos, false, fnc);
}

// Ditto for bosons
template<typename S>
CONSTFNC auto fd_bose(const Spikes_tmpl<S> &s_neg, const Spikes_tmpl<S> &s_pos, double const T) {
  auto fnc = [T](const auto x) { return bose_fnc(x, T); };
  return sum(s_neg, true, fnc) + sum(s_pos, false, fnc);
}

#endif // _spectral_h_
