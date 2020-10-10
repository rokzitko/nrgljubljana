// spectral.h - Code for handling spectral function data
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _spectral_h_
#define _spectral_h_

// XXX redundancy: keep one of DELTA & t_delta_peak

// This structure is used in spec*.cc
struct DELTA {
  double energy{};
  t_weight weight; // XXX: {} ?
};

// Container for holding spectral information represented by delta
// peaks. "Weight" is of type t_weight (complex).
typedef pair<double, t_weight> t_delta_peak;
using Spikes = std::vector<t_delta_peak>;

// XXX: move to output.h
// used in matsubata.h and in save_densfunc()
inline void output(ostream &F, const double x, const t_weight y, const bool imagpart, const double clip_tol_imag = 1e-10) {
  const auto [r, i] = reim(y);
  F << x << " " << r;
  if (imagpart) F << " " << (abs(i) > abs(r) * clip_tol_imag ? i: 0.0);
  F << endl;
}

template<typename T>
void save_densfunc(T&& F, const Spikes &xy, const int prec, const bool imagpart) {
  F << setprecision(prec);
  for (const auto &[e, w] : xy) output(F, e, w, imagpart);
}

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

CONSTFNC t_weight sum_weights(const Spikes &s) {
  t_weight sum = 0.0;
  for (const auto &[e, w] : s) sum += w; // XXX stl algo
  return assert_isfinite(sum);
}

// Calculate "moment"-th spectral moment.
CONSTFNC t_weight moment(const Spikes &sneg, const Spikes &spos, int moment) {
  t_weight sumA = 0.0;
  for (const auto &[e, w] : spos) sumA += w * pow(e, moment);
  t_weight sumB = 0.0;
  for (const auto &[e, w] : sneg) sumB += w * pow(-e, moment);
  t_weight sum = sumA + sumB;
  return assert_isfinite(sum);
}

CONSTFNC double fermi_fnc(const double omega, const double T) { return 1 / (1 + exp(-omega / T)); }

CONSTFNC double bose_fnc(const double omega, const double T) { return 1 / (1 - exp(-omega / T)); }

// Integrated spectral function with a kernel as in FDT for fermions
CONSTFNC t_weight fd_fermi(const Spikes &sneg, const Spikes &spos, double const T) {
  t_weight sumA = 0.0;
  for (const auto &[e, w] : spos) sumA += w * fermi_fnc(e, T);
  t_weight sumB = 0.0;
  for (const auto &[e, w] : sneg) sumB += w * fermi_fnc(-e, T);
  t_weight sum = sumA + sumB;
  return assert_isfinite(sum);
}

// Ditto for bosons
CONSTFNC t_weight fd_bose(const Spikes &sneg, const Spikes &spos, double const T) {
  t_weight sumA = 0.0;
  for (const auto &[e, w] : spos) sumA += w * bose_fnc(e, T);
  t_weight sumB = 0.0;
  for (const auto &[e, w] : sneg) sumB += w * bose_fnc(-e, T);
  t_weight sum = sumA + sumB;
  return sum; // removed assertion due to false triggers
}

#endif // _spectral_h_
