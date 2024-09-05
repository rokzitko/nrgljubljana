// bins.h - Binning of the spectral data
// Copyright (C) 2009-2020 Rok Zitko

#ifndef _bins_hpp_
#define _bins_hpp_

#include <cmath>
#include <range/v3/all.hpp>
#include "portabil.hpp"
#include "traits.hpp"
#include "spectral.hpp"
#include "params.hpp"
#include "numerics.hpp"

namespace NRG {

// Binned spectral peaks. P.bins defines the number of bins per energy decade. The lowest and highest energies are
// defined by the zero-th and last NRG energy scale.

// bins[i].first is the energy of the representative point (it is NOT an interval boundary!). Thus bins[i].second
// represents the spectral weight at energy bins[i].first. This implies that when the spectral weight is collected
// [see function add()], the weight of a delta peak with energy between two consecutive representative points is
// divided between the two points proportionally with respect to the energies of the points. This also implies that
// during broadening, we may consider "bins" as delta peaks, rather than as interval representations (true bins). In
// other words, what we are doing here is coarse graining rather than binning!

template<scalar S>
class Bins {
 private:
   using t_weight = weight_traits<S>;
   double emin{}, emax{};
   double log10emin{}, log10emax{}; // base-10 log of the limits
   void setlimits();
   inline void add_std(const double energy, const t_weight weight);
   inline void add_acc(const double energy, const t_weight weight);
   void loggrid_std(); // standard NRG logarithmic mesh for spec. funcions
   void loggrid_acc(); // log grid with shifted accumulation point
   void loggrid();
   const Params &P;

   // These control the energy range of the bins. This is the shift in the (10-base) exponent of the top-most and
   // bottom-most bins.
   inline static const double max_bin_shift = 2.0;
   inline static const double min_bin_shift = 2.0;
   inline static const double base = 10;
   inline static const double discarded_weight_warn_limit = 1e-8;

 public:
   Spikes<S> bins; // Note: Spikes is vector of (t_eigen,t_weight) pairs  // XXX: make private
   operator const Spikes<S> &() const { return bins; }
   operator Spikes<S> &() { return bins; }
   // auto & get() { return bins; }
   explicit Bins(const Params &P) : P(P) { loggrid(); } // default: logarithmic grid
   inline void add(const double energy, const t_weight weight);
   void merge(const Bins<S> &b);
   void trim();
   auto total_weight() const { return bins.sum_weights(); }
};

template<scalar S>
void Bins<S>::setlimits() {
  // NOTE: this will silently discard spectral peaks far outside the conduction band!!
  emax = (P.emax > 0 ? P.emax : P.SCALE(0) * pow(base, max_bin_shift));
  emin = (P.emin > 0 ? P.emin : P.last_step_scale() / pow(base, min_bin_shift));
  // Trick: use ceil/floor to obtain uniform binning grids for different values of the twist parameter z!
  log10emin = floor(log10(emin));
  log10emax = ceil(log10(emax));
}

template<scalar S>
void Bins<S>::loggrid() {
  my_assert(P.bins > 0);
  setlimits();
  if (P.accumulation > 0.0)
    loggrid_acc();
  else
    loggrid_std();
}

template<scalar S>
void Bins<S>::loggrid_acc() {
  const double a = P.accumulation;
  my_assert(a > 0.0);
  bins.resize(0);
  for (auto e = emin; e <= emax; e *= pow(base, 1.0 / P.bins))
    bins.emplace_back((emax - a) / emax * e + a, 0);
  if (P.linstep > 0)
    for (auto e = a; e > 0.0; e -= P.linstep) bins.emplace_back(e, 0);
  bins.emplace_back(DBL_MIN, 0); // add zero point
  ranges::sort(bins, sortfirst());
  my_assert(bins.size() >= 2);
}

template<scalar S>
void Bins<S>::loggrid_std() {
  const auto nrbins = (size_t)((log10emax - log10emin) * P.bins + 1.0);
  bins.resize(nrbins); // Note: Spikes is a vector type!
  for (const auto i : range0(nrbins)) bins[i] = { pow(base, log10emin + (double)i / P.bins), 0 };
}

// Unbiased assignment of the spectral weight to bins.
template<scalar S>
inline void Bins<S>::add(const double energy, const t_weight weight) {
  if (abs(weight) < P.discard_immediately * energy) return;
  if (P.accumulation > 0.0)
    add_acc(energy, weight);
  else
    add_std(energy, weight);
}

template<scalar S>
inline void Bins<S>::add_std(const double energy, const t_weight weight) {
  // Important: if 'energy' is lower than the lower limit of the first interval, the weight is assigned to the first
  // bin. This is especially relevant for collecting the omega=0 data in bosonic correlators. (rz, 25 Oct 2012)
  if (energy < emin) { // handle this special case separately (for reasons of efficiency)
    bins[0].second += weight;
    return;
  }
  const double log10e  = log10(energy);
  const double x = (log10e - log10emin) * P.bins;
  const double int_part = floor(x);
  if (int_part < 0) {
    bins.front().second += weight;
  } else if (const auto index = size_t(int_part) ; index >= bins.size() - 1) {
    bins.back().second += weight;
  } else {
    const double rem = x - int_part;
    bins[index].second += (1.0 - rem) * weight;
    bins[index + 1].second += rem * weight;
  }
}

template<scalar S>
inline void Bins<S>::add_acc(const double energy, const t_weight weight) {
  for (const auto i: range0(bins.size()-1)) {
    auto &[e1, w1] = bins[i]; // non-const
    auto &[e2, w2] = bins[i+1]; // non-const
    my_assert(e1 < e2);
    if (e1 < energy && energy < e2) {
      const auto dx      = e2 - e1;
      const auto reldist = (energy - e1) / dx;
      w1 += (1.0 - reldist) * weight;
      w2 += reldist * weight;
      return;
    }
  }
  // Note: if no suitable interval is found, the weight is discarded!
}

// Merge two bins. They need to agree in the representative energies
// (first element of the pairs).
template<scalar S>
void Bins<S>::merge(const Bins<S> &b) {
  my_assert(bins.size() == b.bins.size());
  for (const auto i: range0(bins.size())) {
    auto &[e1, w1] = bins[i];
    const auto &[e2, w2] = b.bins[i];
    my_assert(e1 == e2);
    w1 += w2;
  }
}

// Only keep bins which are "heavy" enough.
template<scalar S>
void Bins<S>::trim() {
  Spikes<S> orig{};
  orig.swap(bins);
  bucket discarded_weight_abs;
  // nr-1, because we need to compute the energy interval size 'ewidth'
  for (const auto i: range0(orig.size()-1)) {
    const auto [e, w] = orig[i];
    const auto e_next = orig[i+1].first;
    my_assert(e_next > e);  // increasing!
    const auto e_width = e_next - e;
    if (abs(w) >= P.discard_trim * e_width)
      bins.push_back(orig[i]);
    else
      discarded_weight_abs += abs(w);
  }
  // Always keep the last one.. This ensures that we keep information about the energy range on which the
  // calculations has been performed.
  bins.push_back(orig.back());
  if (discarded_weight_abs > discarded_weight_warn_limit) std::cout << "WARNING: we are probably discarding too much weight!" << std::endl;
}

template<scalar S, typename t_weight = weight_traits<S>>
class Temp : public Spikes<S> {
 private:
   const Params &P;
 public:
   explicit Temp(const Params &P) : P(P) {}
   void add_value(const double energy, const t_weight &weight) {
     for (auto & [e, w] : *this) {
       if (e == energy) {
         w += weight;
         return;
       }
     }
     // or else
     this->emplace_back(energy, weight);
   }
};

} // namespace

#endif // _bins_hpp_
