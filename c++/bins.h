// bins.h - Binning of the spectral data
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _bins_h_
#define _bins_h_

// Binned spectral peaks. P::bins defines the number of bins per
// energy decade. The lowest and highest energies are defined by the
// zero-th and last NRG energy scale.

// bins[i].first is the energy of the representative point (it is NOT
// an interval boundary!). Thus bins[i].second represents the
// spectral weight at energy bins[i].first. This implies that when
// the spectral weight is collected [see function add()], the weight
// of a delta peak with energy between two consecutive representative
// points is divided between the two points proportionally with
// respect to the energies of the points. This also implies that
// during broadening, we may consider "bins" as delta peaks, rather
// than as interval representations (true bins). In other words, what
// we are doing here is coarse graining rather than binning!

class Bins
{
 private:
   size_t binsperdecade;
   size_t nrbins;
   double emin, emax;
   double log10emin, log10emax;  // base-10 log of the limits
   void setlimits();
   inline void add_std(double energy, t_weight weight);
   inline void add_acc(double energy, t_weight weight);
   void loggrid_std(); // standard NRG logarithmic mesh for spec. funcions
   void loggrid_acc(); // log grid with shifted accumulation point
   
 public:
   Spikes bins; // Note: Spikes is vector of (t_eigen,t_weight) pairs
   void loggrid(); 
   Bins() { loggrid(); } // default: logarithmic grid
   inline void add(double energy, t_weight weight);
   void merge(const Bins &b);
   void trim();
   t_weight total_weight() const { return sum_weights(bins); }
};

// These control the energy range of the bins. This is the shift in
// the (10-base) exponent of the top-most and bottom-most bins.
const double MAX_BIN_SHIFT = 2.0;
const double MIN_BIN_SHIFT = 2.0;

void Bins::setlimits()
{
   // NOTE: this will silently discard spectral peaks far outside the
   // conduction band!!
   emax = (P::emax > 0 ? P::emax : SCALE(0) * pow(10, MAX_BIN_SHIFT));
   emin = (P::emin > 0 ? P::emin : LAST_STEP_SCALE() / pow(10, MIN_BIN_SHIFT));
   // Trick: use ceil/floor to obtain uniform binning grids for
   // different values of the twist parameter z! 
   log10emin = floor(log10(emin));
   log10emax = ceil(log10(emax));
} 

void Bins::loggrid()
{
   my_assert(P::bins > 0);
   binsperdecade = P::bins;
   setlimits();
   if (P::accumulation > 0.0) 
      loggrid_acc();
   else 
      loggrid_std();
}

void Bins::loggrid_acc()
{
   const double a = P::accumulation;
   my_assert(a > 0.0);
   bins.resize(0);
   for (double e = emin; e <= emax; e *= pow(10.0, 1.0/binsperdecade)) {
      double x = (emax-a)/emax * e + a;
      bins.push_back(make_pair(x, 0.0));
   }
   if (P::linstep != 0.0) {
      my_assert(P::linstep > 0.0);
      for (double e = a; e > 0.0; e -= P::linstep)
	bins.push_back(make_pair(e, 0.0));
   }   
   bins.push_back(make_pair(DBL_MIN, 0.0)); // add zero point
   sort(begin(bins), end(bins), sortfirst());
   nrbins = bins.size();
   my_assert(nrbins >= 2);
}

void Bins::loggrid_std()
{
  nrbins = (size_t)((log10emax-log10emin)*binsperdecade+1.0);
  bins.resize(nrbins); // Note: Spikes is a vector type!
  for (size_t i = 0; i < nrbins; i++)
    bins[i] = make_pair(pow(10.0, log10emin + (double)i/binsperdecade), 0.0);
}

// Unbiased assignment of the spectral weight to bins.
inline void Bins::add(double energy, t_weight weight)
{
   if (abs(weight) < P::DISCARD_IMMEDIATELY * energy)
      return;
   if (P::accumulation > 0.0) 
      add_acc(energy, weight);
   else 
      add_std(energy, weight);
}

inline void Bins::add_std(double energy, t_weight weight)
{
   const double lg = log10(energy);
   const double pos = (lg-log10emin)*binsperdecade;
   const long index = long(pos); // must be signed!
   // Important: if 'energy' is lower than the lower limit of the
   // first interval, the weight is assigned to the first bin. This
   // is especially relevant for collecting the omega=0 data in
   // bosonic correlators. (rz, 25 Oct 2012)
   if (index < 0) {
      bins[0].second += weight;
      return;
   } else if (size_t(index) >= nrbins-1) {
      bins[nrbins-1].second += weight;
      return;
   }
   const double rem = pos-index;
   bins[index].second += (1.0-rem) * weight;
   bins[index+1].second += rem * weight;
}

inline void Bins::add_acc(double energy, t_weight weight)
{
   const size_t len = bins.size();
   for (size_t i = 0; i < len-1; i++) {
      const double e1 = bins[i].first;
      const double e2 = bins[i+1].first;
      my_assert (e1 < e2);
      if (e1 < energy && energy < e2) {
	 const double dx = e2-e1;
	 const double reldist = (energy-e1)/dx;
	 bins[i].second += (1.0-reldist) * weight;
	 bins[i+1].second += reldist * weight;
	 return;
      }
   }
   // Note: if no suitable interval is found, the weight is
   // discarded! XXX : should this be considered an error?
}

// Merge two bins. They need to agree in the representative energies
// (first element of the pairs).
void Bins::merge(const Bins &b)
{
   my_assert(bins.size() == b.bins.size());
   my_assert(bins.size() == nrbins);
   for (size_t i = 0; i < nrbins; i++) {
      my_assert(bins[i].first == b.bins[i].first);
      bins[i].second += b.bins[i].second;
   }
}

// Only keep bins which are "heavy" enough.
void Bins::trim()
{
   Spikes bins2;
   bucket discarded_weight_abs;
   size_t nr = bins.size();
   bins2.reserve(nr);
   // nr-1, because we need to compute the energy interval size 'ewidth'
   for (size_t i = 0; i < nr-1; i++) {
      const double e = bins[i].first;
      const t_weight w = bins[i].second;
      const double enext = bins[i+1].first; // increasing!
      my_assert(enext > e);
      const double ewidth = enext-e;
      if (abs(w) < P::DISCARD_TRIM * ewidth)
         discarded_weight_abs += abs(w);
      else 
         bins2.push_back(bins[i]);
   }
   // Always keep the last one.. This ensures that we keep
   // information about the energy range on which the calculations
   // has been performed.
   bins2.push_back(bins[nr-1]);
   bins.swap(bins2);
   if (discarded_weight_abs > 1e-8)
      cout << "WARNING: we are probably discarding too much weight!" << endl;
}

class Temp
{
public:
   Spikes v;
   void add_value(double energy, t_weight weight); // for ChainSpectrumTemp
   operator Spikes & () { return v; }
   typedef Spikes::iterator iterator;
   typedef Spikes::const_iterator const_iterator;
   iterator begin() { return v.begin(); }
   iterator end() { return v.end(); }
};
   
inline void Temp::add_value(double energy, t_weight weight)
{
   for (size_t i = 0; i < v.size(); i++) {
      if (v[i].first == energy) {
	 v[i].second += weight;
	 return;
      }
   }
   // or else 
   v.push_back(make_pair(energy, weight));
}

#endif // _bins_h_
