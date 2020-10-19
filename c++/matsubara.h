#ifndef _matsubara_h_
#define _matsubara_h_

// Note: range limited to short numbers.
inline double ww(short n, gf_type mt, const double T)
{
  switch (mt) {
    case gf_type::bosonic: return T * M_PI * (2 * n);
    case gf_type::fermionic: return T * M_PI * (2 * n + 1);
    default: my_assert_not_reached();
  }
}
inline double wb(short n, const double T) { return T * M_PI * (2 * n); }
inline double wf(short n, const double T) { return T * M_PI * (2 * n + 1); }

class Matsubara {
 private:
   using matsgf = std::vector<pair<double, cmpl>>;
   matsgf v;
   gf_type mt;
   double T;
 public:
   Matsubara() = delete;
   Matsubara(size_t mats, gf_type mt, double T) : mt(mt), T(T) {
     my_assert(mt == gf_type::bosonic || mt == gf_type::fermionic);
     for (const auto n: range0(mats)) v.emplace_back(ww(n, mt, T), 0);
   }
   void add(size_t n, cmpl w) { v[n].second += w; }
   void merge(const Matsubara &m2) {
     my_assert(v.size() == m2.v.size());
     for (const auto n: range0(v.size())) {
       my_assert(v[n].first == m2.v[n].first);
       v[n].second += m2.v[n].second;
     }
   }
   template <typename T> void save(T && F, int prec) const {
     F << setprecision(prec); // prec_xy
     for (const auto &[e, w] : v) outputxy(F, e, w, true);
   }
   cmpl total_weight() const {
     weight_bucket sum(v); // XXX sum function?
     return sum;
   }
};

#endif // _matsubara_h_
