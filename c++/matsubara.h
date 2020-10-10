#ifndef _matsubara_h_
#define _matsubara_h_

enum class matstype : unsigned char { bosonic, fermionic, bb, bf, fb }; // note the explicit type

string matstypestring(matstype mt) {
  switch (mt) {
    case matstype::bosonic: return "bosonic";
    case matstype::fermionic: return "fermionic";
    case matstype::bb: return "2*bosonic";
    case matstype::bf: return "bosonic+fermionic";
    case matstype::fb: return "fermionic+bosonic";
    default: my_assert_not_reached();
  }
}

// Note: range limited to short numbers.
inline double ww(short n, matstype mt, const double T)
{
  switch (mt) {
    case matstype::bosonic: return T * M_PI * (2 * n);
    case matstype::fermionic: return T * M_PI * (2 * n + 1);
    default: my_assert_not_reached();
  }
}
inline double wb(short n, const double T) { return T * M_PI * (2 * n); }
inline double wf(short n, const double T) { return T * M_PI * (2 * n + 1); }

class Matsubara {
 private:
   using matsgf = std::vector<pair<double, t_weight>>;
   matsgf v;
   matstype mt;
   double T;
 public:
   Matsubara() = delete;
   Matsubara(size_t mats, matstype mt, double T) : mt(mt), T(T) {
     my_assert(mt == matstype::bosonic || mt == matstype::fermionic);
     for (size_t n = 0; n < mats; n++) v.emplace_back(ww(n, mt, T), 0);
   }
   void add(size_t n, t_weight w) { v[n].second += w; }
   template <typename T> void save(T && F, int prec) const {
     F << setprecision(prec); // prec_xy
     for (const auto &[e, w] : v) output(F, e, w, true);
   }
   t_weight total_weight() const {
     weight_bucket sum(v);
     return sum;
   }
   friend class SpectrumMatsubara;
};

#endif // _matsubara_h_
