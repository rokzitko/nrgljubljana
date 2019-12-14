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

// Note: range limited to short numbers. P::piT is const double.
inline double ww(short n, matstype mt)
{
  switch (mt) {
    case matstype::bosonic: return P::Tpi * (2 * n);
    case matstype::fermionic: return P::Tpi * (2 * n + 1);
    default: my_assert_not_reached();
  }
}
inline double wb(short n) { return P::Tpi * (2 * n); }
inline double wf(short n) { return P::Tpi * (2 * n + 1); }

class Matsubara {
  private:
  typedef std::vector<pair<double, t_weight>> matsgf;
  matsgf v;
  matstype mt;
  public:
  Matsubara() = default;
  Matsubara(size_t mats, matstype _mt) : mt(_mt) {
    my_assert(mt == matstype::bosonic || mt == matstype::fermionic);
    for (size_t n = 0; n < mats; n++) v.push_back(make_pair(ww(n, mt), 0.0));
  }
  void add(size_t n, t_weight w) { v[n].second += w; }
  void save(ostream &F) const {
    F << setprecision(P::prec_xy);
    for (const auto &i : v) output(F, i.first, i.second, true);
  }
  t_weight total_weight() const {
    weight_bucket sum(v);
    return sum;
  }
  friend class SpectrumMatsubara;
};

#endif // _matsubara_h_
