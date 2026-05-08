#ifndef _matsubara_hpp_
#define _matsubara_hpp_

#include <vector>
#include <utility>
#include <cmath>
#include <iomanip>
#include "traits.hpp"

namespace NRG {

enum class gf_type { bosonic, fermionic }; // req'd in matsubara.h

// Sign factor in GFs for bosonic/fermionic operators
inline constexpr auto S_BOSONIC   = +1;
inline constexpr auto S_FERMIONIC = -1;
inline int gf_sign(const gf_type gt) { return gt == gf_type::bosonic ? S_BOSONIC : S_FERMIONIC; }

[[nodiscard]] inline auto ww(const size_t n, const gf_type mt, const double T)
{
  const auto nd = static_cast<double>(n);
  switch (mt) {
    case gf_type::bosonic: return T * M_PI * (2.0 * nd);
    case gf_type::fermionic: return T * M_PI * (2.0 * nd + 1.0);
    default: my_assert_not_reached();
  }
}
[[nodiscard]] inline auto wb(const size_t n, const double T) { return T * M_PI * (2.0 * static_cast<double>(n)); }
[[nodiscard]] inline auto wf(const size_t n, const double T) { return T * M_PI * (2.0 * static_cast<double>(n) + 1.0); }

template<scalar S, typename t_weight = weight_traits<S>>
class Matsubara {
 private:
   using t_temp = typename traits<S>::t_temp;
   using matsgf = std::vector<std::pair<t_temp, t_weight>>;
   matsgf v;
   gf_type mt;
   t_temp T;
 public:
   Matsubara() = delete;
   Matsubara(const size_t mats, const gf_type mt_, const t_temp T_) : mt(mt_), T(T_) {
     my_assert(mt_ == gf_type::bosonic || mt_ == gf_type::fermionic);
     for (const auto n: range0(mats)) v.emplace_back(ww(n, mt_, T_), 0);
   }
   void add(const size_t n, const t_weight &w) { v[n].second += w; }
   void merge(const Matsubara &m2) {
     my_assert(v.size() == m2.v.size());
     for (const auto n: range0(v.size())) {
       my_assert(v[n].first == m2.v[n].first);
       v[n].second += m2.v[n].second;
     }
   }
   template <typename T> void save(T && F, const int prec) const {
     F << std::setprecision(prec); // prec_xy
     for (const auto &[e, w] : v) outputxy(F, e, w, true);
   }
};

} // namespace NRG

#endif // _matsubara_hpp_
