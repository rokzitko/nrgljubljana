#ifndef _matsubara_hpp_
#define _matsubara_hpp_

#include <algorithm>
#include <array>
#include <cstddef>
#include <vector>
#include <utility>
#include <cmath>
#include <iomanip>
#include <ostream>

#include "traits.hpp"

#if NRG_ENABLE_APP_OPENMP
#include <omp.h>
#endif

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

namespace detail {

// Keep the transition order at each frequency. In parallel, stage one tile at
// a time so pair factors are evaluated once without storing the full block.
template<typename MakeTransition, typename Accumulate>
void accumulate_matsubara_transitions(const size_t frequencies, const size_t count,
                                      const MakeTransition &make_transition, const Accumulate &accumulate) {
  if (!frequencies || !count) return;
#if NRG_ENABLE_APP_OPENMP
  constexpr size_t tile_capacity = 256;
  if (count >= tile_capacity && frequencies > 1) {
    const auto threads = static_cast<int>(std::min(frequencies, static_cast<size_t>(omp_get_max_threads())));
    if (threads > 1) {
      std::array<decltype(make_transition(size_t{})), tile_capacity> transitions;
#pragma omp parallel num_threads(threads)
      {
        for (size_t first = 0; first < count; first += tile_capacity) {
          const auto length = std::min(tile_capacity, count - first);
#pragma omp for schedule(static)
          for (size_t offset = 0; offset < length; ++offset)
            transitions[offset] = make_transition(first + offset);
#pragma omp for schedule(static)
          for (size_t n = 0; n < frequencies; ++n)
            for (size_t offset = 0; offset < length; ++offset)
              accumulate(n, transitions[offset]);
        }
      }
      return;
    }
  }
#endif
  for (size_t pair = 0; pair < count; ++pair) {
    const auto transition = make_transition(pair);
    for (size_t n = 0; n < frequencies; ++n)
      accumulate(n, transition);
  }
}

} // namespace detail

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
   [[nodiscard]] auto frequency(const size_t n) const { return v[n].first; }
   void add(const size_t n, const t_weight &w) { v[n].second += w; }
   void merge(const Matsubara &m2) {
     my_assert(v.size() == m2.v.size());
     for (const auto n: range0(v.size())) {
       my_assert(v[n].first == m2.v[n].first);
       v[n].second += m2.v[n].second;
     }
   }
   template <typename T> void save(T && F, const int prec, const double clip_tol_imag) const {
     F << std::setprecision(prec); // prec_xy
     for (const auto &[e, w] : v) outputxy(F, e, w, true, clip_tol_imag);
   }
};

} // namespace NRG

#endif // _matsubara_hpp_
