// Copyright (C) 2020 Rok Zitko

#ifndef _mp_hpp_
#define _mp_hpp_

#include <cstddef>
#include <gmp.h>
#include <vector>

namespace NRG {

inline constexpr mp_bitcnt_t FDM_MPF_PRECISION = 400; // bits, not decimal digits

// Wrapper class for arbitrary precision numbers
class my_mpf {
  private:
  mpf_t val{};

   public:
   my_mpf() { mpf_init(val); }
   explicit my_mpf(const mp_bitcnt_t precision) { mpf_init2(val, precision); }
   // Copy constructor is mandatory!
   my_mpf(const my_mpf &x) {
     mpf_init2(val, mpf_get_prec(x.val));
     mpf_set(val, x.val);
   }
  my_mpf(my_mpf && x) = delete;
  my_mpf & operator=(const my_mpf & x) = delete;
  my_mpf & operator=(my_mpf && x) = delete;
  ~my_mpf() { mpf_clear(val); }
  inline operator mpf_t &() { return val; }
};

using vmpf = std::vector<my_mpf>;

inline vmpf make_vmpf(const std::size_t size, const mp_bitcnt_t precision) {
  return vmpf(size, my_mpf(precision));
}

} // namespace NRG

#endif
