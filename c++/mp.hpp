// Copyright (C) 2020 Rok Zitko

#ifndef _mp_hpp_
#define _mp_hpp_

#include <gmp.h>
#include <vector>

namespace NRG {

// Wrapper class for arbitrary precision numbers
class my_mpf {
  private:
  mpf_t val{};

  public:
  my_mpf() { mpf_init(val); }
  // Copy constructor is mandatory!
  my_mpf(const my_mpf &x) {
    mpf_init(val);
    mpf_set(val, x.val);
  }
  my_mpf(my_mpf && x) = delete;
  my_mpf & operator=(const my_mpf & x) = delete;
  my_mpf & operator=(my_mpf && x) = delete;
  ~my_mpf() { mpf_clear(val); }
  inline operator mpf_t &() { return val; }
};

using vmpf = std::vector<my_mpf>;

} // namespace NRG

#endif
