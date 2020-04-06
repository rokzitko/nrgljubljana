// Copyright (C) 2020 Rok Zitko

#ifndef _mp_h_
#define _mp_h_

#include <gmp.h>

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
  ~my_mpf() { mpf_clear(val); }
  inline operator mpf_t &() { return val; }
};

using vmpf = std::vector<my_mpf>;

#endif // _tridiag_h_
