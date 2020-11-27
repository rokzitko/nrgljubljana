// debug.h - auxiliary debugging macros
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _debug_hpp_
#define _debug_hpp_

#include <iostream>
#include <cstdlib> // exit

namespace NRG {

// character c defines what to log
#define nrglog(c, ...)   if (P.logletter(c)) { std::cout << __VA_ARGS__ << std::endl; }
#define nrglogdp(c, ...) if (DP.logletter(c)) { std::cout << __VA_ARGS__ << std::endl; }

//#define MPI_DEBUG

#ifdef MPI_DEBUG
 #define mpilog(...) { std::cout << __VA_ARGS__ << " [rank " << myrank() << "]" << std::endl; }
#else
 #define mpilog(...) {}
#endif

#define nrgdump3(STR1, STR2, STR3) #STR1 << "=" << (STR1) << " " << #STR2 << "=" << (STR2) << " " << #STR3 << "=" << (STR3) << " "

#define nrgdump4(STR1, STR2, STR3, STR4)                                                                                                     \
  #STR1 << "=" << (STR1) << " " << #STR2 << "=" << (STR2) << " " << #STR3 << "=" << (STR3) << " " << #STR4 << "=" << (STR4) << " "

#define nrgdump5(STR1, STR2, STR3, STR4, STR5)                                                                                               \
  #STR1 << "=" << (STR1) << " " << #STR2 << "=" << (STR2) << " " << #STR3 << "=" << (STR3) << " " << #STR4 << "=" << (STR4) << " " << #STR5  \
       << "=" << (STR5) << " "

} // namespace

#endif
