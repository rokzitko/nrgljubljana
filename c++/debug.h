// debug.h - auxiliary debugging macros
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _debug_h_
#define _debug_h_

// Non-recoverable error or explicit request to stop at some calculation stage
#define exit1(...)                                                                                                                                   \
  {                                                                                                                                                  \
    std::cout << std::endl << __VA_ARGS__ << std::endl;                                                                                              \
    exit(1);                                                                                                                                         \
  }

#define debug(...) std::cout << __VA_ARGS__ << std::endl;

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

#endif
