// debug.h - auxiliary debugging macros
// Copyright (C) 2005-2019 Rok Zitko

#ifndef _debug_h_
#define _debug_h_

// Non-recoverable error or explicit request to stop at some calculation stage
#define exit1(...) \
   { cout << endl << __VA_ARGS__ << endl; exit(1); }

#define debug(...) \
  cout << __VA_ARGS__ << endl;

// character c defines what to log
#define nrglog(c, ...) \
  if (logletter(c)) { cout << __VA_ARGS__ << endl; }

// Dump the value of variable STR to the standard output.
#define nrgdump(STR) \
   cout << #STR << "=" << (STR) << " "

#define nrgdump2(STR1, STR2) \
   cout << #STR1 << "=" << (STR1) << " " \
        << #STR2 << "=" << (STR2) << " "

#define nrgdump3(STR1, STR2, STR3) \
   cout << #STR1 << "=" << (STR1) << " " \
        << #STR2 << "=" << (STR2) << " " \
        << #STR3 << "=" << (STR3) << " "

#define nrgdump4(STR1, STR2, STR3, STR4) \
   cout << #STR1 << "=" << (STR1) << " " \
        << #STR2 << "=" << (STR2) << " " \
        << #STR3 << "=" << (STR3) << " " \
        << #STR4 << "=" << (STR4) << " "

#define nrgdump5(STR1, STR2, STR3, STR4, STR5) \
   cout << #STR1 << "=" << (STR1) << " " \
        << #STR2 << "=" << (STR2) << " " \
        << #STR3 << "=" << (STR3) << " " \
        << #STR4 << "=" << (STR4) << " " \
        << #STR5 << "=" << (STR5) << " "

#define nrgdump6(STR1, STR2, STR3, STR4, STR5, STR6) \
   cout << #STR1 << "=" << (STR1) << " " \
        << #STR2 << "=" << (STR2) << " " \
        << #STR3 << "=" << (STR3) << " " \
        << #STR4 << "=" << (STR4) << " " \
        << #STR5 << "=" << (STR5) << " " \
        << #STR6 << "=" << (STR6) << " "

#endif
