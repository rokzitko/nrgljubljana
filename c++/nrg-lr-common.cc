// This file is part of "NRG Ljubljana".
//   
// NRG Ljubljana is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the
// Free Software Foundation; either version 2 of the License, or (at your
// option) any later version.
//	
// NRG Ljubljana is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//			
// You should have received a copy of the GNU General Public License along
// with NRG Ljubljana; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA


// Routines for the parity quantum number
// Rok Zitko, rok.zitko@ijs.si, Aug 2006

#ifndef _nrg_lr_common_
#define _nrg_lr_common_

// Channel (left-right) symmetry support. We use class inheritance to
// solve this problem semi-generically for all symmetry types, as
// defined in nrg-invar-*.cc.

// We could just use an integer to store the parity information,
// however using a new enum type helps to ensure type-correctness (in
// C++ there is an implicit conversion from enum to int, but not vice
// versa).

// Even parity corresponds to eigenvalue 1, odd parity corresponds to
// eigenvalue -1.
enum eLRSYM { LRODD = -1, LREVEN = 1 };

inline eLRSYM int2lr(int ilr)
{
   my_assert(ilr == -1 || ilr == 1);
   return (ilr == -1 ? LRODD : LREVEN);
}

inline int lr2int(eLRSYM lr)
{
   return (lr == LRODD ? -1 : 1);
}

inline eLRSYM invert(eLRSYM lr)
{
   if (lr == LRODD)
     return LREVEN;
   else
     return LRODD;
}

// Parity "addition" law
inline eLRSYM lrcombine(const eLRSYM lr1, const eLRSYM lr2) 
{
   if (lr1 == lr2)
     return LREVEN;
   else
     return LRODD;
}

inline bool operator<(const eLRSYM lr1, const eLRSYM lr2)
{
  return lr2int(lr1) < lr2int(lr2);
}

#endif
