// *** WARNING!!! Modify nrg-recalc-NONE.cc.m4, not nrg-recalc-NONE.cc !!!

// Quantum number dependent recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, April 2010
// This file pertains to the case with no symmetry.

namespace NONE {
#include "none/none-1ch-def.dat"
#include "none/none-2ch-def.dat"
}

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers

















// Driver routine for recalc_f()
void SymmetryNONE::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Invar I1 = Invar();

    switch (channels) {
case 1: {
{
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-1ch-a-CR-DO.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_1CH);
    recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, NONE::LENGTH_I_1CH);
 }
};
            {
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-1ch-a-CR-UP.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_1CH);
    recalc_f(diag, a.opch[0][1], Ip, I1, recalc_table, NONE::LENGTH_I_1CH);
 }
};
}
break;
case 2: {
{
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-2ch-a-CR-DO.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_2CH);
    recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, NONE::LENGTH_I_2CH);
 }
};
	    {
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-2ch-b-CR-DO.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_2CH);
    recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, NONE::LENGTH_I_2CH);
 }
};
            {
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-2ch-a-CR-UP.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_2CH);
    recalc_f(diag, a.opch[0][1], Ip, I1, recalc_table, NONE::LENGTH_I_2CH);
 }
};
            {
 if (diag.count(I1)) {
    struct Recalc_f recalc_table[]={
       #include "none/none-2ch-b-CR-UP.dat"
    };
    BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_I_2CH);
    recalc_f(diag, a.opch[1][1], Ip, I1, recalc_table, NONE::LENGTH_I_2CH);
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetryNONE::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Invar Ip = Invar();

    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "none/none-1ch-doublet.dat" << ", len=" << NONE::LENGTH_D_1CH << ", Iop=" << Invar() << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "none/none-1ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, NONE::LENGTH_D_1CH, Invar());
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "none/none-2ch-doublet.dat" << ", len=" << NONE::LENGTH_D_2CH << ", Iop=" << Invar() << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "none/none-2ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == NONE::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, NONE::LENGTH_D_2CH, Invar());
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}

#undef SPINX
#define SPINX(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
#undef SPINZ
#define SPINZ(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

// Isospin operator need an appropriate phase factor (bipartite
// sublattice index)
#define USEISOFACTOR

#if defined(USEISOFACTOR)
 #define ISOFACTOR psgn(getnn()+1)
#else
 #define ISOFACTOR 1
#endif

#ifdef NRG_COMPLEX
 #undef SPINY
 #define SPINY(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)
 
 #undef ISOSPINY
 #define ISOSPINY(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * complex<double>(ISOFACTOR))

 #undef Complex
 #define Complex(x,y) cmpl(x,y)
#endif // NRG_COMPLEX

#undef CHARGE
#define CHARGE(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * ISOFACTOR)

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * ISOFACTOR)

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) \
  recalc1_global(diag, I1, cn, i1, ip, value * ISOFACTOR)

void SymmetryNONE::recalc_global(DiagInfo &diag,
                                 string name,
                                 MatrixElements &cnew)
{
   if (name == "SZtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "none/none-1ch-spinz.dat"
         break;
        case 2:
#include "none/none-2ch-spinz.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

#ifdef NRG_COMPLEX
   if (name == "SYtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "none/none-1ch-spiny.dat"
         break;
        case 2:
#include "none/none-2ch-spiny.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
#endif

  if (name == "SXtot") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 1:
#include "none/none-1ch-spinx.dat"
         break;
        case 2:
#include "none/none-2ch-spinx.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

 if (name == "Qtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Qtot.dat"
        break;
      case 2:
#include "none/none-2ch-Qtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

  if (name == "Iztot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Iztot.dat"
        break;
      case 2:
#include "none/none-2ch-Iztot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

 if (name == "Ixtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Ixtot.dat"
        break;
      case 2:
#include "none/none-2ch-Ixtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

#ifdef NRG_COMPLEX
 if (name == "Iytot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Iytot.dat"
        break;
      case 2:
#include "none/none-2ch-Iytot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }
#endif

if (name == "Iptot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Iptot.dat"
        break;
      case 2:
#include "none/none-2ch-Iptot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

 if (name == "Imtot") {
    LOOP(diag, is1) {
      Invar I1 = INVAR(is1);
      const Twoinvar II = make_pair(I1, I1);
      Matrix & cn = cnew[II];
      switch (channels) {
      case 1:
#include "none/none-1ch-Imtot.dat"
        break;
      case 2:
#include "none/none-2ch-Imtot.dat"
        break;
      default:
        my_assert_not_reached();
      }
    } // LOOP
  }

}
