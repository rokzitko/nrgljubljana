// *** WARNING!!! Modify nrg-recalc-SL.cc.m4, not nrg-recalc-SL.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, June 2006, Nov 2007
// This file pertains to the spinless-fermions code.

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















namespace SL {
#include "sl/sl-1ch-def.dat"
#include "sl/sl-2ch-def.dat"
#include "sl/sl-3ch-def.dat"
}

// Recalculate matrix elements of a doublet tensor operator
void SymmetrySL::recalc_doublet(DiagInfo &diag,
                    MatrixElements &cold,
                    MatrixElements &cnew)
{
  LOOP(diag, is1) {
    Invar I1 = INVAR(is1);
    Number q1 = I1.get("Q");
    Invar Ip = Invar(q1-1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC(fn=" << "sl/sl-1ch-doublet.dat" << ", len=" << SL::LENGTH_D_1CH << ", Iop=" << Invar(1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "sl/sl-1ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_D_1CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SL::LENGTH_D_1CH, Invar(1));
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC(fn=" << "sl/sl-2ch-doublet.dat" << ", len=" << SL::LENGTH_D_2CH << ", Iop=" << Invar(1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "sl/sl-2ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_D_2CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SL::LENGTH_D_2CH, Invar(1));
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC(fn=" << "sl/sl-3ch-doublet.dat" << ", len=" << SL::LENGTH_D_3CH << ", Iop=" << Invar(1) << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include "sl/sl-3ch-doublet.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_D_3CH);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, SL::LENGTH_D_3CH, Invar(1));
 }
}
}
break;
default:
my_assert_not_reached();
};
  }
}

// Driver routine for recalc_f()
void SymmetrySL::recalc_irreduc(const DiagInfo &diag)
{
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Number qp = Ip.get("Q");
    Invar I1 = Invar(qp+1);
    switch (channels) {
case 1: {
{
 nrglog('f', "RECALC_F(fn=" << "sl/sl-1ch-a.dat" << ", ch=" << 0 << ", len=" << SL::LENGTH_I_1CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-1ch-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_1CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SL::LENGTH_I_1CH);
 }
}
}
break;
case 2: {
{
 nrglog('f', "RECALC_F(fn=" << "sl/sl-2ch-a.dat" << ", ch=" << 0 << ", len=" << SL::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-2ch-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_2CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SL::LENGTH_I_2CH);
 }
};
	   {
 nrglog('f', "RECALC_F(fn=" << "sl/sl-2ch-b.dat" << ", ch=" << 1 << ", len=" << SL::LENGTH_I_2CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-2ch-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_2CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SL::LENGTH_I_2CH);
 }
}
}
break;
case 3: {
{
 nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-a.dat" << ", ch=" << 0 << ", len=" << SL::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-3ch-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SL::LENGTH_I_3CH);
 }
};
    	   {
 nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-b.dat" << ", ch=" << 1 << ", len=" << SL::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-3ch-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SL::LENGTH_I_3CH);
 }
};
	   {
 nrglog('f', "RECALC_F(fn=" << "sl/sl-3ch-c.dat" << ", ch=" << 2 << ", len=" << SL::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "sl/sl-3ch-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SL::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SL::LENGTH_I_3CH);
 }
}
	   
}
break;
default:
my_assert_not_reached();
};
  }
}

#undef QDIFF
#define QDIFF(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef QTOT
#define QTOT(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N1
#define N1(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N2
#define N2(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

#undef N3
#define N3(i1, ip, ch, value) recalc1_global(diag, I1, cn, i1, ip, value)

void SymmetrySL::recalc_global(DiagInfo &diag,
                               string name,
                               MatrixElements &cnew)
{
   if (name == "Qdiff") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 2:
#include "sl/sl-2ch-qdiff.dat"
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
        case 2:
#include "sl/sl-2ch-qtot.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

   if (name == "N1") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N1.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }

   if (name == "N2") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N2.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
   
   if (name == "N3") {
     LOOP(diag, is1) {
       Invar I1 = INVAR(is1);
       const Twoinvar II = make_pair(I1, I1);
       Matrix & cn = cnew[II];
       switch (channels) {
        case 3:
#include "sl/sl-3ch-N3.dat"
         break;
        default:
         my_assert_not_reached();
       }
     } // LOOP
   }
}
