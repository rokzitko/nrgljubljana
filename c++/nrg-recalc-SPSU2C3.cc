// *** WARNING!!! Modify nrg-recalc-SPSU2C3.cc.m4, not nrg-recalc-SPSU2C3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (S,P) subspaces

namespace SPSU2C3 {
#include "spsu2c3/spsu2c3-def.dat"
}

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers















#define xRECALC_F_TAB(a,b,c) 0;

// Driver routine for recalc_f()
void SymmetrySPSU2C3::recalc_irreduc(const DiagInfo &diag)
{
#ifdef NRG_COMPLEX
  // CONVENTION: primed indeces are on the right side (ket)
  LOOP_const(diag, isp) {
    Invar Ip = INVAR(isp);
    Sspin ssp = Ip.get("SS");
    int p = Ip.get("P");
    
    Invar I1;
    
// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(ssp+1, (p+0)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup0-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup0-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup0-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(ssp-1, (p+0)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown0-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown0-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown0-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};

    I1 = Invar(ssp+1, (p+1)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup1-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup1-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup1-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(ssp-1, (p+1)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown1-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown1-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown1-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};

    I1 = Invar(ssp+1, (p+2)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup2-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup2-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spinup2-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    
    I1 = Invar(ssp-1, (p+2)%3);
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-a.dat" << ", ch=" << 0 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown2-a.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[0][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-b.dat" << ", ch=" << 1 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown2-b.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[1][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
    {
 nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-c.dat" << ", ch=" << 2 << ", len=" << SPSU2C3::LENGTH_I_3CH << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include "spsu2c3/spsu2c3-spindown2-c.dat"
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == SPSU2C3::LENGTH_I_3CH);
        recalc_f(diag, a.opch[2][0], Ip, I1, recalc_table, SPSU2C3::LENGTH_I_3CH);
 }
};
#undef sqrt
#undef Power

}
#endif
}
