// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2015

// m4 comment: $2 is length, $3,... are quantum numbers
define(`RECALC_TAB',
{
 nrglog('f', "RECALC(fn=" << $1 << ", len=" << $2 << ", Iop=" << $3 << ")");
 if (diag.count(Ip)) {
  struct Recalc recalc_table[]={
#include $1
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == $2);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, $2, $3);
 }
})

define(`RECALC_TAB_NUM',
{
 nrglog('f', "RECALC(fn=" << $1 << ", len=" << $2 << ", Iop=" << $3 << ")");
 if (diag.count(Ip)) {
  struct Recalc_num recalc_num_table[]={
#include $1
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_num_table) == $2);
	struct Recalc recalc_table[$2];
	lookup($2, recalc_num_table, recalc_table, I1);
        recalc_general(diag, cold, cnew, I1, Ip, 
                       recalc_table, $2, $3);
 }
})

define(`RECALC_F_TAB',
{
 nrglog('f', "RECALC_F(fn=" << $1 << ", ch=" << $2 << ", len=" << $3 << ")");
 if (diag.count(I1)) {
  struct Recalc_f recalc_table[]={
#include $1
        };
	BOOST_STATIC_ASSERT(ARRAYLENGTH(recalc_table) == $3);
        recalc_f(diag, a.opch[$2][0], Ip, I1, recalc_table, $3);
 }
})

define(`ONETWO',
switch (channels) {
case 1: {
$1
}
break;
case 2: {
$2
}
break;
default:
my_assert_not_reached();
})

define(`ONE23',
switch (channels) {
case 1: {
$1
}
break;
case 2: {
$2
}
break;
case 3: {
$3
}
break;
default:
my_assert_not_reached();
})

define(`ONE234',
switch (channels) {
case 1: {
$1
}
break;
case 2: {
$2
}
break;
case 3: {
$3
}
break;
case 4: {
$4
}
break;
default:
my_assert_not_reached();
})

define(`ONE2345',
switch (channels) {
case 1: {
$1
}
break;
case 2: {
$2
}
break;
case 3: {
$3
}
break;
case 4: {
$4
}
break;
case 5: {
$5
}
break;
default:
my_assert_not_reached();
})
