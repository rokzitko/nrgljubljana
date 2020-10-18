// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

define(`RECALC_TAB', {
  nrglog('f', "RECALC(fn=" << $1 << ", Iop=" << $2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc recalc_table[] = {
#include $1
        };
        cnew[II] = recalc_general(diag, qsrmax, cold, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table), $2);
      }(); // immediately executed lambda
    }
  }
})

define(`RECALC_F_TAB', {
  nrglog('f', "RECALC_F(fn=" << $1 << ", ch=" << $2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include $1
        };
        opch[$2][0][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
})

define(`RECALC_F_TAB_N', {
  nrglog('f', "RECALC_F(fn=" << $1 << ", ch=" << $2 << ", n=" << $3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && recalc_f_coupled(I1, Ip, Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        Recalc_f recalc_table[] = {
#include $1
        };
        opch[$2][$3][II] = recalc_f(diag, qsrmax, I1, Ip, recalc_table, ARRAYLENGTH(recalc_table));
      }();
    }
  }
})

define(`ONETWO', 
  switch (P.channels) {
  case 1: { $1 } break;
  case 2: { $2 } break;
  default: my_assert_not_reached();
  })

define(`ONE23', 
  switch (P.channels) {
  case 1: { $1 } break;
  case 2: { $2 } break;
  case 3: { $3 } break;
  default: my_assert_not_reached();
  })
  
define(`ONE234', 
  switch (P.channels) {
  case 1: { $1 } break;
  case 2: { $2 } break;
  case 3: { $3 } break;
  case 4: { $4 } break;
  default: my_assert_not_reached();
})

define(`ONE2345', 
  switch (P.channels) {
  case 1: { $1 } break;
  case 2: { $2 } break;
  case 3: { $3 } break;
  case 4: { $4 } break;
  case 5: { $5 } break;
  default: my_assert_not_reached();
})
