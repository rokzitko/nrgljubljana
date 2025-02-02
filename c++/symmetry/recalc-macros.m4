namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020

define(`RECALC_TAB', {
  nrglog('f', "RECALC(fn=" << $1 << ", Iop=" << $2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include $1
      };
      auto cn = this->recalc_general(diag, cold, I1, Ip, recalc_table, $2);
      if (cn) cnew[II] = *cn;
    }
  }
})

define(`RECALC_F_TAB', {
  nrglog('f', "RECALC_F(fn=" << $1 << ", ch=" << $2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include $1
      };
      opch[$2][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
})

define(`RECALC_F_TAB_N', {
  nrglog('f', "RECALC_F(fn=" << $1 << ", ch=" << $2 << ", n=" << $3 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include $1
      };
      opch[$2][$3][II] = this->recalc_f(diag, I1, Ip, recalc_table);
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

}
