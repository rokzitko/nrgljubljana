// *** WARNING!!! Modify nrg-recalc-SPSU2C3.cc.m4, not nrg-recalc-SPSU2C3.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (S,P) subspaces

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  





#define xRECALC_F_TAB(a, b, c) 0;

template<typename SC>
Opch<SC> SymmetrySPSU2C3<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag, const QSrmax &qsrmax) {
  Opch<SC> opch = newopch<SC>(P);

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
 
  // CONVENTION: primed indeces are on the right side (ket)
  for(const auto &[Ip, eig]: diag) {
    Sspin ssp = Ip.get("SS");
    int p     = Ip.get("P");

    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(ssp + 1, (p + 0) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup0-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup0-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup0-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup0-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ssp - 1, (p + 0) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown0-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown0-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown0-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown0-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ssp + 1, (p + 1) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup1-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup1-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup1-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup1-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ssp - 1, (p + 1) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown1-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown1-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown1-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown1-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ssp + 1, (p + 2) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup2-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup2-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spinup2-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spinup2-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};

    I1 = Invar(ssp - 1, (p + 2) % 3);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-a.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown2-a.dat"
        };
        opch[0][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-b.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown2-b.dat"
        };
        opch[1][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2c3/spsu2c3-spindown2-c.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      [&]() ATTRIBUTE_NO_SANITIZE_DIV_BY_ZERO { 
        std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2c3/spsu2c3-spindown2-c.dat"
        };
        opch[2][0][II] = this->recalc_f(diag, qsrmax, I1, Ip, recalc_table);
      }();
    }
  }
};
#undef sqrt
#undef Power
  }
  
  }
  return opch;
}
