#include <complex>
#include <vector>
#include <type_traits>

namespace NRG {

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Oct 2015
// This file pertains to (S,P) subspaces

#define xRECALC_F_TAB(a, b, c) 0;

template<typename SC>
Opch<SC> SymmetrySPSU2C3<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  Opch<SC> opch(P);

  if constexpr (std::is_same_v<SC, std::complex<double>>) {
 
  // CONVENTION: primed indeces are on the right side (ket)
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    int p     = Ip.get("P");

    Invar I1;

// TRICK: ensure we are evaluating the expressions in the complex plane
#undef Power
#define Power(x, y) pow(cmpl(x), cmpl(y))

#undef sqrt
#define sqrt(x) csqrt(x)

    I1 = Invar(ssp + 1, (p + 0) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup0-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup0-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup0-c.dat"
    }; });

    I1 = Invar(ssp - 1, (p + 0) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown0-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown0-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown0-c.dat"
    }; });

    I1 = Invar(ssp + 1, (p + 1) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup1-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup1-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup1-c.dat"
    }; });

    I1 = Invar(ssp - 1, (p + 1) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown1-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown1-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown1-c.dat"
    }; });

    I1 = Invar(ssp + 1, (p + 2) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup2-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup2-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spinup2-c.dat"
    }; });

    I1 = Invar(ssp - 1, (p + 2) % 3);
    this->recalc_f_table(diag, I1, Ip, opch[0][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown2-a.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[1][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown2-b.dat"
    }; });
    this->recalc_f_table(diag, I1, Ip, opch[2][0], [&] { return RecalcFTable<SC>{
#include "spsu2c3/spsu2c3-spindown2-c.dat"
    }; });
#undef sqrt
#undef Power
  }
  
  }
  return opch;
}

}
