namespace NRG {

// *** WARNING!!! Modify nrg-recalc-SPSU2.cc.m4, not nrg-recalc-SPSU2.cc !!!

// Quantum number dependant recalculation routines
// Rok Zitko, rok.zitko@ijs.si, Feb 2006, Dec 2007
// This file pertains to (S) subspaces

namespace NRG {

// m4 macros for nrg-recalc-*.cc files
// Rok Zitko, rok.zitko@ijs.si, 2007-2020










  




}


template<typename SC>
MatrixElements<SC> SymmetrySPSU2<SC>::recalc_doublet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-doubletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-doubletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-doubletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-doubletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-doubletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-doubletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1 + 1);
      {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletp.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-doubletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
};

      Ip = Invar(ss1 - 1);
      {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-doubletm.dat" << ", Iop=" << Invar(2) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-doubletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(2));
    }
  }
};
    }
  }
  return cnew;
}

template<typename SC>
Opch<SC> SymmetrySPSU2<SC>::recalc_irreduc(const Step &step, const DiagInfo<SC> &diag) const {
  my_assert(!P.substeps);
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-spinupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-spinupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
	         {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-spinupb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupa.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spinupa.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
       	   {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spinupb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
	         {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spinupc.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spinupc.dat"
      };
      opch[2][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    I1 = Invar(ssp-1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-spindowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-spindowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-2ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-spindownb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindowna.dat" << ", ch=" << 0 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spindowna.dat"
      };
      opch[0][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindownb.dat" << ", ch=" << 1 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spindownb.dat"
      };
      opch[1][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
           {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-3ch-spindownc.dat" << ", ch=" << 2 << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-spindownc.dat"
      };
      opch[2][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    // Note: for 3ch cases, the lengths in the three channels are not the same!
    // The same thing occurs for all SU(2)_spin cases, for instance for symtype=QS.
    // RZ, oct 2015
  }
  return opch;
}

template<typename SC>
OpchChannel<SC> SymmetrySPSU2<SC>::recalc_irreduc_substeps(const Step &step, const DiagInfo<SC> &diag, const int M) const {
  my_assert(P.substeps);
  Opch<SC> opch(P);
  for(const auto &[Ip, eig]: diag) {
    int ssp = Ip.get("SS");
    Invar I1;

    I1 = Invar(ssp + 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spinupa.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-spinupa.dat"
      };
      opch[M][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};

    I1 = Invar(ssp - 1);
    {
  nrglog('f', "RECALC_F(fn=" << "spsu2/spsu2-1ch-spindowna.dat" << ", ch=" << M << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip) && this->recalc_f_coupled(I1, Ip, this->Invar_f)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc_f<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-spindowna.dat"
      };
      opch[M][0][II] = this->recalc_f(diag, I1, Ip, recalc_table);
    }
  }
};
  }
  return opch[M];
}

template<typename SC>
MatrixElements<SC> SymmetrySPSU2<SC>::recalc_triplet(const DiagInfo<SC> &diag, const MatrixElements<SC> &cold) const {
  MatrixElements<SC> cnew;
  if (!P.substeps) {
    for(const auto &[I1, eig]: diag) {
      int ss1 = I1.get("SS");
      Invar Ip;

      Ip = Invar(ss1);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-triplets.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-triplets.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-triplets.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-triplets.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-triplets.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-triplets.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1+2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-tripletp.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-tripletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-tripletp.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-tripletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  case 3: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-3ch-tripletp.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-3ch-tripletp.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  default: my_assert_not_reached();
  };

    Ip = Invar(ss1-2);
    switch (P.channels) {
  case 1: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-1ch-tripletm.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-1ch-tripletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  case 2: { {
  nrglog('f', "RECALC(fn=" << "spsu2/spsu2-2ch-tripletm.dat" << ", Iop=" << Invar(3) << ")");
  auto II = Twoinvar(I1, Ip);
  if (diag.count(I1) && diag.count(Ip)) {
    if (diag.at(I1).getnrstored() && diag.at(Ip).getnrstored()) {
      std::initializer_list<Recalc<SC>> recalc_table = {
#include "spsu2/spsu2-2ch-tripletm.dat"
      };
      cnew[II] = this->recalc_general(diag, cold, I1, Ip, recalc_table, Invar(3));
    }
  }
} } break;
  default: my_assert_not_reached();
  };
    }
  } else my_assert_not_reached();
  return cnew;
}

#undef CHARGE
#define CHARGE(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef QDIFF
#define QDIFF(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q1
#define Q1(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef Q2
#define Q2(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

#undef ISOSPINZ
#define ISOSPINZ(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value)

// NOTE: the transverse components of the isospin depend on the site
// index! This is taken into account by appropriately multiplying 'value'
// by (-1)^N.

#undef ISOSPINX
#define ISOSPINX(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINP
#define ISOSPINP(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

#undef ISOSPINM
#define ISOSPINM(i1, ip, ch, value) this->recalc1_global(diag, I1, cn, i1, ip, value *psgn(step.getnn() + 1))

template<typename SC>
void SymmetrySPSU2<SC>::recalc_global(const Step &step, const DiagInfo<SC> &diag, const std::string name, MatrixElements<SC> &cnew) const {
  // NOTE: none of these are implemented for substeps==true.

  if (name == "Qtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Qtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Qtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Qdiff") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-qdiff.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q1") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q1.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Q2") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 2:
#include "spsu2/spsu2-2ch-q2.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iztot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iztot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iztot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Ixtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Ixtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Ixtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Iptot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Iptot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Iptot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  if (name == "Imtot") {
    for(const auto &[I1, eig]: diag) {
      const Twoinvar II = {I1, I1};
      Matrix &cn        = cnew[II];
      switch (P.channels) {
        case 1:
#include "spsu2/spsu2-1ch-Imtot.dat"
          break;
        case 2:
#include "spsu2/spsu2-2ch-Imtot.dat"
          break;
        default: my_assert_not_reached();
      }
    }
    return;
  }

  my_assert_not_reached();
}

}
