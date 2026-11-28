#include <string>
#include <vector>

#include "sym-QSZ-impl.hpp"

namespace NRG {

// U(1) charge and Z_4 spin symmetry: Q and 2S_z modulo 4 conserved, for normal-state problems with magnetic
// anisotropy D S_z^2 + E (S_x^2 - S_y^2) (see sym-SPZ4-impl.hpp). Implemented as QSZ with the spin label taken modulo 4
// (sym-QSZ-impl.hpp). <S_z> and <S_z^2> are not determined by the label and are not computed from it.
template<typename SC>
class SymmetryQSZ4 : public SymmetryQSZ<SC> {
 public:
   SymmetryQSZ4(const Params &P) : SymmetryQSZ<SC>(P, std::vector{"<Q>", "<Q^2>"}) {
     this->kind = AbelianShell::Kind::QSZ4;
     this->modulus = 4;
     initInvar({
       {"Q", additive}, // charge
       {"SSZ", mod4}    // 2S_z modulo 4
     });
   }

   void calculate_TD(const Step &step, const DiagInfo<SC> &diag, Stats<SC> &stats, const double factor) const override {
     bucket trQ, trQ2; // Tr[Q], Tr[Q^2]
     for (const auto &[I, eig]: diag) {
       const int q    = I.get("Q");
       const double sumZ = this->calculate_Z(I, eig, factor);
       trQ += sumZ * q;
       trQ2 += sumZ * q * q;
     }
     stats.td.set("<Q>",   trQ / stats.Z);
     stats.td.set("<Q^2>", trQ2 / stats.Z);
   }
};

}
