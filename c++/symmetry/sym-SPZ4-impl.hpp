#include <string>
#include <vector>

#include "sym-SPU1-impl.hpp"

namespace NRG {

// Z_4 spin symmetry: 2S_z conserved modulo 4, for problems with superconducting pairing and magnetic anisotropy.
// Pairing, D S_z^2, B_z, and the exchange coupling conserve S_z, while E (S_x^2 - S_y^2) = E/2 (S_+^2 + S_-^2) changes
// it by +-2; the label SSZ (2S_z + const) modulo 4 is conserved by all of these terms. Within each of the four
// subspaces the fermion parity is fixed, so the decomposition is finer than for the P symmetry type. Implemented as
// SPU1 with the label taken modulo 4 (sym-SPU1-impl.hpp). The expectation values <S_z>, <S_z^2> are not determined by
// the label and are not computed from it.
template<typename SC>
class SymmetrySPZ4 : public SymmetrySPU1<SC> {
 public:
   SymmetrySPZ4(const Params &P) : SymmetrySPU1<SC>(P, std::vector<std::string>{}) {
     this->kind = AbelianShell::Kind::SPZ4;
     this->modulus = 4;
     initInvar({
        {"SSZ", mod4} // 2S_z modulo 4
     });
   }
   void calculate_TD(const Step &, const DiagInfo<SC> &, Stats<SC> &, const double) const override {}
};

}
