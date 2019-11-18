// symmetry.cc - Classes representing various symmetry types
// Copyright (C) 2009-2019 Rok Zitko

#ifndef _symmetry_cc_
#define _symmetry_cc_

// some forward declarations
CONSTFNC double calculate_Z(const DiagInfo::value_type &is, double factor);
namespace STAT {
  double Z; // statistical sum (at shell n)
}

// In and QN are global variables. They contain information about
// how the invariant subspaces at consecutive iteration steps are
// combined.

/* In is the array of quantum number DIFFERENCES used in the construction
   of the basis. QN is the array of the conserved quantum numbers
   corresponding to the states being added. For example, in case of SU(2)
   symmetry, In will include S_z, while QN will include S. In other words,
   QN involves those quantum numbers that we need to retain in the
   calculation, while In involves those quantum numbers that "drop out" of
   the problem due to the symmetry. */

std::vector<Invar> In, QN;

// Check if the triangle inequality is satisfied (i.e. if
// Clebsch-Gordan coefficient can be different from zero). This is
// important, for example, for triplet operators, which are zero when
// evaluated between two singlet states.
// Arguments ss1, ss2, ss3 are spin multiplicities.
// Returns true if the inequality is satisfied, false otherwise.
bool su2_triangle_inequality(int ss1, int ss2, int ss3) {
  return (abs(ss1 - ss2) <= ss3 - 1) && (abs(ss2 - ss3) <= ss1 - 1) && (abs(ss3 - ss1) <= ss2 - 1);
}

// The equality for U(1) symmetry.
bool u1_equality(int q1, int q2, int q3) { return q1 == q2 + q3; }

bool z2_equality(int p1, int p2, int p3) { return p1 == p2 * p3; }

// C_3 quantum number: Equality modulo 3
bool c3_equality(int p1, int p2, int p3) { return p1 == (p2 + p3) % 3; }

class Symmetry;

// List of all symmetries that are compiled-in, indexed by the
// symmtry name string.
typedef map<string, Symmetry *> sym_map;
sym_map all_syms;

class Symmetry {
  protected:
  size_t channels{};
  bool substeps{};
  int combs{};

  public:
  Symmetry() = default;
  virtual ~Symmetry()= default;;

  virtual void init() { my_error("Bug: Initializer must be defined!"); };

  Invar InvarSinglet; // QNs for singlet operator
  Invar Invar_f;      // QNs for f operator

  void set_channels(size_t ch) { channels = ch; }
  void set_substeps(bool sb) { substeps = sb; }

  // For some symmetry types with two-channels we distinguish between
  // even and odd parity with respect to the channel-interchange
  // operation.
  virtual bool islr() { return false; }

  // Ditto for 3 channels: C_3 symmetry.
  virtual bool isc3() { return false; }

  // For some symmetry types, we may distinguish between spin-up
  // and spin-down quantities (in particular spin-up and spin-down
  // spectral functions).
  virtual bool isfield() { return false; }

  // Multiplicity of the states in the invariant subspace
  virtual int mult(const Invar &) { // XXX: should be unsigned, probably size_t
    return 1;
  };

  // Does the combination of invariant subspaces I1 and I2 contribute
  // to the spectral function corresponding to spin SPIN?
  virtual bool check_SPIN(const Invar &I1, const Invar &I2, const int &SPIN) {
    my_assert(SPIN == 0);
    return true;
  }

  // Is the triangle inequality satisfied for I1, I2, and I3 (i.e. if
  // Clebsch-Gordan coefficient can be different from zero).  This is
  // important, for example, for triplet operators which are zero
  // when evaluated between two singlet states.
  virtual bool triangle_inequality(const Invar &I1, const Invar &I2, const Invar &I3) { return true; }

  void set_combs(int _combs) {
    combs = _combs;
    In.resize(combs + 1);
    QN.resize(combs + 1);
  }

  // Setup the combinations of quantum numbers that are used in the
  // construction of the Hamiltonian matrix.
  virtual void load() = 0;

  void report() {
    if (!logletter('Q')) return;
    for (size_t i = 1; i <= P::combs; i++) cout << "In[" << i << "]=(" << In[i] << ")" << endl;
    for (size_t i = 1; i <= P::combs; i++) cout << "QN[" << i << "]=(" << QN[i] << ")" << endl;
  }

  // Is an invariant subspace with given quantum numbers allowed?
  virtual bool Invar_allowed(const Invar &I) { return true; }

  virtual void makematrix(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In) = 0;

  // Called from recalc_dynamicsusceptibility().  This is the factor due
  // to the spin degeneracy when calculating the trace of Sz.Sz.
  virtual double dynamicsusceptibility_factor(const Invar &Ip, const Invar &I1) { return 1.0; }

  // Called from recalc_dynamic_orb_susceptibility().  This is the factor due
  // to the orbital moment degeneracy when calculating the trace of Tz.Tz.
  virtual double dynamic_orb_susceptibility_factor(const Invar &Ip, const Invar &I1) { return 1.0; }

  // Called from calc_specdens().
  // See spectral_density_clebschgordan.nb and DMNRG_clebschgordan.nb.
  virtual double specdens_factor(const Invar &Ip, const Invar &I1) { return 1.0; }

  virtual double specdensquad_factor(const Invar &Ip, const Invar &I1) { return 1.0; }

  virtual void calculate_TD(const DiagInfo &diag, double factor) = 0;

  virtual void recalc_irreduc(const DiagInfo &diag) { my_error("Not implemented."); }

  virtual void recalc_irreduc_substeps(const DiagInfo &diag, int M) { my_error("Not implemented."); }

  virtual void recalc_doublet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) { my_error("Not implemented."); }

  virtual void recalc_triplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) { my_error("Not implemented."); }

  virtual void recalc_orb_triplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) { my_error("Not implemented."); }

  virtual void recalc_quadruplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) { my_error("Not implemented."); }

  virtual void recalc_global(DiagInfo &diag, string name, MatrixElements &cnew) { my_error("Not implemented."); }

  virtual void show_coefficients() {}
};

Symmetry *Sym = nullptr;

inline size_t mult(const Invar &I) { return Sym->mult(I); }

// Add DECL declaration in each symmetry class
#define DECL                                                                                                                                         \
  void makematrix(Matrix &h, const Rmaxvals &qq, const Invar &I, const InvarVec &In) override;                                                                \
  void recalc_irreduc(const DiagInfo &diag) override

// Optional declaration
#define HAS_DOUBLET void recalc_doublet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) override
#define HAS_TRIPLET void recalc_triplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) override
#define HAS_ORB_TRIPLET void recalc_orb_triplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) override
#define HAS_QUADRUPLET void recalc_quadruplet(DiagInfo &diag, MatrixElements &cold, MatrixElements &cnew) override
#define HAS_GLOBAL void recalc_global(DiagInfo &diag, string name, MatrixElements &cnew) override
#define HAS_SUBSTEPS void recalc_irreduc_substeps(const DiagInfo &diag, int M) override

class SymField : public Symmetry {
  public:
  SymField() : Symmetry(){};
  bool isfield() override { return true; }
};

class SymLR : public Symmetry {
  public:
  SymLR() : Symmetry(){};
  bool islr() override { return true; }
};

class SymC3 : public Symmetry {
  public:
  SymC3() : Symmetry(){};
  bool isc3() override { return true; }
};

class SymFieldLR : public Symmetry {
  public:
  SymFieldLR() : Symmetry(){};
  bool isfield() override { return true; }
  bool islr() override { return true; }
};

// Helper functions
void check_abs_diff(const Invar &Ip, const Invar &I1, const string &what, int diff) {
  const int a = Ip.get(what);
  const int b = I1.get(what);
  my_assert(abs(b - a) == diff);
}

void check_diff(const Invar &Ip, const Invar &I1, const string &what, int diff) {
  const int a = Ip.get(what);
  const int b = I1.get(what);
  my_assert(b - a == diff);
}

#endif // _symmetry_cc_
