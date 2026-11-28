// Pairing (delta) and anomalous-hopping (kappa) Wilson-chain coefficients enter the Hamiltonian only for the
// superconducting symmetry types. Nonzero values for other types would be ignored, so they are rejected when the
// input is read (validate_pairing_coefficients in read-input.hpp).

#include <gtest/gtest.h>

#include <complex>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>

#include <coef.hpp>
#include <mk_sym.hpp>
#include <params.hpp>
#include <read-input.hpp>
#include <symmetry_traits.hpp>

using namespace NRG;

namespace {

const std::set<std::string> pairing_types{"SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2T", "SPSU2C3", "P", "PP", "NONE"};
const std::set<std::string> anomalous_types{"SPSU2", "SPU1", "SPU1LR", "SPSU2LR", "SPSU2T"};

// One coefficient set (one channel) with values for indexes 0 and 1, in the data-file format: max index, values.
template <typename S> void read_table(set_of_coef_tables<S> &table, const std::string &values) {
  std::istringstream s("1 " + values);
  table.read(s, 1);
}

// Run validate_pairing_coefficients for symmetry type 'name' with the given delta and kappa tables.
// Returns false if the symmetry type is not compiled in.
template <typename S> bool validate(const std::string &name, const std::string &delta, const std::string &kappa) {
  Params P;
  P.symtype = name;
  std::shared_ptr<Symmetry<S>> Sym;
  try {
    Sym = get<S>(name, P);
  } catch (const std::runtime_error &) { return false; }
  Coef<S> coef(P);
  read_table(coef.delta, delta);
  read_table(coef.kappa, kappa);
  validate_pairing_coefficients(coef, *Sym, P);
  return true;
}

template <typename S> int check_hooks() {
  int checked = 0;
  for (const auto &t : symmetry_traits) {
    const std::string name(t.name);
    Params P;
    try {
      const auto Sym = get<S>(name, P);
      EXPECT_EQ(Sym->uses_pairing(), pairing_types.contains(name)) << name;
      EXPECT_EQ(Sym->uses_anomalous_hopping(), anomalous_types.contains(name)) << name;
      EXPECT_EQ(Sym->complex_pairing_supported(), name != "SPSU2T") << name;
      ++checked;
    } catch (const std::runtime_error &) {} // not compiled in
  }
  return checked;
}

} // namespace

TEST(pairing_coefficients, hooks) {
  EXPECT_GE(check_hooks<double>() + check_hooks<cmpl>(), 4);
}

TEST(pairing_coefficients, zero_tables_accepted_for_all_types) {
  for (const auto &t : symmetry_traits) {
    validate<double>(std::string(t.name), "0 0", "0 0");
    validate<cmpl>(std::string(t.name), "0 0", "0 0");
  }
}

TEST(pairing_coefficients, no_pairing_block) {
  Params P;
  P.symtype = "QS";
  const auto Sym = get<double>("QS", P);
  Coef<double> coef(P); // no 'Z' block in the data file: empty tables
  EXPECT_NO_THROW(validate_pairing_coefficients(coef, *Sym, P));
}

TEST(pairing_coefficients, delta_rejected_without_pairing_term) {
  EXPECT_THROW(validate<double>("QS", "0 0.1", "0 0"), std::invalid_argument);
  EXPECT_THROW(validate<double>("QSZ", "1e-300 0", "0 0"), std::invalid_argument);
}

TEST(pairing_coefficients, kappa_rejected_without_anomalous_term) {
  EXPECT_THROW(validate<double>("QS", "0 0", "0.123 0"), std::invalid_argument);
#ifdef NRG_SYM_ALL
  EXPECT_NO_THROW(validate<double>("P", "0.1 0.1", "0 0"));
  EXPECT_THROW(validate<double>("P", "0.1 0.1", "0.123 0"), std::invalid_argument);
  EXPECT_THROW(validate<double>("NONE", "0.1 0.1", "0 0.123"), std::invalid_argument);
  EXPECT_THROW(validate<double>("PP", "0.1 0.1", "0.123 0"), std::invalid_argument);
  EXPECT_THROW(validate<cmpl>("SPSU2C3", "(0.1,0.1) (0.1,0.1)", "0.123 0"), std::invalid_argument);
#endif
}

TEST(pairing_coefficients, superconducting_types_accept_delta_and_kappa) {
#ifdef NRG_SYM_MORE
  EXPECT_NO_THROW(validate<double>("SPSU2", "0.1 0.1", "0.123 0"));
  EXPECT_NO_THROW(validate<cmpl>("SPU1", "(0.1,0.1) (0.1,0.1)", "0.123 0"));
#endif
}

TEST(pairing_coefficients, SPSU2T_requires_real_delta) {
#ifdef NRG_SYM_ALL
  EXPECT_NO_THROW(validate<cmpl>("SPSU2T", "(0.1,0) (0.1,0)", "0.123 0"));
  EXPECT_THROW(validate<cmpl>("SPSU2T", "(0.1,0) (0.1,0.1)", "0 0"), std::invalid_argument);
#endif
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
