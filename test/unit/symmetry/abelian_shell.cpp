// The coefficient tables generated at run time for the abelian symmetry types (c++/abelian_shell.hpp) must reproduce
// the tables generated with Mathematica (c++/symmetry/coefnew), entry by entry, including combination order and
// signs. Exceptions, which are corrections of the Mathematica tables, are listed in the tests below.

#include <gtest/gtest.h>

#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

#include <cmake_configure.hpp>
#include <abelian_shell.hpp>
#include <invar.hpp>

using namespace NRG;
using namespace NRG::AbelianShell;

namespace {

using Key = std::vector<long>;
using Table = std::map<Key, cmpl>;

const std::filesystem::path coefnew = std::filesystem::path(PROJECT_SOURCE_DIR) / "c++" / "symmetry" / "coefnew";

cmpl parse_value(std::string v) {
  v = std::regex_replace(v, std::regex("\\s+"), "");
  std::smatch m;
  if (std::regex_match(v, m, std::regex("Complex\\(([^,]+),([^)]+)\\)"))) return {std::stod(m[1]), std::stod(m[2])};
  return std::stod(v);
}

std::vector<long> parse_ints(const std::string &s) {
  std::vector<long> r;
  std::stringstream ss(s);
  std::string tok;
  while (std::getline(ss, tok, ',')) r.push_back(std::stol(tok));
  return r;
}

// Parse a Mathematica-generated table. Supported line forms:
//   NAME(i, j, ..., value);                 -> key (i, j, ...)
//   { i, j, value},                         -> key (i, j)              (<||f||> recalculation)
//   { i, j, Invar(..), Invar(..), value},   -> key (i, j)              (operator recalculation)
//   In[k] = Invar(a, b);                    -> key (k), value packs a + i*b (at most two quantum numbers)
Table read_table(const std::string &relpath) {
  const auto path = coefnew / relpath;
  std::ifstream f(path);
  EXPECT_TRUE(f.good()) << path;
  Table t;
  std::string line;
  const std::regex call(R"(^\s*[A-Z_0-9]+\(([0-9, ]+),\s*(Complex\([^)]*\)|[^,]+)\);\s*$)");
  const std::regex brace(R"(^\s*\{\s*([0-9]+)\s*,\s*([0-9]+)\s*,(.*),([^,]+)\},?\s*$)");
  const std::regex brace3(R"(^\s*\{\s*([0-9]+)\s*,\s*([0-9]+)\s*,([^,]+)\},?\s*$)");
  const std::regex inqn(R"(^\s*(In|QN)\[([0-9]+)\]\s*=\s*Invar\(([^)]*)\);\s*$)");
  std::smatch m;
  while (std::getline(f, line)) {
    if (line.empty()) continue;
    if (std::regex_match(line, m, inqn)) {
      const auto qn = parse_ints(m[3]);
      t[{std::stol(m[2])}] = cmpl(qn.at(0), qn.size() > 1 ? qn[1] : 0);
    } else if (std::regex_match(line, m, brace3)) {
      t[{std::stol(m[1]), std::stol(m[2])}] = parse_value(m[3]);
    } else if (std::regex_match(line, m, brace)) {
      t[{std::stol(m[1]), std::stol(m[2])}] = parse_value(m[4]);
    } else if (std::regex_match(line, m, call)) {
      t[parse_ints(m[1])] = parse_value(m[2]);
    } else {
      ADD_FAILURE() << "Unparsed line in " << relpath << ": " << line;
    }
  }
  return t;
}

Table nonzero(const Table &t) {
  Table r;
  for (const auto &[k, v] : t)
    if (std::abs(v) > 1e-14) r[k] = v;
  return r;
}

void expect_equal(const Table &expected, const Table &generated, const std::string &what) {
  const auto e = nonzero(expected);
  const auto g = nonzero(generated);
  EXPECT_EQ(e.size(), g.size()) << what;
  for (const auto &[k, v] : e) {
    const auto it = g.find(k);
    if (it == g.end()) {
      ADD_FAILURE() << what << ": missing entry " << ::testing::PrintToString(k) << " = " << v;
      continue;
    }
    EXPECT_NEAR(std::abs(it->second - v), 0.0, 1e-12) << what << " " << ::testing::PrintToString(k);
  }
  for (const auto &[k, v] : g)
    if (!e.contains(k)) ADD_FAILURE() << what << ": extra entry " << ::testing::PrintToString(k) << " = " << v;
}

std::string prefix(const std::string &sym, const std::size_t K) { return sym + "/" + sym + "-" + std::to_string(K) + "ch-"; }
char chname(const std::size_t ch) { return char('a' + ch); }

// Tables of the form NAME(i, j, ch, value) from per-channel entries.
template<typename F> Table per_channel(const Shell &s, F entries_of_channel) {
  Table t;
  for (std::size_t ch = 0; ch < s.channels(); ++ch)
    for (const auto &e : entries_of_channel(ch)) t[{long(e.i), long(e.j), long(ch)}] = e.value;
  return t;
}
template<typename E> Table pairs(const std::vector<E> &entries) {
  Table t;
  for (const auto &e : entries) t[{long(e.i), long(e.j)}] = e.value;
  return t;
}
Table in_table(const Shell &s, const bool qn) {
  Table t;
  for (std::size_t i = 1; i <= s.size(); ++i) {
    const auto &I = qn ? s.qn(i) : s.in(i);
    t[{long(i)}] = cmpl(I.getqn(0), Invar::invdim > 1 ? I.getqn(1) : 0);
  }
  return t;
}
Table diag_table(const Shell &s) {
  Table t;
  for (std::size_t ch = 0; ch < s.channels(); ++ch)
    for (std::size_t i = 1; i <= s.size(); ++i) t[{long(i), long(ch)}] = s.electrons(i, ch);
  return t;
}
Table diag_spin_table(const Shell &s, const Spin sp) {
  Table t;
  for (std::size_t ch = 0; ch < s.channels(); ++ch)
    for (std::size_t i = 1; i <= s.size(); ++i) t[{long(i), long(i), long(ch)}] = s.electrons(i, ch, sp);
  return t;
}
Operator both(const std::function<Operator(std::size_t, Spin)> &op, const std::size_t ch) { return op(ch, Spin::up) + op(ch, Spin::down); }

void check_common_operators(const Shell &s, const std::string &pre) {
  expect_equal(read_table(pre + "In2.dat"), in_table(s, false), pre + "In2");
  expect_equal(read_table(pre + "QN.dat"), in_table(s, true), pre + "QN");
  expect_equal(read_table(pre + "diag.dat"), diag_table(s), pre + "diag");
  expect_equal(read_table(pre + "diag-UP.dat"), diag_spin_table(s, Spin::up), pre + "diag-UP");
  expect_equal(read_table(pre + "diag-DOWN.dat"), diag_spin_table(s, Spin::down), pre + "diag-DOWN");
}

void init(const Kind kind) {
  switch (kind) {
    case Kind::NONE: initInvar({{"x", additive}}); break;
    case Kind::P: initInvar({{"P", multiplicative}}); break;
    case Kind::U1: initInvar({{"Q", additive}}); break;
    case Kind::QSZ: initInvar({{"Q", additive}, {"SSZ", additive}}); break;
    case Kind::SPU1: initInvar({{"SSZ", additive}}); break;
  }
}

void check_P_NONE(const Kind kind, const std::string &sym, const std::size_t K) {
  init(kind);
  const Shell s(kind, K);
  const auto pre = prefix(sym, K);
  check_common_operators(s, pre);
  for (const auto &[name, sp] : {std::pair{"UP", Spin::up}, std::pair{"DO", Spin::down}}) {
    expect_equal(read_table(pre + "offdiag-CR-" + name + ".dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, sp), true); }), pre + "offdiag-CR");
    if (kind == Kind::P)
      expect_equal(read_table(pre + "offdiag-AN-" + name + ".dat"), per_channel(s, [&](auto ch) { return s.hop_an(cdag(ch, sp), true); }), pre + "offdiag-AN");
    for (std::size_t ch = 0; ch < K; ++ch) {
      expect_equal(read_table(pre + chname(ch) + "-CR-" + name + ".dat"), pairs(s.entries(cdag(ch, sp))), pre + "f CR");
      if (kind == Kind::P) expect_equal(read_table(pre + chname(ch) + "-AN-" + name + ".dat"), pairs(s.entries(c(ch, sp))), pre + "f AN");
    }
  }
  const std::vector<std::pair<std::string, Operator (*)(std::size_t)>> ops{
    {"Ixtot", isospinx}, {"Iytot", isospiny}, {"Iztot", isospinz}, {"Iptot", isospinplus}, {"Imtot", isospinminus},
    {"Qtot", charge}, {"spinx", spinx}, {"spiny", spiny}, {"spinz", spinz}};
  for (const auto &[name, op] : ops) expect_equal(read_table(pre + name + ".dat"), per_channel(s, [&](auto ch) { return s.entries(op(ch)); }), pre + name);
  if (kind == Kind::NONE) expect_equal(read_table(pre + "doublet.dat"), pairs(s.recalc_rows(true)), pre + "doublet");
}

} // namespace

TEST(AbelianShell, NONE) {
  for (const std::size_t K : {1, 2}) check_P_NONE(Kind::NONE, "none", K);
}

TEST(AbelianShell, P) {
  for (const std::size_t K : {1, 2}) check_P_NONE(Kind::P, "p", K);
}

// The Mathematica doublet tables for P lack the fermionic sign (-1)^n_i (all factors are +1), unlike those of NONE,
// QSZ, SPU1, and U1. The run-time tables include the sign; with it, the spectral function of d agrees with QSZ.
TEST(AbelianShell, P_doublet_sign_corrected) {
  for (const std::size_t K : {1, 2}) {
    init(Kind::P);
    const Shell s(Kind::P, K);
    const auto legacy = nonzero(read_table(prefix("p", K) + "doublet.dat"));
    const auto rows = s.recalc_rows(true);
    ASSERT_EQ(legacy.size(), rows.size());
    for (const auto &[i, j, v] : rows) {
      EXPECT_EQ(legacy.at({long(i), long(j)}), 1.0);
      EXPECT_EQ(v, s.sign(i));
    }
  }
}

TEST(AbelianShell, SPU1) {
  init(Kind::SPU1);
  for (const std::size_t K : {1, 2}) {
    const Shell s(Kind::SPU1, K);
    const auto pre = prefix("spu1", K);
    check_common_operators(s, pre);
    expect_equal(read_table(pre + "offdiag.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(both(c, ch)); }), pre + "offdiag");
    expect_equal(read_table(pre + "offdiag-UP.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::up)); }), pre + "offdiag-UP");
    expect_equal(read_table(pre + "offdiag-DOWN.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::down)); }), pre + "offdiag-DOWN");
    expect_equal(read_table(pre + "anomalous.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(cdag(ch, Spin::down) + (-1.0) * cdag(ch, Spin::up)); }), pre + "anomalous");
    expect_equal(read_table(pre + "isospinx.dat"), per_channel(s, [&](auto ch) { return s.entries(isospinx(ch), true); }), pre + "isospinx");
    const std::vector<std::pair<std::string, Operator (*)(std::size_t)>> ops{
      {"Ixtot", isospinx}, {"Iztot", isospinz}, {"Iptot", isospinplus}, {"Imtot", isospinminus}, {"Qtot", charge}};
    for (const auto &[name, op] : ops) expect_equal(read_table(pre + name + ".dat"), per_channel(s, [&](auto ch) { return s.entries(op(ch)); }), pre + name);
    if (K == 2) {
      const auto only = [&](std::size_t which, Operator op) { return per_channel(s, [&](auto ch) { return ch == which ? s.entries(op) : std::vector<Entry<cmpl>>{}; }); };
      expect_equal(read_table(pre + "q1.dat"), only(0, number(0)), pre + "q1");
      expect_equal(read_table(pre + "q2.dat"), only(1, number(1)), pre + "q2");
      expect_equal(read_table(pre + "q1up.dat"), only(0, number(0, Spin::up)), pre + "q1up");
      expect_equal(read_table(pre + "q2up.dat"), only(1, number(1, Spin::up)), pre + "q2up");
      expect_equal(read_table(pre + "q1do.dat"), only(0, number(0, Spin::down)), pre + "q1do");
      expect_equal(read_table(pre + "q2do.dat"), only(1, number(1, Spin::down)), pre + "q2do");
      expect_equal(read_table(pre + "qdiff.dat"), per_channel(s, [&](auto ch) { return s.entries((ch == 1 ? 1.0 : -1.0) * number(ch)); }), pre + "qdiff");
    }
    for (std::size_t ch = 0; ch < K; ++ch) {
      expect_equal(read_table(pre + "spinup" + chname(ch) + ".dat"), pairs(s.entries(cdag(ch, Spin::up))), pre + "spinup");
      expect_equal(read_table(pre + "spindown" + chname(ch) + ".dat"), pairs(s.entries(cdag(ch, Spin::down))), pre + "spindown");
    }
    for (const auto *n : {"doubletp", "doubletm"}) expect_equal(read_table(pre + n + ".dat"), pairs(s.recalc_rows(true)), pre + n);
    for (const auto *n : {"triplets", "tripletp", "tripletm"}) expect_equal(read_table(pre + n + ".dat"), pairs(s.recalc_rows(false)), pre + n);
  }
}

TEST(AbelianShell, QSZ) {
  init(Kind::QSZ);
  for (const std::size_t K : {1, 2, 3}) {
    const Shell s(Kind::QSZ, K);
    const auto pre = prefix("qsz", K);
    check_common_operators(s, pre);
    expect_equal(read_table(pre + "offdiag.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(both(c, ch)); }), pre + "offdiag");
    expect_equal(read_table(pre + "offdiag-UP.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::up)); }), pre + "offdiag-UP");
    expect_equal(read_table(pre + "offdiag-DOWN.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::down)); }), pre + "offdiag-DOWN");
    if (K == 1) expect_equal(read_table(pre + "spinz.dat"), per_channel(s, [&](auto ch) { return s.entries(spinz(ch), true); }), pre + "spinz");
    if (K == 1) {
      const auto only0 = [&](Operator op) { return per_channel(s, [&](auto ch) { return ch == 0 ? s.entries(op, true) : std::vector<Entry<cmpl>>{}; }); };
      expect_equal(read_table(pre + "q1u.dat"), only0(number(0, Spin::up) + (-0.5) * unity()), pre + "q1u");
      expect_equal(read_table(pre + "q1d.dat"), only0(number(0, Spin::down) + (-0.5) * unity()), pre + "q1d");
    }
    if (K == 2) {
      expect_equal(read_table(pre + "offdiag-mix.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(both(c, 1 - ch)); }), pre + "offdiag-mix");
      expect_equal(read_table(pre + "runghop.dat"), pairs(s.entries(hop(0, 1))), pre + "runghop");
    }
    for (std::size_t ch = 0; ch < K; ++ch) {
      expect_equal(read_table(pre + "spinup" + chname(ch) + ".dat"), pairs(s.entries(cdag(ch, Spin::up))), pre + "spinup");
      expect_equal(read_table(pre + "spindown" + chname(ch) + ".dat"), pairs(s.entries(cdag(ch, Spin::down))), pre + "spindown");
    }
    for (const auto *n : {"doubletp", "doubletm"}) expect_equal(read_table(pre + n + ".dat"), pairs(s.recalc_rows(true)), pre + n);
    for (const auto *n : {"triplets", "tripletp", "tripletm"}) expect_equal(read_table(pre + n + ".dat"), pairs(s.recalc_rows(false)), pre + n);
  }
}

TEST(AbelianShell, U1) {
  init(Kind::U1);
  for (const std::size_t K : {1, 2, 3}) {
    const Shell s(Kind::U1, K);
    const auto pre = prefix("u1", K);
    check_common_operators(s, pre);
    expect_equal(read_table(pre + "offdiag-UP.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::up)); }), pre + "offdiag-UP");
    expect_equal(read_table(pre + "offdiag-DO.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::down)); }), pre + "offdiag-DO");
    expect_equal(read_table(pre + "offdiag-UPDO.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::down)); }), pre + "offdiag-UPDO");
    expect_equal(read_table(pre + "offdiag-DOUP.dat"), per_channel(s, [&](auto ch) { return s.hop_cr(c(ch, Spin::up)); }), pre + "offdiag-DOUP");
    expect_equal(read_table(pre + "diag-DOUP.dat"), per_channel(s, [&](auto ch) { return s.entries(product(cdag(ch, Spin::down), c(ch, Spin::up)), true); }), pre + "diag-DOUP");
    for (const auto &[name, op] : std::vector<std::pair<std::string, Operator (*)(std::size_t)>>{{"spinx", spinx}, {"spiny", spiny}, {"spinz", spinz}})
      expect_equal(read_table(pre + name + ".dat"), per_channel(s, [&](auto ch) { return s.entries(op(ch)); }), pre + name);
    for (std::size_t ch = 0; ch < K; ++ch) {
      expect_equal(read_table(pre + chname(ch) + "-UP.dat"), pairs(s.entries(cdag(ch, Spin::up))), pre + "f UP");
      expect_equal(read_table(pre + chname(ch) + "-DO.dat"), pairs(s.entries(cdag(ch, Spin::down))), pre + "f DO");
    }
    expect_equal(read_table(pre + "doublet.dat"), pairs(s.recalc_rows(true)), pre + "doublet");
  }
}

// Shells with more channels than the Mathematica tables: consistency checks.
TEST(AbelianShell, many_channels) {
  for (const auto kind : {Kind::NONE, Kind::P, Kind::U1, Kind::QSZ, Kind::SPU1}) {
    init(kind);
    for (const std::size_t K : {3, 4}) {
      const Shell s(kind, K);
      ASSERT_EQ(s.size(), std::size_t(1) << (2 * K));
      // anticommutation {c_m, c^dag_m'} = delta_mm' in the combination basis
      for (std::size_t ch = 0; ch < K; ++ch)
        for (const auto sp : both_spins) {
          const auto anti = product(c(ch, sp), cdag(ch, sp)) + product(cdag(ch, sp), c(ch, sp));
          for (std::size_t i = 1; i <= s.size(); ++i) EXPECT_EQ(s.me(i, anti, i), cmpl(1.0));
          const auto other = product(c(0, Spin::up), cdag(K - 1, Spin::down)) + product(cdag(K - 1, Spin::down), c(0, Spin::up));
          for (std::size_t i = 1; i <= s.size(); ++i)
            for (std::size_t j = 1; j <= s.size(); ++j) EXPECT_EQ(s.me(i, other, j), cmpl(0.0));
        }
      // quantum numbers: ancestor relation In = inverse of QN
      for (std::size_t i = 1; i <= s.size(); ++i) {
        auto x = s.qn(i);
        x.combine(s.in(i));
        if (kind == Kind::P) EXPECT_EQ(x, Invar(1)); else EXPECT_EQ(x, s.kind() == Kind::QSZ ? Invar(0, 0) : Invar(0));
      }
    }
  }
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
