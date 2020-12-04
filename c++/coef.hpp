#ifndef _coef_hpp_
#define _coef_hpp_

#include <vector>
#include <fstream>
#include "traits.hpp"
#include "numerics.hpp" // read_vector
#include "params.hpp"

namespace NRG {

// Table of Wilson chain coefficients
template <scalar S, typename t_coef = coef_traits<S>>
class coef_table {
private:
  using t_coef_v = std::vector<t_coef>;
  t_coef_v table;
public:
  using t_ndx = typename t_coef_v::size_type;
  // Read values from a stream f
  void read_values(std::istream &f) {
    table = read_std_vector<t_coef>(f, true); // len=nr+1
  }
  [[nodiscard]] auto coef(const t_ndx n) const {
    my_assert(n < table.size());
    return table[n];
  }
  // Index of the last coefficient still included in the table
  [[nodiscard]] auto max() const {
    my_assert(table.size() >= 1);
    return table.size()-1;
  }
  void setvalue(const t_ndx n, const t_coef val) {
    if (n+1 > table.size()) table.resize(n+1);
    table[n] = val;
  }
};

// One table of discretization coefficients for each channel
template <scalar S, typename t_coef = coef_traits<S>>
class set_of_coef_tables {
private:
  using t_coef_vv = std::vector<coef_table<S>>;
  t_coef_vv tabs;
public:
  using t_ndx = typename coef_table<S>::t_ndx;
  using t_ch = typename t_coef_vv::size_type;
  [[nodiscard]] auto nr_tabs() const { return tabs.size(); }
  void read(std::istream &fdata, const t_ch coefchannels) {
    tabs.resize(coefchannels);
    for (auto &i : tabs) i.read_values(fdata);
  }
  [[nodiscard]] auto operator()(const t_ndx N, const t_ch alpha) const {
    my_assert(alpha < tabs.size());
    return tabs[alpha].coef(N);
  }
  [[nodiscard]] auto max(const t_ch alpha) const {
    my_assert(alpha < tabs.size());
    return tabs[alpha].max();
  }
  void setvalue(const t_ndx N, const t_ch alpha, const t_coef val) { // used in tridiag.hpp
    my_assert(alpha < tabs.size());
    tabs[alpha].setvalue(N, val);
  }
};

template <scalar S> 
class Coef {
private:
  const Params &P;

public:
  explicit Coef(const Params &P) : P(P) {}
  using t_ndx = typename set_of_coef_tables<S>::t_ndx;
  using t_ch  = typename set_of_coef_tables<S>::t_ch;

  set_of_coef_tables<S> xi;   // f^dag_N f_N+1 terms
  set_of_coef_tables<S> zeta; // f^dag_N f_N terms

  // channel-mixing chains
  set_of_coef_tables<S> xiR;
  set_of_coef_tables<S> zetaR;

  // superconducting chains
  set_of_coef_tables<S> delta; // f^dag_up,N f^dag_down,N terms
  set_of_coef_tables<S> kappa; // f^dag_N f^dag_down,N+1 terms

  // star-representation coefficients
  set_of_coef_tables<S> ep, em;   // e_n coefficients
  set_of_coef_tables<S> u0p, u0m; // u_{0,m} coefficients

  // Support for spin-polarized conduction bands. See also P.polarized. Hack: the total number of channels is
  // doubled, the index runs from 0...2*P.channels-1. Numbers 0...P.channels-1 correspond to spin up, while
  // P.channels...2*P.channels-1 correspond to spin down. Compare P.channels and P.coefchannels (which reflects the
  // same convention in initial.m, i.e. CHANNELS vs. COEFCHANNELS).
  [[nodiscard]] auto xiUP    (const t_ndx N, const t_ch ch) const { return xi  (N, ch); }
  [[nodiscard]] auto xiDOWN  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + P.channels); }
  [[nodiscard]] auto zetaUP  (const t_ndx N, const t_ch ch) const { return zeta(N, ch); }
  [[nodiscard]] auto zetaDOWN(const t_ndx N, const t_ch ch) const { return zeta(N, ch + P.channels); }

  // Support for conduction bands with full 2x2 matrix structure, a generalization of P.polarized. The total number
  // of "channels" is here multiplied by 4, i.e., the index runs from 0 to 4*P.channels-1. Numbers
  // 2*P.channels...3*P.channels-1 correspond to UP/DO, 3*P.channels...4*P.channels-1 correspond to DO/UP.
  [[nodiscard]] auto xiUPDO  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + 2 * P.channels); }
  [[nodiscard]] auto xiDOUP  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + 3 * P.channels); }
  [[nodiscard]] auto zetaUPDO(const t_ndx N, const t_ch ch) const { return zeta(N, ch + 2 * P.channels); }
  [[nodiscard]] auto zetaDOUP(const t_ndx N, const t_ch ch) const { return zeta(N, ch + 3 * P.channels); }
};

} // namespace

#endif
