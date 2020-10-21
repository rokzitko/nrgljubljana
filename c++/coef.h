#ifndef _coef_h_
#define _coef_h_

// Table of Wilson chain coefficients
template <typename S> class coef_table : traits<S> {
private:
  using t_coef = typename traits<S>::t_coef;
  using t_storage = ublas::vector<t_coef>;
  t_storage table;

public:
  using t_ndx = typename t_storage::size_type;
  // Read values from a stream f
  void read_values(ifstream &f) {
    table = read_vector<t_coef>(f, true); // len=nr+1
  }
  auto coef(const t_ndx n) const {
    my_assert(n < table.size());
    return table[n];
  }
  // Index of the last coefficient still included in the table
  auto max() const {
    my_assert(table.size() >= 1);
    return table.size()-1;
  }
  void setvalue(const t_ndx n, const t_coef val) {
    if (n+1 > table.size()) table.resize(n+1);
    table[n] = val;
  }
};

// One table of discretization coefficients for each channel
template <typename S> class set_of_coef_tables : traits<S> {
private:
  using t_coef = typename traits<S>::t_coef;
  using t_storage = std::vector<coef_table<S>>;
  t_storage tabs;

public:
  using t_ndx = typename coef_table<S>::t_ndx;
  using t_ch = typename t_storage::size_type;
  auto nr_tabs() const { return tabs.size(); }
  void read(ifstream &fdata, const t_ch coefchannels) {
    tabs.resize(coefchannels);
    for (auto &i : tabs) i.read_values(fdata);
  }
  auto operator()(const t_ndx N, const t_ch alpha) const {
    my_assert(alpha < tabs.size());
    return tabs[alpha].coef(N);
  }
  auto max(const t_ch alpha) const {
    my_assert(alpha < tabs.size());
    return tabs[alpha].max();
  }
  void setvalue(const t_ndx N, const t_ch alpha, const t_coef val) {
    my_assert(alpha < tabs.size());
    tabs[alpha].setvalue(N, val);
  }
};

template <typename S> class Coef : traits<S> {
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
  auto xiUP    (const t_ndx N, const t_ch ch) const { return xi  (N, ch); }
  auto xiDOWN  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + P.channels); }
  auto zetaUP  (const t_ndx N, const t_ch ch) const { return zeta(N, ch); }
  auto zetaDOWN(const t_ndx N, const t_ch ch) const { return zeta(N, ch + P.channels); }

  // Support for conduction bands with full 2x2 matrix structure, a generalization of P.polarized. The total number
  // of "channels" is here multiplied by 4, i.e., the index runs from 0 to 4*P.channels-1. Numbers
  // 2*P.channels...3*P.channels-1 correspond to UP/DO, 3*P.channels...4*P.channels-1 correspond to DO/UP.
  auto xiUPDO  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + 2 * P.channels); }
  auto xiDOUP  (const t_ndx N, const t_ch ch) const { return xi  (N, ch + 3 * P.channels); }
  auto zetaUPDO(const t_ndx N, const t_ch ch) const { return zeta(N, ch + 2 * P.channels); }
  auto zetaDOUP(const t_ndx N, const t_ch ch) const { return zeta(N, ch + 3 * P.channels); }
};

#endif
