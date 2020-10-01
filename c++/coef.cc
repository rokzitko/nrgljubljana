// This class holds table of generalized xi/zeta/etc. coefficients
class coef_table {
  private:
   ublas::vector<t_coef> t;

  public:
  // Read values from a stream f
  void read_values(ifstream &f) {
    t = read_vector<t_coef>(f, true); // len=nr+1
  }
  t_coef coef(size_t n) const {
    my_assert(n < t.size());
    return t[n];
  }
  // Returns the index of the last coefficient still included in the table.
  size_t max() const {
    my_assert(t.size() >= 1);
    return t.size() - 1;
  }
  void setvalue(size_t n, t_coef val) {
    if (n + 1 > t.size()) t.resize(n + 1);
    t[n] = val;
  }
};

// NOTE: One table of discretization coefficients for each channel
class set_of_tables {
  private:
  std::vector<coef_table> tabs;

  public:
  size_t nr_tabs() const { return tabs.size(); }
  void read(ifstream &fdata) {
    tabs.resize(P.coefchannels);
    for (auto &i : tabs) i.read_values(fdata);
  }
  t_coef operator()(size_t N, size_t alpha) const {
    P.allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size());
    return tabs[alpha].coef(N);
  }
  size_t max(size_t alpha) const {
    P.allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size());
    return tabs[alpha].max();
  }
  void setvalue(size_t N, size_t alpha, t_coef val) {
    P.allowed_coefchannel(alpha);
    my_assert(alpha < tabs.size() && N <= P.Nmax);
    tabs[alpha].setvalue(N, val);
  }
};

set_of_tables xi;   // f^dag_N f_N+1 terms
set_of_tables zeta; // f^dag_N f_N terms

// Support for spin-polarized conduction bands. See also P.polarized.
// Hack: the total number of channels is doubled, the index runs from
// 0...2*P.channels-1. Numbers 0...P.channels-1 correspond to spin up,
// while P.channels...2*P.channels-1 correspond to spin down.
// Compare P.channels and P.coefchannels (which reflects the same
// convention in initial.m, i.e. CHANNELS vs. COEFCHANNELS).
t_coef xiUP(size_t N, size_t ch) { return xi(N, ch); }
t_coef xiDOWN(size_t N, size_t ch) { return xi(N, ch + P.channels); }
t_coef zetaUP(size_t N, size_t ch) { return zeta(N, ch); }
t_coef zetaDOWN(size_t N, size_t ch) { return zeta(N, ch + P.channels); }

// Support for conduction bands with full 2x2 matrix structure, a
// generalization of P.polarized. The total number of "channels" is
// here multiplied by 4, i.e., the index runs from 0 to
// 4*P.channels-1. Numbers 2*P.channels...3*P.channels-1 correspond
// to UP/DO, 3*P.channels...4*P.channels-1 correspond to DO/UP.
t_coef xiUPDO(size_t N, size_t ch) { return xi(N, ch + 2 * P.channels); }
t_coef xiDOUP(size_t N, size_t ch) { return xi(N, ch + 3 * P.channels); }
t_coef zetaUPDO(size_t N, size_t ch) { return zeta(N, ch + 2 * P.channels); }
t_coef zetaDOUP(size_t N, size_t ch) { return zeta(N, ch + 3 * P.channels); }

// Support for channel-mixing Wilson chains
set_of_tables xiR;
set_of_tables zetaR;
// Support for superconducting bands.
set_of_tables delta; // f^dag_up,N f^dag_down,N terms
set_of_tables kappa; // f^dag_N f^dag_down,N+1 terms

set_of_tables ep, em;   // e_n coefficients
set_of_tables u0p, u0m; // u_{0,m} coefficients
