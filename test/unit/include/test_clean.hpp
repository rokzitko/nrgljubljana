#pragma once

// test0_clean

auto setup_P_clean(Params &P) {
  EXPECT_EQ(P.Lambda.value(), 2.0); // defaults
  EXPECT_EQ(P.discretization.value(), "Z"s);
  EXPECT_EQ(P.Ninit.value(), 0);
  P.Nmax = P.Nlen = 1;
  P.keep = 100;
  P.ops = "n_f n_f^2";
  P.specs = "n_f-n_f";
  P.finite = true;
  P.cfs = true;
  P.fdm = true;
  P.dmnrg = true;
  P.finitemats = true;
  P.fdmmats = true;
  P.dmnrgmats = true;
  P.mats = 10;
  P.T = 0.1; // temperature
}

template<typename S>
auto setup_diag_clean(Params &P, [[maybe_unused]] Symmetry<S> *Sym)
{
  std::string data =
      "-1 1\n"
      "1\n"
      "0.\n"
      "0 2\n"
      "1\n"
      "0.\n"
      "1 1\n"
      "1\n"
      "0.\n";
  std::istringstream ss(data);
  P.absolute = false; // NOTE!!
  DiagInfo<S> diag(ss, 3, P);
  return diag;
}

template<typename S>
auto setup_opch_clean(const Params &P, [[maybe_unused]] Symmetry<S> *Sym, const DiagInfo<S> &diag) {
  std::string str =
    "f 0 0\n"
    "2\n"
    "1 1 0 2\n"
    "1.4142135623730951\n"
    "0 2 -1 1\n"
    "1.\n";
  std::istringstream ss(str);
  return Opch<double>(ss, diag, P);
}

template<typename S>
auto setup_operators_clean(const DiagInfo<S> &diag)
{
  auto op = Operators<S>();
  std::string str_f = 
    "3\n"
    "-1 1 -1 1\n"
    "0.\n"
    "0 2 0 2\n"
    "1.\n"
    "1 1 1 1\n"
    "2.";
  std::istringstream ss_f(str_f);
  op.ops["n_f"] = MatrixElements<double>(ss_f, diag);

  std::string str_f2 =
    "3\n"
    "-1 1 -1 1\n"
    "0.\n"
    "0 2 0 2\n"
    "1.\n"
    "1 1 1 1\n"
    "4.";
  std::istringstream ss_f2(str_f2);
  op.ops["n_f^2"] = MatrixElements<double>(ss_f2, diag);

  return op;
}

template<typename S>
auto setup_coef_clean(const Params &P) {
  Coef<double> coef(P);
  std::string str =
    "1\n"
    "0.54528747084262258072\n"
    "0.41550946829175445321\n"
    "1\n"
    "0.e-999\n"
    "0.e-998\n";
  std::istringstream ss(str);
  coef.xi.read(ss, P.coefchannels);
  coef.zeta.read(ss, P.coefchannels);
  return coef;
}
