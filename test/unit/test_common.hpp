#pragma once

#include <tuple>

#include <workdir.hpp>
#include <params.hpp>
#include <mk_sym.hpp>
#include <invar.hpp>
#include <eigen.hpp>

using namespace NRG;

template<typename S>
auto setup_Sym(Params &P)
{
  P.symtype.setvalue("QS");
  P.set_channels(1);
  return set_symmetry<S>(P);
}

template<typename S>
auto setup_diag(Params &P, Symmetry<S> *Sym)
{
  std::string data =
      "0 1\n"
      "2 1 2\n"
      "1 2\n"
      "3 4 5 6\n";
  std::istringstream ss(data);
  P.absolute.setvalue(true); // disable any energy rescaling
  DiagInfo<S> diag(ss, 2, P);
  return diag;
}

template<typename S>
auto setup_diag3(Params &P, Symmetry<S> *Sym)
{
  std::string data =
      "0 1\n"
      "2 1 2\n"
      "1 2\n"
      "3 4 5 6\n"
      "2 1\n"
      "4 1 2 3 4\n";
  std::istringstream ss(data);
  P.absolute.setvalue(true);
  DiagInfo<S> diag(ss, 3, P);
  return diag;
}

// test0_clean

template<typename S>
auto setup_diag_clean(Params &P, Symmetry<S> *Sym)
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
  P.absolute.setvalue(false); // NOTE!!
  DiagInfo<S> diag(ss, 3, P);
  return diag;
}

template<typename S>
auto setup_opch_clean(const Params &P, Symmetry<S> *Sym, const DiagInfo<S> &diag) {
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
void setup_operators_clean(Operators<S> &op, const DiagInfo<S> &diag)
{
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
