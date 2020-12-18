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
  return set_symmetry<S>(P, "QS", 1);
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
  P.absolute = true; // disable any energy rescaling
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
  P.absolute = true;
  DiagInfo<S> diag(ss, 3, P);
  return diag;
}
