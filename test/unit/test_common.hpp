#pragma once

#include <tuple>

#include <workdir.hpp>
#include <params.hpp>
#include <mk_sym.hpp>
#include <invar.hpp>

using namespace NRG;

// Minimal working instantiations of increasing complexity

auto test_setup_P()
{
  Workdir workdir(".");
  Params P("", "", workdir, true);
  return P;
}

auto test_setup_basic()
{
  auto P = test_setup_P();
  P.symtype.setvalue("QS");
  P.set_channels(1);
  auto sym = set_symmetry<double>(P);
  return std::make_tuple(P, sym);
}
