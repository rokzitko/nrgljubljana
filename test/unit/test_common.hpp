#pragma once

#include <tuple>

#include <workdir.hpp>
#include <params.hpp>
#include <mk_sym.hpp>
#include <invar.hpp>

using namespace NRG;

auto test_setup_basic()
{
  // Minimal working instantiation
  Workdir workdir("test_workdir");
  Params P("", "", workdir, true);
  P.symtype.setvalue("QS");
  P.set_channels(1);
  auto sym = set_symmetry<double>(P);
  return std::make_tuple(P, sym);
}
