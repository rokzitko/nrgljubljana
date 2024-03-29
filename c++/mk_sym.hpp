#ifndef _mk_sym_hpp_
#define _mk_sym_hpp_

#include <memory>
#include <string>
#include <stdexcept>

#include "params.hpp"
#include "symmetry.hpp"
#include "outfield.hpp"

#include "sym-QS.hpp"
#include "sym-QSZ.hpp"

#ifdef NRG_SYM_MORE
#include "sym-ISO.hpp"
#include "sym-ISOSZ.hpp"
#include "sym-SPSU2.hpp"
#include "sym-SPU1.hpp"
#endif

#ifdef NRG_SYM_ALL
#include "sym-DBLQSZ.hpp"
#include "sym-DBLSU2.hpp"
#include "sym-DBLISOSZ.hpp"
#include "sym-ISOLR.hpp"
#include "sym-ISOSZLR.hpp"
#include "sym-NONE.hpp"
#include "sym-P.hpp"
#include "sym-PP.hpp"
#include "sym-SL.hpp"
#include "sym-SL3.hpp"
#include "sym-SPSU2LR.hpp"
#include "sym-SPSU2T.hpp"
#include "sym-SPU1LR.hpp"
#include "sym-SU2.hpp"
#include "sym-QSLR.hpp"
#include "sym-QST.hpp"
#include "sym-QSTZ.hpp"
#include "sym-QSZTZ.hpp"
#include "sym-QSZLR.hpp"
#include "sym-QJ.hpp"
#include "sym-U1.hpp"
#include "sym-SPSU2C3.hpp"
#include "sym-QSC3.hpp"
#endif

namespace NRG {

template<scalar S>
auto get(const std::string &sym_string, const Params &P)
{
  if (sym_string == "QS") return mk_QS<S>(P);
  if (sym_string == "QSZ") return mk_QSZ<S>(P);
#ifdef NRG_SYM_MORE
  if (sym_string == "ISO") return mk_ISO<S>(P);
  if (sym_string == "ISO2") return mk_ISO2<S>(P);
  if (sym_string == "ISOSZ") return mk_ISOSZ<S>(P);
  if (sym_string == "SPSU2") return mk_SPSU2<S>(P);
  if (sym_string == "SPU1") return mk_SPU1<S>(P);
#endif
#ifdef NRG_SYM_ALL
  if (sym_string == "DBLQSZ") return mk_DBLQSZ<S>(P);
  if (sym_string == "DBLSU2") return mk_DBLSU2<S>(P);
  if (sym_string == "DBLISOSZ") return mk_DBLISOSZ<S>(P);
  if (sym_string == "ISOLR") return mk_ISOLR<S>(P);
  if (sym_string == "ISO2LR") return mk_ISO2LR<S>(P);
  if (sym_string == "ISOSZLR") return mk_ISOSZLR<S>(P);
  if (sym_string == "NONE") return mk_NONE<S>(P);
  if (sym_string == "P") return mk_P<S>(P);
  if (sym_string == "PP") return mk_PP<S>(P);
  if (sym_string == "SL") return mk_SL<S>(P);
  if (sym_string == "SL3") return mk_SL3<S>(P);
  if (sym_string == "SPSU2LR") return mk_SPSU2LR<S>(P);
  if (sym_string == "SPSU2T") return mk_SPSU2T<S>(P);
  if (sym_string == "SPU1LR") return mk_SPU1LR<S>(P);
  if (sym_string == "SU2") return mk_SU2<S>(P);
  if (sym_string == "QSLR") return mk_QSLR<S>(P);
  if (sym_string == "QST") return mk_QST<S>(P);
  if (sym_string == "QSTZ") return mk_QSTZ<S>(P);
  if (sym_string == "QSZTZ") return mk_QSZTZ<S>(P);
  if (sym_string == "QSZLR") return mk_QSZLR<S>(P);
  if (sym_string == "QJ") return mk_QJ<S>(P);
  if (sym_string == "U1") return mk_U1<S>(P);
  if constexpr (std::is_same_v<S, std::complex<double>>) {
    if (sym_string == "SPSU2C3") return mk_SPSU2C3<S>(P);
    if (sym_string == "QSC3") return mk_QSC3<S>(P);
  }
#endif 
  throw std::runtime_error("Unknown symmetry " + sym_string);
}

// Called immediately after parsing the information about the number of channels from the data file. This ensures
// that Invar can be parsed correctly.
template <scalar S>
std::shared_ptr<Symmetry<S>> set_symmetry(Params &P, const std::string &sym_string, const size_t channels) {
  P.symtype = sym_string;
  P.set_channels_and_combs(channels);
  std::cout << "\nSYMMETRY TYPE: " << sym_string << std::endl;
  auto Sym = get<S>(sym_string, P);
  Sym->load();
  Sym->erase_first();
  return Sym;
}

} // namespace

#endif
