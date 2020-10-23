#ifndef _mk_sym_h_
#define _mk_sym_h_

#include <memory>
#include <string>
#include "params.h"
#include "symmetry.h"

#include "sym-QS.h"
#include "sym-QSZ.h"

#ifdef NRG_SYM_MORE
#include "sym-ISO.h"
#include "sym-ISOSZ.h"
#include "sym-SPSU2.h"
#include "sym-SPU1.h"
#endif

#ifdef NRG_SYM_ALL
#include "sym-DBLSU2.h"
#include "sym-DBLISOSZ.h"
#include "sym-ISOLR.h"
#include "sym-ISOSZLR.h"
#include "sym-NONE.h"
#include "sym-P.h"
#include "sym-PP.h"
#include "sym-SL.h"
#include "sym-SL3.h"
#include "sym-SPSU2LR.h"
#include "sym-SPSU2T.h"
#include "sym-SPU1LR.h"
#include "sym-SU2.h"
#include "sym-QSLR.h"
#include "sym-QST.h"
#include "sym-QSTZ.h"
#include "sym-QSZTZ.h"
#include "sym-QSZLR.h"
#include "sym-QJ.h"
#include "sym-U1.h"
#include "sym-SPSU2C3.h"
#include "sym-QSC3.h"
#endif

template<typename S>
std::unique_ptr<Symmetry<S>> get(const std::string &sym_string, const Params &P, Allfields &allfields)
{
  if (sym_string == "QS") return mk_QS<S>(P, allfields);
  if (sym_string == "QSZ") return mk_QSZ<S>(P, allfields);
#ifdef NRG_SYM_MORE
  if (sym_string == "ISO") return mk_ISO<S>(P, allfields);
  if (sym_string == "ISO2") return mk_ISO2<S>(P, allfields);
  if (sym_string == "ISOSZ") return mk_ISOSZ<S>(P, allfields);
  if (sym_string == "SPSU2") return mk_SPSU2<S>(P, allfields);
  if (sym_string == "SPU1") return mk_SPU1<S>(P, allfields);
#endif
#ifdef NRG_SYM_ALL
  if (sym_string == "DBLSU2") return mk_DBLSU2<S>(P, allfields);
  if (sym_string == "DBLISOSZ") return mk_DBLISOSZ<S>(P, allfields);
  if (sym_string == "ISOLR") return mk_ISOLR<S>(P, allfields);
  if (sym_string == "ISO2LR") return mk_ISO2LR<S>(P, allfields);
  if (sym_string == "ISOSZLR") return mk_ISOSZLR<S>(P, allfields);
  if (sym_string == "NONE") return mk_NONE<S>(P, allfields);
  if (sym_string == "P") return mk_P<S>(P, allfields);
  if (sym_string == "PP") return mk_PP<S>(P, allfields);
  if (sym_string == "SL") return mk_SL<S>(P, allfields);
  if (sym_string == "SL3") return mk_SL3<S>(P, allfields);
  if (sym_string == "SPSU2LR") return mk_SPSU2LR<S>(P, allfields);
  if (sym_string == "SPSU2T") return mk_SPSU2T<S>(P, allfields);
  if (sym_string == "SPU1LR") return mk_SPU1LR<S>(P, allfields);
  if (sym_string == "SU2") return mk_SU2<S>(P, allfields);
  if (sym_string == "QSLR") return mk_QSLR<S>(P, allfields);
  if (sym_string == "QST") return mk_QST<S>(P, allfields);
  if (sym_string == "QSTZ") return mk_QSTZ<S>(P, allfields);
  if (sym_string == "QSZTZ") return mk_QSZTZ<S>(P, allfields);
  if (sym_string == "QSZLR") return mk_QSZLR<S>(P, allfields);
  if (sym_string == "QJ") return mk_QJ<S>(P, allfields);
  if (sym_string == "U1") return mk_U1<S>(P, allfields);
  if constexpr (std::is_same_v<S, std::complex<double>>) {
    if (sym_string == "SPSU2C3") return mk_SPSU2C3<S>(P, allfields);
    if (sym_string == "QSC3") return mk_QSC3<S>(P, allfields);
  }
#endif 
  throw std::runtime_error("Unknown symmetry " + sym_string);
}

// Called immediately after parsing the information about the number of channels from the data file. This ensures
// that Invar can be parsed correctly.
template <typename S>
std::shared_ptr<Symmetry<S>> set_symmetry(const Params &P, Stats<S> &stats) {
  my_assert(P.channels > 0 && P.combs > 0); // must be set at this point
  std::cout << "SYMMETRY TYPE: " << P.symtype.value() << std::endl;
  auto Sym = get<S>(P.symtype.value(), P, stats.td.allfields);
  Sym->load();
  Sym->erase_first();
  return Sym;
}

#endif
