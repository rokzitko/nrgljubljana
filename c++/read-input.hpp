#ifndef _read_input_hpp_
#define _read_input_hpp_

#include <memory>
#include <string>
#include <iostream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "traits.hpp"
#include "constants.hpp"
#include "symmetry.hpp"
#include "params.hpp"
#include "eigen.hpp"
#include "operators.hpp"
#include "coef.hpp"
#include "mk_sym.hpp"

#include <fmt/format.h>

namespace NRG {

// Determine Nmax & Nlen, taking into account P.substeps. Must be called after the
// tridiagonalization routines if not using the tables in the data file.
template<scalar S>
void determine_Nmax_Nlen(const Coef<S> &coef, const size_t Nmax0, Params &P) { // Params is non-const !
  const auto length_coef_table = coef.xi.max(0); // all channels have the same nr. of coefficients
  my_assert(length_coef_table == Nmax0); // check consistency
  P.Nmax = !P.substeps ? Nmax0 : P.channels * Nmax0;
  P.Nlen = P.Nmax;
  if (!P.silent) std::cout << "\nlength_coef_table=" << length_coef_table << " Nmax0=" << Nmax0 << " Nmax=" << P.Nmax << "\n\n";
  my_assert(P.Nlen < MAX_NDX);
}

template<scalar S>
class InputData {
private:
  struct HeaderInfo {
    std::string sym_string;
    size_t channels;
    size_t Nmax;
    size_t nsubs;
  };

  static auto read_header(std::istream &fdata, const int expected_version = 9) {
    std::string sym_string;
    int dataversion = -1;
    while (fdata.peek() == '#') {
      fdata.ignore();
      if (fdata.peek() == '!') {
        fdata.ignore();
        fdata >> dataversion;
      } else {
        std::string line;
        std::getline(fdata, line);
        if (const auto pos = line.find("symtype", 1); pos != std::string::npos) {
          if (const auto p = line.find_last_of(" \t"); p != std::string::npos && p < line.size() - 1)
            sym_string = line.substr(p + 1);
        }
      }
      if (fdata.peek() == '\n') fdata.ignore();
    }
    my_assert(dataversion == expected_version);
    return HeaderInfo{sym_string, read_one<size_t>(fdata), read_one<size_t>(fdata), read_one<size_t>(fdata)};
  }

  static auto read_next_block(std::istream &fdata) {
    while (!fdata.eof() && std::isspace(fdata.peek())) fdata.get();
    char ch = char(fdata.get());
    std::string opname;
    std::getline(fdata, opname);
    return std::make_pair(ch, opname);
  }

  void initialize_symmetry(Params &P, const HeaderInfo &header) {
    my_assert(header.sym_string == P.symtype.value());
    Sym = set_symmetry<S>(P, header.sym_string, header.channels);
  }

  void read_seed_data(std::istream &fdata, const HeaderInfo &header, Params &P) {
    diag = DiagInfo<S>(fdata, header.nsubs, P);
    if (!P.silent)
      diag.states_report(Sym->multfnc());
    operators.opch = Opch<S>(fdata, diag, P);
  }

  void read_blocks(std::istream &fdata, Params &P) {
    while (true) {
      const auto [ch, opname] = read_next_block(fdata);
      if (fdata.eof()) break;
      if (ch != '#' && !P.silent) std::cout << "Reading <||" << opname << "||> (" << ch << ")" << std::endl;
      switch (ch) {
      case '#':
        break;
      case 'e': GS_energy = read_one<double>(fdata); break;
      case 's': operators.ops[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'p': operators.opsp[opname] = MatrixElements<S>(fdata, diag); break;
      case 'g': operators.opsg[opname] = MatrixElements<S>(fdata, diag); break;
      case 'd': operators.opd[opname]  = MatrixElements<S>(fdata, diag); break;
      case 't': operators.opt[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'o': operators.opot[opname] = MatrixElements<S>(fdata, diag); break;
      case 'q': operators.opq[opname]  = MatrixElements<S>(fdata, diag); break;
      case 'z':
        coef.xi.read(fdata, P.coefchannels);
        coef.zeta.read(fdata, P.coefchannels);
        break;
      case 'Z':
        coef.delta.read(fdata, P.coefchannels);
        coef.kappa.read(fdata, P.coefchannels);
        break;
      case 'X':
        coef.xiR.read(fdata, P.coefchannels);
        coef.zetaR.read(fdata, P.coefchannels);
        break;
      case 'T':
        coef.ep.read(fdata, P.coefchannels);
        coef.em.read(fdata, P.coefchannels);
        coef.u0p.read(fdata, P.coefchannels);
        coef.u0m.read(fdata, P.coefchannels);
        break;
#ifdef CHAIN_COEF
      case 'W':
        coef.eps.read(fdata);
        coef.t.read(fdata);
        break;
#endif
      default: throw std::invalid_argument(fmt::format("Unknown block {} in data file.", ch));
      }
    }
  }

  void finalize_coefficients(Params &P, const HeaderInfo &header) {
    if (std::string(P.tri) == "cpp") Tridiag<S>(coef, header.Nmax, P);
    determine_Nmax_Nlen(coef, header.Nmax, P);
    P.validate_after_data_file();
  }
public:
  std::shared_ptr<Symmetry<S>> Sym;
  DiagInfo<S> diag;
  Operators<S> operators;
  Coef<S> coef;
  double GS_energy = 0.0;
  InputData(Params &P, const std::string &filename = "data") : coef(P) {
    std::ifstream fdata(filename);
    if (!fdata) throw std::runtime_error("Can't load initial data.");
    const auto header = read_header(fdata);
    initialize_symmetry(P, header);
    read_seed_data(fdata, header, P);
    read_blocks(fdata, P);
    finalize_coefficients(P, header);
  };
};

} // namespace

#endif
