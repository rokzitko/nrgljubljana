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

namespace NRG {

// Parse the header of the data file, check the version, determine the symmetry type.
inline auto parse_datafile_header(std::istream &fdata, const int expected_version = 9)
{
  std::string sym_string = "";
  int dataversion = -1;
  while (fdata.peek() == '#') {
    fdata.ignore(); // ignore '#'
    if (fdata.peek() == '!') {
      fdata.ignore(); // ignore '!'
      fdata >> dataversion;
    } else {
      std::string line;
      std::getline(fdata, line);
      auto pos = line.find("symtype", 1);
      if (pos != std::string::npos) {
        // Symmetry type declaration
        auto p = line.find_last_of(" \t");
        if (p != std::string::npos && p < line.size() - 1)
          sym_string = line.substr(p + 1); // global variable
      }
    }
    if (fdata.peek() == '\n') fdata.ignore();
  }
  my_assert(dataversion == expected_version);
  return sym_string;
}

// Determine Nmax & Nlen, taking into account P.substeps. Must be called after the
// tridiagonalization routines if not using the tables in the data file.
template<scalar S>
void determine_Nmax_Nlen(const Coef<S> &coef, const size_t Nmax0, Params &P) { // Params is non-const !
  const auto length_coef_table = coef.xi.max(0); // all channels have the same nr. of coefficients
  my_assert(length_coef_table == Nmax0); // check consistency
  P.Nmax = !P.substeps ? Nmax0 : P.channels * Nmax0;
  P.Nlen = !P.ZBW() ? P.Nmax : P.Nmax+1; // an additional element in the tables for ZBW
  if (P.ZBW()) std::cout << "\nZBW=true -> zero-bandwidth calculation\n";
  std::cout << "\nlength_coef_table=" << length_coef_table << " Nmax0=" << Nmax0 << " Nmax=" << P.Nmax << "\n\n";
  my_assert(P.Nlen < MAX_NDX);
  if (P.ZBW()) my_assert(P.substeps == false);
}

inline auto get_next_block(std::istream &fdata) {
  while (!fdata.eof() && std::isspace(fdata.peek())) fdata.get();       // skip white space
  char ch = fdata.get();
  std::string opname;
  std::getline(fdata, opname);
  return std::make_pair(ch, opname);
}

template<scalar S>
class InputData {
private:
  std::string sym_string;
  size_t channels;
  size_t Nmax;
  size_t nsubs;
public:
  std::shared_ptr<Symmetry<S>> Sym;
  DiagInfo<S> diag;
  Operators<S> operators;
  Coef<S> coef;
  double GS_energy = 0.0;
  InputData(Params &P, const std::string &filename = "data") : coef(P) {
    std::ifstream fdata(filename);
    if (!fdata) throw std::runtime_error("Can't load initial data.");
    sym_string = parse_datafile_header(fdata);
    channels   = read_one<size_t>(fdata); // Number of channels in the bath
    Nmax       = read_one<size_t>(fdata); // Length of the Wilson chain
    nsubs      = read_one<size_t>(fdata); // Number of invariant subspaces
    my_assert(sym_string == P.symtype.value()); // Check consistency between 'param' and 'data' file
    Sym = set_symmetry<S>(P, sym_string, channels);
    diag = DiagInfo<S>(fdata, nsubs, P); // 0-th step of the NRG iteration
    operators.opch = Opch<S>(fdata, diag, P);
    while (true) {
      const auto [ch, opname] = get_next_block(fdata);
      if (fdata.eof()) break;
      if (ch != '#') std::cout << "Reading <||" << opname << "||> (" << ch << ")" << std::endl;
      switch (ch) {
        case '#':
          break; // ignore embedded comment lines
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
      default: throw std::invalid_argument(fmt::format("Unknown block {} in data file.", ch));
      }
    }
    if (std::string(P.tri) == "cpp") Tridiag<S>(coef, Nmax, P); // before calling determine_Nmax_Nlen()
    determine_Nmax_Nlen(coef, Nmax, P);
  };
};

} // namespace

#endif
