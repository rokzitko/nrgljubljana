#ifndef _read_input_hpp_
#define _read_input_hpp_

#include <memory>
#include <string>
#include <iostream>
#include <tuple>

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

// Read the number of channels from data file. Also sets P.combs accordingly, depending on the spin of the conduction
// band electrons.
inline void read_nr_channels(std::istream &fdata, const std::string &sym_string, Params &P) {
  const auto channels = read_one<size_t>(fdata);
  my_assert(channels >= 1);
  P.channels = channels;
  // Number of tables of coefficients. It is doubled in the case of spin-polarized conduction bands. The first half
  // corresponds to spin-up, the second half to spin-down. It is quadrupled in the case of full 2x2 matrix structure
  // in the spin space.
  if (P.pol2x2)
    P.coeffactor = 4;
  else if (P.polarized)
    P.coeffactor = 2;
  else
    P.coeffactor = 1;
  P.coefchannels = P.coeffactor * P.channels;
  nrglog('!', "coefchannels=" << P.coefchannels);
  if (sym_string == "U1" || sym_string == "SU2" || sym_string == "DBLSU2") {
    P.perchannel = 2; // We distinguish spin-up and spin-down operators.
  } else if (sym_string == "NONE" || sym_string == "P" || sym_string == "PP") {
    P.perchannel = 4; // We distinguish CR/AN and spin UP/DO.
  } else {
    P.perchannel = 1;
  }
  nrglog('!', "perchannel=" << P.perchannel);
  my_assert(P.perchannel >= 1);
  P. spin = sym_string == "SL" || sym_string == "SL3" ? 1 : 2;
  const int statespersite = intpow(2, P.spin);
  P.combs = !P.substeps ? intpow(statespersite, P.channels) : statespersite;
  nrglog('!', "combs=" << P.combs);
}

// Determine Nmax from the length of the coefficient tables! Modify it for substeps==true. Call after
// tridiagonalization routines (if not using the tables computed by initial.m).
template<typename S>
inline void determine_Nmax(const Coef<S> &coef, Params &P) { // Params is non-const !
  const auto length_coef_table = coef.xi.max(0); // all channels have same nr. of coefficients
  std::cout << std::endl << "length_coef_table=" << length_coef_table << " Nmax(0)=" << P.Nmax << std::endl << std::endl;
  my_assert(length_coef_table == P.Nmax);
  if (P.substeps) P.Nmax = P.channels * P.Nmax;
  P.Nlen = P.Nmax;       // this is the usual situation
  if (P.Nmax == P.Ninit) {
    std::cout << std::endl << "ZBW=true -> zero-bandwidth calculation" << std::endl;
    P.ZBW  = true;
    P.Nlen = P.Nmax + 1; // an additional element in the tables for ZBW=true
  }
  my_assert(P.Nlen < MAX_NDX);
  std::cout << std::endl << "length_coef_table=" << length_coef_table << " Nmax=" << P.Nmax << std::endl << std::endl;
}

inline void skipline(std::ostream &F = std::cout) { F << std::endl; }

// Read all initial energies and matrix elements
template<typename S> 
inline auto read_data(Params &P, std::string filename = "data") {
  skipline();
  std::ifstream fdata(filename);
  if (!fdata) throw std::runtime_error("Can't load initial data.");
  const auto sym_string = parse_datafile_header(fdata);
  my_assert(sym_string == P.symtype.value());
  read_nr_channels(fdata, sym_string, P);
  auto Sym = set_symmetry<S>(P);
  P.Nmax = read_one<size_t>(fdata); // Length of the Wilson chain
  const auto nsubs = read_one<size_t>(fdata); // Number of invariant subspaces
  skip_comments(fdata);
  DiagInfo<S> diag0(fdata, nsubs, P); // 0-th step of the NRG iteration
  skip_comments(fdata);
  Operators<S> operators0;
  operators0.opch = Opch<S>(fdata, diag0, P);
  Coef<S> coef(P);
  double GS_energy = 0.0;
  while (true) {
    /* skip white space */
    while (!fdata.eof() && std::isspace(fdata.peek())) fdata.get();
    if (fdata.eof()) break;
    char ch = fdata.get();
    std::string opname;
    std::getline(fdata, opname);
    if (ch != '#') debug("Reading <||" << opname << "||> (" << ch << ")");
    switch (ch) {
      case '#':
        // ignore embedded comment lines
        break;
      case 'e': GS_energy = read_one<double>(fdata); break;
      case 's': operators0.ops[opname]  = MatrixElements<S>(fdata, diag0); break;
      case 'p': operators0.opsp[opname] = MatrixElements<S>(fdata, diag0); break;
      case 'g': operators0.opsg[opname] = MatrixElements<S>(fdata, diag0); break;
      case 'd': operators0.opd[opname]  = MatrixElements<S>(fdata, diag0); break;
      case 't': operators0.opt[opname]  = MatrixElements<S>(fdata, diag0); break;
      case 'o': operators0.opot[opname] = MatrixElements<S>(fdata, diag0); break;
      case 'q': operators0.opq[opname]  = MatrixElements<S>(fdata, diag0); break;
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
  if (std::string(P.tri) == "cpp") Tridiag<S>(coef, P); // before calling determine_Nmax()
  determine_Nmax(coef, P);
  return std::make_tuple(Sym, diag0, operators0, coef, GS_energy);
}

} // namespace

#endif
