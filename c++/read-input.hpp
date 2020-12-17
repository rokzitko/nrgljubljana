#ifndef _read_input_hpp_
#define _read_input_hpp_

#include <memory>
#include <string>
#include <iostream>
#include <stdexcept>
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
  std::cout << "\nSymmetry: " << sym_string << std::endl;
  return sym_string;
}

// Determine Nmax & Nlen, taking into account P.substeps. Must be called after the
// tridiagonalization routines if not using the tables in the data file.
template<scalar S>
inline void determine_Nmax_Nlen(const Coef<S> &coef, Params &P) { // Params is non-const !
  const auto length_coef_table = coef.xi.max(0); // all channels have the same nr. of coefficients
  std::cout << std::endl << "length_coef_table=" << length_coef_table << " Nmax(0)=" << P.Nmax << std::endl << std::endl;
  my_assert(length_coef_table == P.Nmax);
  if (P.ZBW()) my_assert(P.substeps == false);
  if (P.substeps) P.Nmax = P.channels * P.Nmax; // XXX: P.Nmax_steps
  P.Nlen = P.Nmax;       // this is the usual situation
  if (P.ZBW()) {
    std::cout << std::endl << "ZBW=true -> zero-bandwidth calculation" << std::endl;
    P.Nlen = P.Nmax + 1; // an additional element in the tables for ZBW
  }
  my_assert(P.Nlen < MAX_NDX);
  std::cout << std::endl << "length_coef_table=" << length_coef_table << " Nmax=" << P.Nmax << std::endl << std::endl;
}

template<scalar S>
struct initial_data {
public:
  std::string sym_string;
  size_t channels;
  size_t Nmax;
  size_t nsubs;
  DiagInfo<S> diag;
  Operators<S> operators;
  Coef<S> coef;
  double GS_energy = 0.0;
  initial_data(Params &P) : coef(P) {};
};
   
// Read all initial energies and matrix elements
template<scalar S> 
auto read_data(Params &P, const std::string &filename = "data") {
  initial_data<S> init(P);
  std::ifstream fdata(filename);
  if (!fdata) throw std::runtime_error("Can't load initial data.");
  init.sym_string = parse_datafile_header(fdata);
  init.channels   = read_one<size_t>(fdata); // Number of channels in the bath
  init.Nmax       = read_one<size_t>(fdata); // Length of the Wilson chain
  init.nsubs      = read_one<size_t>(fdata); // Number of invariant subspaces
  my_assert(init.sym_string == P.symtype.value());
  P.set_channels(init.channels);
  std::shared_ptr<Symmetry<S>> Sym = set_symmetry<S>(P);
  P.Nmax = init.Nmax;
  init.diag = DiagInfo<S>(fdata, init.nsubs, P); // 0-th step of the NRG iteration
  skip_comments(fdata);
  init.operators.opch = Opch<S>(fdata, init.diag, P);
  while (true) {
    /* skip white space */
    while (!fdata.eof() && std::isspace(fdata.peek())) fdata.get();
    if (fdata.eof()) break;
    char ch = fdata.get();
    std::string opname;
    std::getline(fdata, opname);
    if (ch != '#') std::cout << "Reading <||" << opname << "||> (" << ch << ")" << std::endl;
    switch (ch) {
      case '#':
        // ignore embedded comment lines
        break;
      case 'e': init.GS_energy = read_one<double>(fdata); break;
      case 's': init.operators.ops[opname]  = MatrixElements<S>(fdata, init.diag); break;
      case 'p': init.operators.opsp[opname] = MatrixElements<S>(fdata, init.diag); break;
      case 'g': init.operators.opsg[opname] = MatrixElements<S>(fdata, init.diag); break;
      case 'd': init.operators.opd[opname]  = MatrixElements<S>(fdata, init.diag); break;
      case 't': init.operators.opt[opname]  = MatrixElements<S>(fdata, init.diag); break;
      case 'o': init.operators.opot[opname] = MatrixElements<S>(fdata, init.diag); break;
      case 'q': init.operators.opq[opname]  = MatrixElements<S>(fdata, init.diag); break;
      case 'z':
        init.coef.xi.read(fdata, P.coefchannels);
        init.coef.zeta.read(fdata, P.coefchannels);
        break;
      case 'Z':
        init.coef.delta.read(fdata, P.coefchannels);
        init.coef.kappa.read(fdata, P.coefchannels);
        break;
      case 'X':
        init.coef.xiR.read(fdata, P.coefchannels);
        init.coef.zetaR.read(fdata, P.coefchannels);
        break;
      case 'T':
        init.coef.ep.read(fdata, P.coefchannels);
        init.coef.em.read(fdata, P.coefchannels);
        init.coef.u0p.read(fdata, P.coefchannels);
        init.coef.u0m.read(fdata, P.coefchannels);
        break;
    default: throw std::invalid_argument(fmt::format("Unknown block {} in data file.", ch));
    }
  }
  if (std::string(P.tri) == "cpp") Tridiag<S>(init.coef, P); // before calling determine_Nmax()
  determine_Nmax_Nlen(init.coef, P);
  return std::make_pair(Sym, init);
}

} // namespace

#endif
