#ifndef _read_input_cc_
#define _read_input_cc_

class Stats;

shared_ptr<Symmetry> set_symmetry(const Params &P, Stats &stats);

// Parse the header of the data file, check the version, determine the symmetry type.
std::string parse_datafile_header(istream &fdata, const int expected_version = 9)
{
  std::string sym_string = "";
  int dataversion = -1;
  while (fdata.peek() == '#') {
    fdata.ignore(); // ignore '#'
    if (fdata.peek() == '!') {
      fdata.ignore(); // ignore '!'
      fdata >> dataversion;
    } else {
      string line;
      getline(fdata, line);
      string::size_type pos = line.find("symtype", 1);
      if (pos != string::npos) {
        // Symmetry type declaration
        string::size_type p = line.find_last_of(" \t");
        if (p != string::npos && p < line.size() - 1)
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
void read_nr_channels(ifstream &fdata, std::string sym_string, Params &P) {
  size_t channels;
  fdata >> channels;
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
  if (sym_string == "SL" || sym_string == "SL3")
    P.spin = 1;
  else
    P.spin = 2;
  const int statespersite = pow(2, P.spin);
  if (!P.substeps)
    P.combs = pow(statespersite, P.channels);
  else
    P.combs = statespersite;
  nrglog('!', "combs=" << P.combs);
}

// Read the length of the Wilson chain
void read_Nmax(ifstream &fdata, Params &P) {
  size_t nmax;
  fdata >> nmax;
  P.Nmax = nmax;
}

size_t read_nsubs(ifstream &fdata)
{
  size_t nsubs; // Number of invariant subspaces
  fdata >> nsubs;
  my_assert(nsubs > 0);
  return nsubs;
}

// Read the ground state energy from data file ('e' flag)
void read_gs_energy(ifstream &fdata, Stats &stats) {
  fdata >> stats.total_energy;
}

// Determine Nmax from the length of the coefficient tables! Modify it for substeps==true. Call after
// tridiagonalization routines (if not using the tables computed by initial.m).
void determine_Nmax(const Coef &coef, Params &P) {
  size_t length_coef_table = coef.xi.max(0); // all channels have same nr. of coefficients
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax(0)=" << P.Nmax << endl << endl;
  my_assert(length_coef_table == P.Nmax);
  if (P.substeps) P.Nmax = P.channels * P.Nmax;
  P.Nlen = P.Nmax;       // this is the usual situation
  if (P.Nmax == P.Ninit) {
    cout << endl << "ZBW=true -> zero-bandwidth calculation" << endl;
    P.ZBW  = true;
    P.Nlen = P.Nmax + 1; // an additional element in the tables for ZBW=true
  }
  my_assert(P.Nlen < MAX_NDX);
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax=" << P.Nmax << endl << endl;
}

inline void skipline(ostream &F = std::cout) { F << std::endl; }

// Read all initial energies and matrix elements
std::tuple<DiagInfo, IterInfo, Coef, shared_ptr<Symmetry>> read_data(Params &P, Stats &stats) {
  skipline();
  ifstream fdata("data");
  if (!fdata) throw std::runtime_error("Can't load initial data.");
  auto sym_string = parse_datafile_header(fdata);
  my_assert(sym_string == P.symtype.value());
  read_nr_channels(fdata, sym_string, P);
  auto Sym = set_symmetry(P, stats);
  read_Nmax(fdata, P);
  size_t nsubs = read_nsubs(fdata);
  skip_comments(fdata);
  DiagInfo diag0(fdata, nsubs, P); // 0-th step of the NRG iteration
  skip_comments(fdata);
  IterInfo iterinfo0;
  iterinfo0.opch = Opch(fdata, diag0, P);
  Coef coef(P);
  while (true) {
    /* skip white space */
    while (!fdata.eof() && isspace(fdata.peek())) fdata.get();
    if (fdata.eof()) break;
    char ch = fdata.get();
    string opname;
    getline(fdata, opname);
    if (ch != '#') debug("Reading <||" << opname << "||> (" << ch << ")");
    switch (ch) {
      case '#':
        // ignore embedded comment lines
        break;
      case 'e': read_gs_energy(fdata, stats); break;
      case 's': iterinfo0.ops[opname]  = MatrixElements(fdata, diag0); break;
      case 'p': iterinfo0.opsp[opname] = MatrixElements(fdata, diag0); break;
      case 'g': iterinfo0.opsg[opname] = MatrixElements(fdata, diag0); break;
      case 'd': iterinfo0.opd[opname]  = MatrixElements(fdata, diag0); break;
      case 't': iterinfo0.opt[opname]  = MatrixElements(fdata, diag0); break;
      case 'o': iterinfo0.opot[opname] = MatrixElements(fdata, diag0); break;
      case 'q': iterinfo0.opq[opname]  = MatrixElements(fdata, diag0); break;
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
  if (string(P.tri) == "cpp") Tridiag(coef, P); // before calling determine_Nmax()
  determine_Nmax(coef, P);
  return {diag0, iterinfo0, coef, Sym};
}

#endif
