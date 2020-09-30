#ifndef _read_input_cc_
#define _read_input_cc_

void set_symmetry(const string &sym_string);

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

// Read the number of channels from data file. Also sets P.combs
// accordingly, depending on the spin of the conduction band electrons.
void read_nr_channels(ifstream &fdata, std::string sym_string, Params &P) {
  size_t channels;
  fdata >> channels;
  my_assert(channels >= 1);
  P.channels = channels;
  // Number of tables of coefficients. It is doubled in the case of
  // spin-polarized conduction bands. The first half corresponds to
  // spin-up, the second half to spin-down. It is quadrupled in the
  // case of full 2x2 matrix structure in the spin space.
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
void read_gs_energy(ifstream &fdata) {
  fdata >> STAT::total_energy;
  assert_isfinite(STAT::total_energy);
}

// Read energies of initial states. nsubs is the number of subspaces.
DiagInfo read_energies(ifstream &fdata, size_t nsubs, Params &P) {
  DiagInfo diag;
  for (size_t i = 1; i <= nsubs; i++) {
    Invar I;
    size_t nrr;
    fdata >> I >> nrr;
    my_assert(nrr > 0);
    EVEC energies = EVEC(nrr);
    read_vector(fdata, energies, nrr);
    if (!P.data_has_rescaled_energies)
      energies /= SCALE(P.Ninit); // rescale to the suitable energy scale
    diag[I] = Eigen(nrr, nrr);
    diag[I].diagonal(energies);
  }
  my_assert(diag.size() == nsubs);
  return diag;
}

// Read irreducible matrix elements from stream fdata and store them in a map of matrices.
MatrixElements read_matrix_elements(ifstream &fdata, const DiagInfo &diag) {
  MatrixElements m;
  size_t nf; // Number of I1 x I2 combinations
  fdata >> nf;
  for (size_t i = 1; i <= nf; i++) {
    Invar I1, I2;
    fdata >> I1 >> I2;
    if (const auto it1 = diag.find(I1), it2 = diag.find(I2); it1 != diag.end() && it2 != diag.end())
      read_matrix(fdata, m[{I1, I2}], it1->second.getnr(), it2->second.getnr());
    else
      my_error("Corrupted input file. Stopped in read_matrix_elements()");
  }
  my_assert(m.size() == nf);
  return m;
}

void read_irreduc_f(ifstream &fdata, const DiagInfo &diag, Opch &opch, Params &P) {
  nrglog('@', "read_irreduc_f()");
  opch = Opch(P.channels);
  for (size_t i = 0; i < P.channels; i++) {
    opch[i] = OpchChannel(P.perchannel);
    for (size_t j = 0; j < P.perchannel; j++) {
      char ch;
      size_t iread, jread;
      fdata >> ch >> iread >> jread;
      my_assert(ch == 'f' && i == iread && j == jread);
      opch[i][j] = read_matrix_elements(fdata, diag);
    }
  }
}

// Determine Nmax from the length of the coefficient tables! Modify it for substeps==true. Call after
// tridiagonalization routines (if not using the tables computed by initial.m).
void determine_Nmax(Params &P) {
  size_t length_coef_table = xi.max(0); // all channels have same nr. of coefficients
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax(0)=" << P.Nmax << endl << endl;
  my_assert(length_coef_table == P.Nmax);
  if (P.substeps) P.Nmax = P.channels * P.Nmax;
  P.Nlen = P.Nmax;
  if (P.Nmax == P.Ninit) {
    cout << endl << "ZBW=true -> zero-bandwidth calculation" << endl;
    P.ZBW  = true;
    P.Nlen = P.Nmax + 1; // an additional element in the tables!
  }
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax=" << P.Nmax << endl << endl;
}

inline void skipline(ostream &F = std::cout) { F << std::endl; }

// Read all initial energies and matrix elements
std::tuple<DiagInfo, IterInfo> read_data(Params &P) {
  skipline();
  ifstream fdata("data");
  if (!fdata) my_error("Can't load initial data.");
  auto sym_string = parse_datafile_header(fdata);
  read_nr_channels(fdata, sym_string, P);
  set_symmetry(sym_string);
  read_Nmax(fdata, P);
  size_t nsubs = read_nsubs(fdata);
  skip_comments(fdata);
  DiagInfo diag0; // 0-th step of the NRG iteration
  diag0 = read_energies(fdata, nsubs, P);
  skip_comments(fdata);
  IterInfo iterinfo0;
  read_irreduc_f(fdata, diag0, iterinfo0.opch, P);
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
      case 'e': read_gs_energy(fdata); break;
      case 's': iterinfo0.ops[opname]  = read_matrix_elements(fdata, diag0); break;
      case 'p': iterinfo0.opsp[opname] = read_matrix_elements(fdata, diag0); break;
      case 'g': iterinfo0.opsg[opname] = read_matrix_elements(fdata, diag0); break;
      case 'd': iterinfo0.opd[opname]  = read_matrix_elements(fdata, diag0); break;
      case 't': iterinfo0.opt[opname]  = read_matrix_elements(fdata, diag0); break;
      case 'o': iterinfo0.opot[opname] = read_matrix_elements(fdata, diag0); break;
      case 'q': iterinfo0.opq[opname]  = read_matrix_elements(fdata, diag0); break;
      case 'z':
        xi.read(fdata);
        zeta.read(fdata);
        break;
      case 'Z':
        delta.read(fdata);
        kappa.read(fdata);
        break;
      case 'X':
        xiR.read(fdata);
        zetaR.read(fdata);
        break;
      case 'T':
        ep.read(fdata);
        em.read(fdata);
        u0p.read(fdata);
        u0m.read(fdata);
        break;
      default: my_error("Unknown block %c in data file.", ch);
    }
  }
  if (string(P.tri) == "cpp") tridiag(); // before determine_Nmax()
  determine_Nmax(P);
  return {diag0, iterinfo0};
}

#endif
