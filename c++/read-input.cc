#ifndef _read_input_cc_
#define _read_input_cc_

void set_symmetry(const string &sym_string);

// Parse the header of the data file, check the version, determine the symmetry type.
void parse_datafile_header(istream &fdata, int expected_version = 9)
{
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
}

// Read the number of channels from data file. Also sets P::combs
// accordingly, depending on the spin of the conduction band electrons.
void read_nr_channels(ifstream &fdata) {
  size_t channels;
  fdata >> channels;
  my_assert(channels >= 1);
  P::channels = channels;
  // Number of tables of coefficients. It is doubled in the case of
  // spin-polarized conduction bands. The first half corresponds to
  // spin-up, the second half to spin-down. It is quadrupled in the
  // case of full 2x2 matrix structure in the spin space.
  if (P::pol2x2)
    P::coeffactor = 4;
  else if (P::polarized)
    P::coeffactor = 2;
  else
    P::coeffactor = 1;
  P::coefchannels = P::coeffactor * P::channels;
  nrglog('!', "coefchannels=" << P::coefchannels);
  if (sym_string == "U1" || sym_string == "SU2" || sym_string == "DBLSU2") {
    P::perchannel = 2; // We distinguish spin-up and spin-down operators.
  } else if (sym_string == "NONE" || sym_string == "P" || sym_string == "PP") {
    P::perchannel = 4; // We distinguish CR/AN and spin UP/DO.
  } else {
    P::perchannel = 1;
  }
  nrglog('!', "perchannel=" << P::perchannel);
  my_assert(P::perchannel >= 1);
  if (sym_string == "SL" || sym_string == "SL3")
    P::spin = 1;
  else
    P::spin = 2;
  const int statespersite = pow(2, P::spin);
  if (!P::substeps)
    P::combs = pow(statespersite, P::channels);
  else
    P::combs = statespersite;
  nrglog('!', "combs=" << P::combs);
}

// Read the length of the Wilson chain
void read_Nmax(ifstream &fdata) {
  size_t nmax;
  fdata >> nmax;
  P::Nmax = nmax;
}

size_t read_nsubs(ifstream &fdata)
{
  size_t nsubs; // Number of invariant subspaces
  fdata >> nsubs;
  my_assert(nsubs > 0);
  return nsubs;
}

// Check if the parameters make sense.
// Also initialize maps related to parameters (switch for string type).
void validateparameters() {
  my_assert(P::keep > 1);
  if (P::keepenergy > 0.0) my_assert(P::keepmin <= P::keep);
  if (P::dmnrg || cfs_flags()) P::dm.setvalue(true);
  if (cfs_flags() && !P::lastalloverride) P::lastall.setvalue(true);
  my_assert(P::Lambda > 1.0);
  P::diagroutine = undefined;
  if (NRG::v == "real") {
    if (string(P::diag) == "default") P::diag.setvalue("dsyev");
    if (string(P::diag) == "dsyev") P::diagroutine = diagdsyev;
    if (string(P::diag) == "dsyevr") P::diagroutine = diagdsyevr;
  }
  if (NRG::v == "complex") {
    if (string(P::diag) == "default") P::diag.setvalue("zheev");
    if (string(P::diag) == "zheev") P::diagroutine = diagzheev;
    if (string(P::diag) == "zheevr") P::diagroutine = diagzheevr;
  }
  if (P::diagroutine == undefined) my_error("Unknown diagonalization routine.");
  if (P::diagroutine == diagdsyevr || P::diagroutine == diagzheevr) {
    my_assert(0.0 < P::diagratio && P::diagratio <= 1.0);
    if (cfs_flags() && P::diagratio != 1.0) my_error("CFS/FDM is not compatible with partial diagonalisation.");
  }
  // dumpabs=true and dumpscaled=true is a meaningless combination
  my_assert(!(P::dumpabs && P::dumpscaled));
}

// Read the ground state energy from data file ('e' flag)
void read_gs_energy(ifstream &fdata) {
  fdata >> STAT::totalenergy;
  assert_isfinite(STAT::totalenergy);
}

// Read energies of initial states. nsubs is the number of subspaces.
void read_energies(ifstream &fdata, DiagInfo &diag, size_t nsubs) {
  diag.clear(); // clean-up first!
  for (size_t i = 1; i <= nsubs; i++) {
    Invar I;
    size_t nrr;
    fdata >> I >> nrr;
    my_assert(nrr > 0);
    EVEC energies = EVEC(nrr);
    read_vector(fdata, energies, nrr);
    diag[I] = Eigen(nrr, nrr);
    diag[I].diagonal(energies);
  }
  my_assert(diag.size() == nsubs);
}

// Read irreducible matrix elements from stream fdata and store them in a
// map of matrices m. The format is specified in file "data.spec"
void read_matrix_elements(ifstream &fdata, MatrixElements &m, DiagInfo &dg) {
  m.clear();
  size_t nf; // Number of I1 x I2 combinations
  fdata >> nf;
  for (size_t i = 1; i <= nf; i++) {
    Invar I1, I2;
    fdata >> I1 >> I2;
    const size_t size1 = dg[I1].getnr();
    const size_t size2 = dg[I2].getnr();
    read_matrix(fdata, m[make_pair(I1, I2)], size1, size2);
  }
  my_assert(m.size() == nf);
}

void read_ireducf(ifstream &fdata, DiagInfo &dg) {
  a.opch = Opch(P::channels);
  for (size_t i = 0; i < P::channels; i++) {
    a.opch[i] = OpchChannel(P::perchannel);
    for (size_t j = 0; j < P::perchannel; j++) {
      char ch;
      size_t iread, jread;
      fdata >> ch >> iread >> jread;
      my_assert(ch == 'f' && i == iread && j == jread);
      read_matrix_elements(fdata, a.opch[i][j], diagprev);
    }
  }
}

// Misc checks for validity of the input parameters and data file.
void check_validity() {
  if (P::substeps) my_assert(sym_string == "QS" || sym_string == "QSZ" || sym_string == "SPSU2" || sym_string == "SPU1");
}

// Determine Nmax from the length of the coefficient tables! Modify
// it for substeps==true. Call after tridiagonalization routines (if
// not using the tables computed by initial.m).
void determine_Nmax() {
  size_t length_coef_table = xi.max(0); // all channels have same nr. of coefficients
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax(0)=" << P::Nmax << endl << endl;
  my_assert(length_coef_table == P::Nmax);
  if (P::substeps) P::Nmax = P::channels * P::Nmax;
  P::Nlen = P::Nmax;
  if (P::Nmax == P::Ninit) {
    cout << endl << "ZBW=true -> zero-bandwidth calculation" << endl;
    P::ZBW  = true;
    P::Nlen = P::Nmax + 1; // an additional element in the tables!
  }
  cout << endl << "length_coef_table=" << length_coef_table << " Nmax=" << P::Nmax << endl << endl;
}

// Read all initial energies and matrix elements
void read_data() {
  cout << endl;
  ifstream fdata("data");
  if (!fdata) my_error("Can't load initial data.");
  parse_datafile_header(fdata);
  read_nr_channels(fdata);
  set_symmetry(sym_string);
  read_Nmax(fdata);
  size_t nsubs = read_nsubs(fdata);
  skip_comments(fdata);
  // Note: we are reading diagprev, not diag, since this is
  // information from the previous (0-th) NRG step
  read_energies(fdata, diagprev, nsubs);
  skip_comments(fdata);
  read_ireducf(fdata, diagprev);
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
      case 's': read_matrix_elements(fdata, a.ops[opname], diagprev); break;
      case 'p': read_matrix_elements(fdata, a.opsp[opname], diagprev); break;
      case 'g': read_matrix_elements(fdata, a.opsg[opname], diagprev); break;
      case 'd': read_matrix_elements(fdata, a.opd[opname], diagprev); break;
      case 't': read_matrix_elements(fdata, a.opt[opname], diagprev); break;
      case 'o': read_matrix_elements(fdata, a.opot[opname], diagprev); break;
      case 'q': read_matrix_elements(fdata, a.opq[opname], diagprev); break;
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
  check_validity();
  if (string(P::tri) == "cpp") tridiag(); // before determine_Nmax()
  determine_Nmax();
}

#endif
