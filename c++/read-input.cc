#ifndef _read_input_cc_
#define _read_input_cc_

// Check if the parameters make sense.
// Also initialize maps related to parameters (switch for string type).
void validateparameters()
{
  my_assert(P::keep > 1);
  if (P::keepenergy > 0.0) 
    my_assert(P::keepmin <= P::keep);
  if (P::dmnrg || cfs_flags())
    P::dm.setvalue(true);
  if (cfs_flags() && ! P::lastalloverride)
    P::lastall.setvalue(true);
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
  if (P::diagroutine == undefined)
    my_error("Unknown diagonalization routine.");
  if (P::diagroutine == diagdsyevr || P::diagroutine == diagzheevr) {
    my_assert(0.0 < P::diagratio && P::diagratio <= 1.0);
    if (cfs_flags() && P::diagratio != 1.0)
       my_error("CFS/FDM is not compatible with partial diagonalisation.");
  }
  // dumpabs=true and dumpscaled=true is a meaningless combination
  my_assert(!(P::dumpabs && P::dumpscaled));
}

// Read the length of the Wilson chain
void read_Nmax(ifstream &fdata)
{
   size_t nmax;
   fdata >> nmax;
   if (P::Nmax == 0)
      P::Nmax = nmax;
}

// Read the number of channels from data file. Also sets P::combs
// accordingly, depending on the spin of the conduction band electrons.
void read_nr_channels(ifstream &fdata) 
{
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
  if (sym_string == "SL" || sym_string == "SL3") {
    P::spin.setvalue(1);
  } else if (sym_string != "ANYJ")
    my_assert(P::spin == 2);
  nrglog('!', "spin=" << P::spin);
  my_assert(P::spin >= 1);
  const int statespersite = pow(2, P::spin);
  if (!P::substeps) 
     P::combs = pow(statespersite, P::channels);
   else
     P::combs = statespersite;
  nrglog('!', "combs=" << P::combs);
}

// Read the ground state energy from data file ('e' flag)
void read_gs_energy(ifstream &fdata)
{
  fdata >> STAT::totalenergy;
  assert_isfinite(STAT::totalenergy);
}

// Read energies of initial states. nsubs is the number of subspaces.
void read_energies(ifstream &fdata, DiagInfo &diag, size_t nsubs)
{
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
  my_assert( diag.size() == nsubs );
}

// Read irreducible matrix elements from stream fdata and store them in a
// map of matrices m. The format is specified in file "data.spec"
void read_matrix_elements(ifstream &fdata, MatrixElements &m, DiagInfo &dg)
{
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
  my_assert( m.size() == nf );
}

void read_ireducf(ifstream &fdata, DiagInfo &dg)
{
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

#endif
