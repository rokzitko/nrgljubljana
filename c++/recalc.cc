#ifndef _recalc_cc_
#define _recalc_cc_

struct Recalc_f {
  size_t i1; // subspace indexes
  size_t ip;
  t_factor factor;
};

// Recalculates irreducible matrix elements <I1|| f || Ip>. Called
// from recalc_irreduc() in nrg-recalc-* files.
void recalc_f(const DiagInfo &dg,
              MatrixElements &ff,
              const Invar &Ip, 
              const Invar &I1, 
              const struct Recalc_f table[],
              size_t jmax)
{
  nrglog('f', "recalc_f() ** f: (" << I1 << ") (" << Ip << ")");
   
  // TO DO: implement and *test* this for other symmetry types. 
  // TO DO: Implement this in Sym class!!
  // TO DO: For sym=QJ two types of operators with different QNs. Thus the QNs should
  // be passed to recalc_f() in order to perform such checks.
  if (sym_string == "QST" || sym_string == "SPSU2T") {
     // SYMMETRY CHECK. Important: If we return at this point, (I1, Ip)
     // combination of subspaces is not created in matrix ff
     if (!Sym->triangle_inequality(I1, Ip, Sym->Invar_f)) {
	nrglog('f', "Does not fulfill the triangle inequalities.");
	return;
     }
  }
   
  const Eigen & dgI1 = dg.find(I1)->second;
  const Eigen & dgIp = dg.find(Ip)->second;
  
  // Number of states in Ip and in I1, i.e. the dimension of the
  // <||f||> matrix of irreducible matrix elements.
  const size_t dim1 = dgI1.getnr();
  const size_t dimp = dgIp.getnr();
  nrglog('f', "dim1=" << dim1 << " dimp=" << dimp);
  if (dim1 == 0 || dimp == 0)
    return; // truncated away! ff[II] is not created.
    
  const Twoinvar II = make_pair(I1, Ip);

  ff[II] = Matrix(dim1, dimp); // new matrix of irreducible matrix elements
  Matrix & f = ff[II];
  f.clear(); // Set it to all zeros.

  // <I1||f||Ip> gets contributions from various |QSr> states. These
  // are given by i1, ip in the Recalc_f type tables.

  for (size_t j = 0; j < jmax; j++) {
      // rmax1, rmaxp are the dimensions of the invariant subspaces
      const size_t rmax1 = qsrmax[I1].rmax(table[j].i1);
      const size_t rmaxp = qsrmax[Ip].rmax(table[j].ip);

      if (!(rmax1 > 0 && rmaxp > 0)) 
	continue;
     
      const Invar & Ianc1 = a.ancestors[I1][table[j].i1];
      const Invar & Iancp = a.ancestors[Ip][table[j].ip];
     
      if (logletter('f')) {
 	nrgdump6(j, table[j].i1, table[j].ip, table[j].factor, rmax1, rmaxp);
 	nrgdump2(Ianc1, Iancp) << endl;
      }
     
      my_assert(Ianc1 == Iancp);
    
      my_assert(isfinite(table[j].factor));
      my_assert(rmax1 == rmaxp);

      const Matrix & U1 = dgI1.blocks[table[j].i1 - 1]; // offset 1.. argh!
      const Matrix & Up = dgIp.blocks[table[j].ip - 1];

      my_assert(U1.size1() == dim1 && Up.size1() == dimp);
      my_assert(U1.size2() == Up.size2()); // rmax1 == rmaxp
      
      if (U1.size2() == 0)
        my_assert_not_reached(); // ??

      my_assert(rmax1 == U1.size2());
      my_assert(rmaxp == Up.size2());
      
      // Additional sanity test: factors are in general order O(1)
      my_assert(abs(table[j].factor) < 1000);

      atlas::gemm(CblasNoTrans, CblasConjTrans, t_factor(table[j].factor),
		  U1, Up, t_factor(1.0), f);
  } // loop over j

  if (logletter('F')) {
     dump_matrix(f);
  }
}

// Structure which holds subspace information and factor for each
// of nonzero irreducible matrix elements. cf. Hofstetter PhD p. 120
// <Q+1 S+-1/2 .. i1 ||f^\dag|| Q S .. ip>_N = factor < IN1 .. ||f^\dag|| INp ..>_{N_1}
struct Recalc { 
  size_t i1; // combination of states
  size_t ip;
  Invar IN1; // subspace in N-1 stage
  Invar INp;
  t_factor factor; // additional multiplicative factor
};

// Used in recalc_general() for debugging purposes.
#define RECALC_GENERAL_DUMP \
   { nrgdump3(j, I1, Ip); \
   nrgdump2(table[j].i1, table[j].ip); \
   nrgdump2(table[j].IN1, table[j].INp); \
   nrgdump(table[j].factor) << endl; }


// We split the matrices of eigenvectors in blocks according to the
// partition into "ancestor subspaces". At the price of some copying,
// this increases memory localisation of data and thus improves
// numerical performence of gemm calls in the recalculation of matrix
// elements. Note that the original (matrix0) data is discarded after
// the splitting had completed!
void split_in_blocks_Eigen(const Invar &I, Eigen &e)
{
   e.blocks.resize(P::combs);
   
   const size_t nr = e.getnr(); // nr. of eigenpairs
   my_assert(nr > 0);

   my_assert(nr <= e.getrmax()); // rmax = length of eigenvectors
   
   for (size_t block = 0; block < P::combs; block++) {
      const size_t rmax = qsrmax[I].rmax(block + 1); // offset 1
      const size_t offset = qsrmax[I].offset(block + 1);

      my_assert(e.matrix0.size1() >= nr);
      my_assert(e.matrix0.size2() >= offset+rmax);

      matrix_range<Matrix> Up(e.matrix0, range(0, nr), 
			                 range(offset, offset+rmax));
      e.blocks[block] = Matrix(Up);
      
      my_assert(e.blocks[block].size1() == nr);
      my_assert(e.blocks[block].size2() == rmax);
   }

   e.matrix0 = Matrix(0, 0); // We don't need matrix0 anymore.
}
   
void split_in_blocks(DiagInfo &diag)
{
   LOOP(diag, i)
     split_in_blocks_Eigen(i.first, i.second);
}

/* Recalculate the (irreducible) matrix elements of various operators. This
 is the most important routine in this program, so it is heavily
 instrumentalized for debugging purposes. It is called from
 recalc_doublet(), recalc_singlet(), and other routines. The inner-most
 for() loops can be found here, so this is the right spot that one should
 try to hand optimize. */

void recalc_general(const DiagInfo &dg,
                    const MatrixElements &cold,
                    MatrixElements &cnew,
                    const Invar &I1, // target subspace (bra)
                    const Invar &Ip, // target subspace (ket)
                    const struct Recalc table[],
                    size_t jmax, // length of table
		    const Invar &Iop) // quantum numbers of the operator
{
  if (logletter('r')) {
     cout << "recalc_general: ";
     nrgdump3(I1, Ip, Iop) << endl;
  }
  
  // SYMMETRY CHECK. Important: If we return at this point, (I1, Ip)
  // combination of subspaces is not created in the matrix "cnew". If
  // triangle_inequality() malfunctions, this will trigger errors in
  // calc_trace_singlet().
  if (!Sym->triangle_inequality(I1, Ip, Iop))
     return;
   
  const Eigen & dgI1 = dg.find(I1)->second;
  const Eigen & dgIp = dg.find(Ip)->second;
  const size_t dim1 = dgI1.getnr();
  const size_t dimp = dgIp.getnr();
  const Twoinvar II = make_pair(I1, Ip);
  Matrix & cn = cnew[II] = Matrix(dim1, dimp);
  cn.clear();

  if (dim1 == 0 || dimp == 0)
    return;

  for (size_t j = 0; j < jmax; j++) { // loop over combinations of i/ip
    if (logletter('r'))
       RECALC_GENERAL_DUMP;
     
    if (!Sym->Invar_allowed(table[j].IN1) || !Sym->Invar_allowed(table[j].INp))
      continue;
       
    const size_t rmax1 = qsrmax[I1].rmax(table[j].i1);
    const size_t rmaxp = qsrmax[Ip].rmax(table[j].ip);

    // Proceed if this combination of i1/ip contributes.
    if (rmax1 == 0 || rmaxp == 0) 
      continue;

    const Invar IN1 = a.ancestors[I1] [ table[j].i1 ];
    const Invar INp = a.ancestors[Ip] [ table[j].ip ];
    my_assert(IN1 == table[j].IN1 && INp == table[j].INp);

    const Twoinvar ININ = make_pair(table[j].IN1, table[j].INp);

    const size_t cnt = cold.count(ININ); // Number of (IN1,INp) subspaces. 
    my_assert(cnt == 0 || cnt == 1); // Anything other than 0 or 1 is
                                     // a severe bug!
            
    /* There are exceptions when a subspace might not contribute. Two are
     known: 1. for triplet operators, two singlet states give no
     contribution [SU(2) symmetry]; 2. some subspaces might not exist at
     low iteration steps. If triangle inequality is not satisfied and there
     are indeed no states for the given subspace pair, this is OK and we
     just skip this case. */
           
    const bool triangle = Sym->triangle_inequality(table[j].IN1, Iop, table[j].INp);
    if (!triangle && cnt == 0)
      continue;

    /* All exceptions should be handled by now. If cnt != 1, this is
     a bug, probably related to symmetry properties and the coefficient
     tables for NRG transformations of matrices. */
    
    if (cnt != 1) { // bug trap
      RECALC_GENERAL_DUMP;
      my_error("cold.count(ININ) != 1");
    }

     // This assertion should be performed after the triangle
     // inequality test above!
     my_assert(isfinite(table[j].factor));
     
     // Additional sanity test: factors are in general order O(1)
     my_assert(abs(table[j].factor) < 1000);
     
     // If we made it this far, this subspace contributes.
     if (logletter('r')) 
	cout << "Contributes: rmax1=" << rmax1 << " rmaxp=" << rmaxp << endl;
     
    /* RECALL: rmax1 - dimension of the subspace of invariant subspace I1
     spanned by the states originating from the combination |I1>_i1, where
     i1=1...P::combs. It is clearly equal to the dimension of the invariant
     subspace IN1 from the previous (N-th) iteration. */

    // m: irreducible elements at previous stage
    const Matrix & m = cold.find(ININ)->second; 
    my_assert_equal(rmax1, m.size1());
    my_assert_equal(rmaxp, m.size2());
     
    const Matrix & U1 = dgI1.blocks[table[j].i1 - 1]; // offset 1.. argh!
    const Matrix & Up = dgIp.blocks[table[j].ip - 1];
    my_assert(U1.size1() == dim1 && U1.size2() == rmax1);
    my_assert(Up.size1() == dimp && Up.size2() == rmaxp);

    // Performace hot-spot. Ensure that you're using highly optimised BLAS
    // library. Beware of the order of library parameters when linking!!
    Matrix temp(rmax1, dimp);
    atlas::gemm(CblasNoTrans, CblasConjTrans, 
		t_factor(1.0), m, Up, t_factor(0.0), temp);
    atlas::gemm(CblasNoTrans, CblasNoTrans, 
		t_factor(table[j].factor), U1, temp, t_factor(1.0), cn);
  } // over table
   
  if (logletter('R')) {
     cout << endl;
     cout << "Matrix dump, I1=" << I1 << " Ip=" << Ip << ":" << endl;
     dump_matrix(cn);      // Dump full matrix
     cout << endl;
  }
}

// This routine is used for recalculation of global operators in
// nrg-recalc-*.cc
void recalc1_global(DiagInfo &dg,
                    const Invar &I,
                    Matrix &m,
                    size_t i1, size_t ip,
                    t_factor value)
{
  nrglog('r', "recalc1_global: " << I);

  const Eigen & dgI = dg.find(I)->second;
  const size_t dim = dgI.getnr();
  
  if (dim == 0)
    return;

  const size_t rmax1 = qsrmax[I].rmax(i1);
  const size_t rmaxp = qsrmax[I].rmax(ip);
  my_assert(rmax1 == rmaxp);

  if (rmax1 == 0 || rmaxp == 0)
     return;

  const Invar IN1 = a.ancestors[I][i1];
  const Invar INp = a.ancestors[I][ip];

  const Matrix & U1 = dgI.blocks[i1-1];
  const Matrix & Up = dgI.blocks[ip-1];
  my_assert(U1.size1() == dim && U1.size2() == rmax1);
  my_assert(Up.size1() == dim && Up.size2() == rmaxp);

  // m = m + value * U1 * Up^trans
  atlas::gemm(CblasNoTrans, CblasConjTrans, 
	      t_factor(value), U1, Up, t_factor(1.0), m);
}

#endif // _recalc_cc_
