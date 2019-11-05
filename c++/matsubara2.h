#ifndef _matsubara2_h_
#define _matsubara2_h_

class Matsubara2 {
private:
   typedef matrix<t_weight> matsgf2;
   matsgf2 v;
   matstype mt;
public:
   Matsubara2() {}
   Matsubara2(size_t mats, matstype _mt) : mt(_mt) {
      my_assert(mt == matstype::bb || mt == matstype::bf || mt == matstype::fb);
      v = matsgf2(mats, mats);
      v.clear();
   }
   void add(size_t m, size_t n, t_weight w) {
      // my_assert(m < v.size1()); // uncomment in production code, significant performance penalty
      // my_assert(n < v.size2()); // idem
      v(m,n) += w;
   }
   void save(ostream &F) const {
      F << setprecision(P::prec_xy);
      matstype mt1 = (mt == matstype::fb ? matstype::fermionic : matstype::bosonic);
      matstype mt2 = (mt == matstype::bf ? matstype::fermionic : matstype::bosonic);
      for (size_t m = 0; m < v.size1(); m++)
	 for (size_t n = 0; n < v.size2(); n++)
	    output(F, w(m,mt1), w(n,mt2), v(m,n), true);
   }
   t_weight total_weight() const { 
      weight_bucket sum;
      for (size_t m = 0; m < v.size1(); ++m)
	 for (size_t n = 0; n < v.size2(); ++n)
	    sum += v(m,n);
      return sum;
   }
   friend class SpectrumMatsubara2;
};

#endif // _matsubara2_h_
