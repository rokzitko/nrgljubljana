#ifndef _matsubara_h_
#define _matsubara_h_

enum class matstype : unsigned char { bosonic, fermionic, bb, bf, fb }; // note the explicit type

string matstypestring(matstype mt)
{
   switch (mt) {
   case matstype::bosonic:
      return "bosonic";
   case matstype::fermionic:
      return "fermionic";
   case matstype::bb:
      return "2*bosonic";
   case matstype::bf:
      return "bosonic+fermionic";
   case matstype::fb:
      return "fermionic+bosonic";
   default:
      my_assert_not_reached();
   }
}

inline double w(short n, matstype mt) // limit range to short numbers
{
   switch (mt) {
   case matstype::bosonic:
      return P::T*2.0*n*M_PI;
   case matstype::fermionic:
      return P::T*(2.0*n+1.0)*M_PI;
   default:
      my_assert_not_reached();
   }
}

inline double wb(short n)
{
   return P::Tpi*(2*n); // P::piT is const double
}

inline double wf(short n)
{
   return P::Tpi*(2*n+1);
}

class Matsubara {
private:
   typedef std::vector<pair<double,t_weight>> matsgf;
   matsgf v;
   matstype mt;
public:
   Matsubara() {}
   Matsubara(size_t mats, matstype _mt) : mt(_mt) {
      my_assert(mt == matstype::bosonic || mt == matstype::fermionic);
      for (size_t n = 0; n < mats; n++)
	 v.push_back(make_pair(w(n, mt), 0.0));
   }
   void add(size_t n, t_weight w) {
      v[n].second += w;
   }
   void save(ostream &F) const {
      F << setprecision(P::prec_xy);
      for (const auto &i : v)
	 output(F, i.first, i.second, true);
   }
   t_weight total_weight() const {
     weight_bucket sum(v);
     return sum;
   }
   friend class SpectrumMatsubara;
};

#endif // _matsubara_h_
