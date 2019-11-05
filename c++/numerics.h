// numerics.h - Miscelaneous numerical routines
// Copyright (C) 2005-2019 Rok Zitko

#ifndef _numerics_h_
#define _numerics_h_

// Accumulator abstraction: automatically initialized to 0, result
// checked for finiteness. Can be constructured from a STL vector of
// pairs, by summing the second elements.
template <typename T>
class generic_bucket
{
 private:
   T value;
 public:
   generic_bucket() { value = 0.0; }
   template <typename T1> generic_bucket(std::vector<pair<T1, T>> v) {
      value = 0.0;
      for (const auto &i : v)
          value += i.second;
   }
   inline T operator+=(T x) { return value += x; }
//   inline operator T () const { return assert_isfinite(value); } // XXXXXXXXXXXXXXXXXXXXX
   inline operator T () const { return value; }
};

typedef generic_bucket<double> bucket;
typedef generic_bucket<t_weight> weight_bucket;
typedef generic_bucket<t_matel> matel_bucket;

#define IS_ODD(n)  ((n) & 1)
#define IS_EVEN(n) (!(IS_ODD(n)))
#define SIGN(x)    ((x) >= 0.0 ? 1 : -1)

// gsl_fcmp replacement. Warning: this is not quite equivalent to the
// functionality of gsl_fcmp, but should be sufficient for our purposes.
CONSTFNC int my_fcmp(double x, double y, double epsilon)
{
   if (x == 0.0 && y == 0.0) // evidently equal
     return 0;
   // If both x and y are small, we ASSUME them to be equivalent. In this
   // context, thus, epsilon is ABSOLUTE error.
   if (abs(x) < epsilon && abs(y) < epsilon)
     return 0;
   // Here epsilon is maximum allowable RELATIVE error.
   if (abs(x-y)/(abs(x)+abs(y)) < epsilon)
     return 0;
   if (x > y)
     return +1;
   else
     return -1;
}

// Test if two numbers are equal to within numerical errors. (Use this for
// comparing values that are expected to be of order 1.)
CONSTFNC bool num_equal(double a, double b, double check_precision = 1.e-12)
{
   return my_fcmp(a, b, check_precision) == 0;
}

CONSTFNC bool num_equal(cmpl a, cmpl b, double check_precision = 1.e-12)
{
   return (my_fcmp(a.real(), b.real(), check_precision) == 0) &&
          (my_fcmp(a.imag(), b.imag(), check_precision) == 0);
}

CONSTFNC bool are_conjugate(double a, double b)
{
   return num_equal(a, b);
}

CONSTFNC bool are_conjugate(cmpl a, cmpl b)
{
   return num_equal(a.real(), b.real()) && num_equal(a.imag(), -b.imag());
}

CONSTFNC double frobenius_norm(const Matrix &m)
{
   bucket sum;
   for (size_t i = 0; i < m.size1(); i++) 
      for (size_t j = 0; j < m.size2(); j++) 
	 sum += sqr(abs(m(i,j))); 
   return sum;
}

void make_identity_matrix(Matrix &m) 
{
   my_assert(m.size1() == m.size2());
   size_t dim = m.size1();
   m.clear();
   for (size_t i = 0; i < dim; i++)
      m(i, i) = 1.0;
}

// Check if the (numeric) matrix m is indeed Hermitian.
// NOTE: current not used (8.10.2009).
CONSTFNC bool check_is_matrix_hermitian(const Matrix &m, 
					bool assert_it = false)
{
   my_assert(m.size2() == m.size1() && m.size2() >= 1);

   for (size_t i = 0; i < m.size2(); i++) 
     for (size_t j = i+1; j < m.size1(); j++) 
       if (!are_conjugate(m(i, j), m(j, i))) {
	  if (assert_it) 
	     exit1("check_is_matrix_hermitian failed");
	  return false;
       }
   return true;
}

void matrix_replicate_l_from_u(Matrix &m)
{
  my_assert(m.size1() == m.size2());
  size_t dim = m.size1();
  for (size_t i = 0; i < dim; i++) 
    for (size_t j = i+1; j < dim; j++)  // j > i
      m(j, i) = m(i, j);
}

// Check if matrix m is upper triangular. In the lower triangle, all
// elements must be 0. NOTE: we store the upper triangular part of the
// symmetric Hamiltonian matrix. In FORTRAN convention, this is the lower
// part !!

void check_is_matrix_upper(const Matrix &m)
{
  my_assert(m.size1() == m.size2() && m.size1() >= 1);
  for (size_t i = 1; i < m.size1(); i++)
    for (size_t j = 0; j < i; j++) // j < i
      my_assert(m(i, j) == 0.);
}

// Assert that two matrices are equal.
void check_are_matrices_equal(const Matrix &m1, const Matrix &m2)
{
   my_assert(m1.size1() == m2.size1());
   my_assert(m1.size2() == m2.size2());
   for (size_t i = 0; i < m1.size1(); i++)
     for (size_t j = 0; j < m1.size2(); j++)
       my_assert(m1(i, j) == m2(i, j));
}

// Dump range [jmin;jmax] of vector v to cout.
void coutloop(const DVEC &v, int jmin, int jmax)
{
  for (int j = jmin; j <= jmax; j++)
    cout << v[j] << " ";
  cout << endl;
}

// x raised to the power of n
CONSTFNC inline int pow(int x, int n)
{
  my_assert(n >= 0);
  int res = 1;
  for (int i = 1; i <= n; i++)
    res *= x;
  return res;
}

// (-1)^n
CONSTFNC inline double psgn(int n)
{
   return (n % 2 == 0 ? 1.0 : -1.0);
}

// Dump a matrix with full numerical precision. The columns
// are aligned for easier inspection. Expect large output!

void dump_matrix(const Matrix &m, ostream &fout = cout)
{
  boost::io::ios_base_all_saver ofs(fout);
  fout << setprecision(16);
  for (size_t r1 = 0; r1 < m.size1(); r1++) {
    for (size_t r2 = 0; r2 < m.size2(); r2++) 
      fout << setw(23) << m(r1, r2) << " ";
    fout << endl;
  }
}

// Chop numerical noise
template <typename T>
CONSTFNC inline T chop(T x, double xlimit = 1.e-8)
{
   return (abs(x) < xlimit ? 0.0 : x);
}

void assert_issquare(const Matrix &m)
{
  my_assert(m.size1() == m.size2());
}

// Powers, such as (-1)^n, appear in the coupling coefficients.
CONSTFNC inline double Power(double i, double nn) 
{
    return pow(i,nn); 
}

// Assert that vector "v" is (approximately) normalized to 1.
void assert_norm1(const DVEC &v)
{
   my_assert(my_fcmp(norm_2(v), 1.0, 1.e-8) == 0);
}

// Read 'len' values of type T into a vector 'vec' (which may be either an
// STL or an ublas vector). Check that values are finite and assert for no
// reading failure. If update=false, we read values without changing the
// values in the vector.
template <typename T>
void read_vector(istream &F, ublas::vector<T> &vec, 
		 int len, bool update = true)
{
  my_assert(F);
  my_assert(len >= 0);
  if (update)
     vec.resize(len);
  for (int j = 0; j < len; j++) {
    T x;
    F >> x;
    if (update) {
       assert_isfinite(x);
       vec[j] = x;
    }
  }
  if (F.fail()) 
     my_error("read_vector() error. 'data' file is corrupted.");
}

// Read 'size1' x 'size2' ublas matrix of type T.
template <typename T>
void read_matrix(istream &F, ublas::matrix<T> &m, int size1, int size2)
{
  my_assert(F);
  my_assert(size1 >= 0);
  my_assert(size2 >= 0);
  m = ublas::matrix<T>(size1, size2);
  for (int j1 = 0; j1 < size1; j1++)
    for (int j2 = 0; j2 < size2; j2++) {
      T x;
      F >> x;
      m(j1, j2) = assert_isfinite(x);
    }
  if (F.fail())
     my_error("read_matrix() error. 'data' file is corrupted.");
}

// Check if the value x is real [for complex number calculations].
const double check_real_TOLERANCE = 1e-8;
CONSTFNC inline bool is_real(t_matel x)
{
#ifdef NRG_REAL
   return true;
#else
   return abs(x.imag()) <= check_real_TOLERANCE;
#endif
}

// Check if x is real and return the real part, i.e. x.real().
CONSTFNC inline double check_real(t_matel x)
{
#ifdef NRG_REAL
   return x;
#else
   if (!is_real(x))
      cout << "Warning: expected real number, but got " << x << endl;
   return x.real();
#endif
}

// Spur (real part)
CONSTFNC double trace_real(const Matrix &m)
{
  assert_issquare(m);
  bucket sum;
  for (size_t i = 0; i < m.size2(); i++)
    sum += check_real( m(i, i) );
  return sum;
}

// As above, no check for finiteness.
CONSTFNC double trace_real_nochecks(const Matrix &m)
{
  assert_issquare(m);  
  double sum = 0.0; // enforce real
  for (size_t i = 0; i < m.size2(); i++)
    sum += check_real( m(i, i) );
  return sum;
}

cmpl csqrt(cmpl z)
{
   return sqrt(z);
}

#endif
