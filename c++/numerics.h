// numerics.h - Miscelaneous numerical routines
// Copyright (C) 2005-2020 Rok Zitko

#ifndef _numerics_h_
#define _numerics_h_

template <typename T>
  using complex_array_ref_t = T(&)[2];

template<typename T>
  complex_array_ref_t<T> reim(std::complex<T>& z) {
    return reinterpret_cast<T(&)[2]>(z);
  }

template <typename T>
  using complex_array_const_ref_t = const T(&)[2];

template<typename T>
  complex_array_const_ref_t<T> reim(const std::complex<T>& z) {
    return reinterpret_cast<const T(&)[2]>(z);
  }

// Accumulator abstraction: automatically initialized to 0, result checked for finiteness. Can be constructured from
// a STL vector of pairs, by summing the second elements.
template <typename T> class generic_bucket {
private:
  T value;
  
public:
  generic_bucket() : value{} {}
  template <typename T1> generic_bucket(std::vector<pair<T1, T>> v) : value{} {
    for (const auto &i : v) value += i.second;
  }
  inline T operator+=(T x) { return value += x; }
  inline operator T() const { return value; }
};

using bucket = generic_bucket<double>;
using weight_bucket = generic_bucket<t_weight>;
using matel_bucket = generic_bucket<t_matel>;

#define IS_ODD(n) ((n)&1)
#define IS_EVEN(n) (!(IS_ODD(n)))

CONSTFNC int my_fcmp(const double x, const double y, const double small_epsilon, const double rel_epsilon) {
  if (x == 0.0 && y == 0.0) return 0.0; // evidently equal
  if (abs(x) < small_epsilon && abs(y) < small_epsilon) return 0; // If both x and y are small, we ASSUME them to be equivalent
  if (abs(x-y) < rel_epsilon * (abs(x)+abs(y))) return 0;
  return boost::math::sign(x-y);
}

CONSTFNC int my_fcmp(const double x, const double y, const double epsilon) { return my_fcmp(x, y, epsilon, epsilon); }

// Test if two numbers are equal to within numerical errors. (Use this for comparing values that are expected to be
// of order 1.)
CONSTFNC bool num_equal(const double a, const double b, const double check_precision = 1.e-12) { 
  return my_fcmp(a, b, check_precision) == 0; 
}

CONSTFNC bool num_equal(const cmpl a, const cmpl b, const double check_precision = 1.e-12) {
  return (my_fcmp(a.real(), b.real(), check_precision) == 0) && (my_fcmp(a.imag(), b.imag(), check_precision) == 0);
}

CONSTFNC bool are_conjugate(const double a, const double b) { return num_equal(a, b); }

CONSTFNC bool are_conjugate(const cmpl a, const cmpl b) { return num_equal(a.real(), b.real()) && num_equal(a.imag(), -b.imag()); }

CONSTFNC double frobenius_norm(const Matrix &m) { // Frobenius norm (without taking the final square root!)
  bucket sum;
  for (auto i = 0; i < m.size1(); i++)
    for (auto j = 0; j < m.size2(); j++) sum += sqr(abs(m(i, j)));
  return sum;
}

// Check if matrix m is upper triangular. In the lower triangle, all elements must be 0. NOTE: we store the upper
// triangular part of the symmetric Hamiltonian matrix. In FORTRAN convention, this is the lower part !!

void check_is_matrix_upper(const Matrix &m) {
  my_assert(m.size1() == m.size2() && m.size1() >= 1);
  for (auto i = 1; i < m.size1(); i++)
    for (auto j = 0; j < i; j++) // j < i
      my_assert(m(i, j) == 0.);
}

// x raised to the power of n
CONSTFNC inline int pow(const int x, const int n) {
  my_assert(n >= 0);
  int res = 1;
  for (int i = 1; i <= n; i++) res *= x;
  return res;
}

// (-1)^n
CONSTFNC inline double psgn(const int n) { return (n % 2 == 0 ? 1.0 : -1.0); }

// Dump a matrix with full numerical precision. The columns
// are aligned for easier inspection. Expect large output!

void dump_matrix(const Matrix &m, ostream &fout = cout) {
  boost::io::ios_base_all_saver ofs(fout);
  fout << setprecision(std::numeric_limits<double>::max_digits10);
  for (auto r1 = 0; r1 < m.size1(); r1++) {
    for (auto r2 = 0; r2 < m.size2(); r2++) fout << setw(23) << m(r1, r2) << " ";
    fout << endl;
  }
}

// Chop numerical noise
template <typename T> CONSTFNC inline T chop(const T x, const double xlimit = 1.e-8) { return abs(x) < xlimit ? 0.0 : x; }

void assert_issquare(const Matrix &m) { my_assert(m.size1() == m.size2()); }

// Powers, such as (-1)^n, appear in the coupling coefficients.
CONSTFNC inline double Power(const double i, const double nn) { return pow(i, nn); }

// Read 'len' values of type T into a ublas vector<T>.
template <typename T> 
  ublas::vector<T> read_vector(istream &F, const bool nr_is_max_index = false) {
    my_assert(F);
    size_t nr;
    F >> nr;
    // nr is either vector dimension or the value of maximum index
    auto len = nr_is_max_index ? nr+1 : nr;
    ublas::vector<T> vec(len);
    for (auto j = 0; j < len; j++)
      F >> vec[j];
    if (F.fail()) throw std::runtime_error("read_vector() error. Input file is corrupted.");
    return vec;
  }

// Read 'size1' x 'size2' ublas matrix of type T.
template <typename T> void read_matrix(istream &F, ublas::matrix<T> &m, const size_t size1, const size_t size2) {
  my_assert(F);
  m = ublas::matrix<T>(size1, size2);
  for (auto j1 = 0; j1 < size1; j1++)
    for (auto j2 = 0; j2 < size2; j2++) {
      T x;
      F >> x;
      m(j1, j2) = assert_isfinite(x);
    }
  if (F.fail()) std::runtime_error("read_matrix() error. 'data' file is corrupted.");
}

// Check if the value x is real [for complex number calculations].
const double check_real_TOLERANCE = 1e-8;
CONSTFNC inline auto is_real(const t_matel x) {
#ifdef NRG_REAL
  return true;
#else
  return abs(x.imag()) <= check_real_TOLERANCE;
#endif
}

// Check if x is real and return the real part, i.e. x.real().
CONSTFNC inline auto check_real(const t_matel x) {
#ifdef NRG_REAL
  return x;
#else
  if (!is_real(x)) cout << "Warning: expected real number, but got " << x << endl;
  return x.real();
#endif
}

// Spur (real part)
CONSTFNC double trace_real(const Matrix &m) {
  assert_issquare(m);
  bucket sum;
  for (auto i = 0; i < m.size2(); i++) sum += check_real(m(i, i));
  return sum;
}

// As above, no check for finiteness.
CONSTFNC auto trace_real_nochecks(const Matrix &m) {
  assert_issquare(m);
  double sum = 0.0; // enforce real
  for (auto i = 0; i < m.size2(); i++) sum += check_real(m(i, i));
  return sum;
}

cmpl csqrt(const cmpl z) { return sqrt(z); }

#endif
