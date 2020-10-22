#ifndef _io_h_
#define _io_h_

#define HIGHPREC(val) std::setprecision(std::numeric_limits<double>::max_digits10) << (val)

inline std::string to_string(const std::complex<double> &z) {
  std::ostringstream s;
  s << z;
  return s.str();
}

template <typename T>
std::ostream & operator<<(std::ostream &os, const std::set<T> &x) {
  std::copy(x.cbegin(), x.cend(), std::ostream_iterator<T>(os, " "));
  return os;
}

// Returns a string with a floating value in fixed (non-exponential) format with N digits of precision after the
// decimal point.
inline std::string prec(const double x, const int N)
{
  std::ostringstream s;
  s << std::fixed << std::setprecision(N) << x;
  return s.str();
}
inline std::string prec3(const double x) { return prec(x, 3); }
       
#endif
