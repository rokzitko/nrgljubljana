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
       
template <typename T>
inline std::string formatted_output(const T x, const Params &P) {
  return fmt::format("{x:>{width}}", "x"_a=x, "width"_a=P.width_custom);
}

inline std::string formatted_output(const double x, const Params &P) {
  return fmt::format("{x:>{width}.{prec}}", "x"_a=x, "prec"_a=P.prec_custom, "width"_a=P.width_custom);
}

template <typename T>
inline bool negligible_imag_part(const std::complex<T> &z, const double output_imag_eps = 1e-13) {
  return abs(z.imag()) < abs(z.real()) * output_imag_eps;
}


// The output format for complex values is X+IY or X-IY, where X and Y are real and imaginary part, respectively. The
// imaginary part is only shown where its value relative to the real part is sufficiently large. No space is used in
// the outputted string in order to simplify parsing.
inline std::string formatted_output(const cmpl z, const Params &P) {
  const auto [r, i] = reim(z);
  const auto str = P.noimag || negligible_imag_part(z) ?
    fmt::format("{r:.{prec}f}", "r"_a=r, "prec"_a=P.prec_custom) :
    fmt::format("{r:.{prec}f}{s}I{absi:.{prec}f}", "r"_a=r, "s"_a=(i>0 ? "+" : "-"), "absi"_a=abs(i), "prec"_a=P.prec_custom);
  return fmt::format("{str:>{width}}", "str"_a=str, "width"_a=P.width_custom); // the width for the whole X+iY string
}

#endif
