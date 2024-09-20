// Computes \int rho(E)/(z-E) dE for given z and tabulated density of states rho(E).
//
// Legacy modes:
// Mode 1: <Re z> <Im z> as input. Returns Im part by default, or both Re and Im parts if the -G switch is used.
// Mode 2: read x and y from a file.
// Mode 3: convert imsigma/resigma.dat to imaw/reaw.dat files in the DMFT loop.
//
// Rok Zitko, rok.zitko@ijs.si, 2009-2020

#ifndef _hilb_hilb_hpp_
#define _hilb_hilb_hpp_

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <utility>
#include <functional>
#include <vector>
#include <string>
#include <map>
#include <optional>
#include <ctime>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <unistd.h>

#include <gsl/gsl_errno.h> // GNU scientific library
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

namespace NRG::Hilb {

inline auto atof(const std::string &s) { return ::atof(s.c_str()); }

// Unwrap a lambda expression and evaluate it at x
// https://martin-ueding.de/articles/cpp-lambda-into-gsl/index.html
inline auto unwrap(const double x, void *p) {
  auto fp = static_cast<std::function<double(double)> *>(p);
  return (*fp)(x);
}

// Wrap around GSL integration routines
class integrator {
  private:
  size_t limit;                              // size of workspace
  bool throw_on_error;                        // if true, integration error will trigger a hard error
  gsl_integration_workspace *work = nullptr; // work space
  gsl_function F;                            // GSL function struct for evaluation of the integrand
  public:
  void initF() noexcept {
    F.function = &unwrap; // the actual function (lambda expression) will be provided as an argument to operator()
    F.params   = nullptr;
  }
  integrator(size_t _limit = 1000, bool _throw_on_error = false) : limit{_limit}, throw_on_error{_throw_on_error} {
    work = gsl_integration_workspace_alloc(limit);
    initF();
  }
  integrator(const integrator &X) : limit{X.limit}, throw_on_error{X.throw_on_error} {
    // keep the same workspace
    initF();
  }
  integrator(integrator &&X) : limit{X.limit}, throw_on_error{X.throw_on_error} {
    work   = X.work; // steal workspace
    X.work = nullptr;
    initF();
  }
  integrator &operator=(const integrator &X) {
    if (this == &X) return *this;
    limit          = X.limit;
    throw_on_error = X.throw_on_error;
    // keep the same workspace
    initF();
    return *this;
  }
  integrator &operator=(integrator &&X) {
    if (this == &X) return *this;
    limit          = X.limit;
    throw_on_error = X.throw_on_error;
    work           = X.work; // steal workspace
    X.work         = nullptr;
    initF();
    return *this;
  }
  ~integrator() {
    if (work) { gsl_integration_workspace_free(work); }
  }

  /**
     * Integrate function f on [a:b].
     *
     * @param f Function to be integrated
     * @param a Lower integration range boundary
     * @param b Upper integration range boundary
     * @param epsabs numeric integration epsilon (absolute)
     * @param epsrel numeric integration epsilon (relative)
     */
  auto operator()(std::function<double(double)> f, const double a, const double b, const double epsabs = 1e-14, const double epsrel = 1e-10) {
    F.params = &f;
    double result, error;
    const auto status = gsl_integration_qag(&F, a, b, epsabs, epsrel, limit, GSL_INTEG_GAUSS15, work, &result, &error);
    if (status && std::abs(result) > epsabs && throw_on_error) throw std::runtime_error("qag error: " + std::to_string(status) + " -- " + gsl_strerror(status));
    return result;
  }
};

// Wrap around GSL interpolation routines
class interpolator {
  private:
  size_t len;                      // number of data points
  std::vector<double> X, Y;        // X and Y tables
  double Xmin, Xmax;               // boundary points
  double oob_value;                // out-of-boundary value
  gsl_interp_accel *acc = nullptr; // workspace
  gsl_spline *spline    = nullptr; // spline data
  public:
  interpolator(const std::vector<double> &_X, const std::vector<double> &_Y, const double _oob_value = 0.0) : X{_X}, Y{_Y}, oob_value{_oob_value} {
    assert(std::is_sorted(X.begin(), X.end()));
    assert(X.size() == Y.size());
    acc    = gsl_interp_accel_alloc();
    len    = X.size();
    spline = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init(spline, X.data(), Y.data(), len);
    Xmin = X.front();
    Xmax = X.back();
  }
  interpolator(const interpolator &I) : len{I.len}, X{I.X}, Y{I.Y}, Xmin{I.Xmin}, Xmax{I.Xmax}, oob_value{I.oob_value} {
    // keep the same accelerator workspace
    if (spline) { gsl_spline_free(spline); } // we need new spline data object
    spline = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init(spline, X.data(), Y.data(), len);
  }
  interpolator(interpolator &&I) : len{I.len}, X{I.X}, Y{I.Y}, Xmin{I.Xmin}, Xmax{I.Xmax}, oob_value{I.oob_value} {
    acc      = I.acc; // steal workspace
    I.acc    = nullptr;
    spline   = I.spline; // steal spline data
    I.spline = nullptr;
  }
  interpolator &operator=(const interpolator &I) {
    if (this == &I) return *this;
    len       = I.len;
    X         = I.X;
    Y         = I.Y;
    Xmin      = I.Xmin;
    Xmax      = I.Xmax;
    oob_value = I.oob_value;
    // keep the same accelerator workspace
    if (spline) { gsl_spline_free(spline); } // we need new spline data object
    spline = gsl_spline_alloc(gsl_interp_cspline, len);
    gsl_spline_init(spline, X.data(), Y.data(), len);
    return *this;
  }
  interpolator &operator=(interpolator &&I) {
    if (this == &I) return *this;
    len       = I.len;
    X         = std::move(I.X);
    Y         = std::move(I.Y);
    Xmin      = I.Xmin;
    Xmax      = I.Xmax;
    oob_value = I.oob_value;
    acc       = I.acc; // steal workspace
    I.acc     = nullptr;
    spline    = I.spline; // steal spline data
    I.spline  = nullptr;
    return *this;
  }
  ~interpolator() {
    if (spline) { gsl_spline_free(spline); }
    if (acc) { gsl_interp_accel_free(acc); }
  }
  auto operator()(const double x) { return (Xmin <= x && x <= Xmax ? gsl_spline_eval(spline, x, acc) : oob_value); }
};

// Square of x
inline auto sqr(const double x) { return x * x; }

// Result of Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}] (atg -> imQ).
inline auto imQ(const double x, const double y, const double B) { return atan((-B + x) / y) - atan((B + x) / y); }

// Result of Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}] (logs -> reQ).
inline auto reQ(const double x, const double y, const double B) { return (-log(sqr(B - x) + sqr(y)) + log(sqr(B + x) + sqr(y))) / 2.0; }

// Calculate the (half)bandwidth, i.e., the size B of the enclosing interval [-B:B].
inline auto bandwidth(const std::vector<double> &X) {
  assert(std::is_sorted(X.begin(), X.end()));
  const auto Xmin = X.front();
  const auto Xmax = X.back();
  return std::max(abs(Xmin), abs(Xmax));
}

/**
   * Calculate the Hilbert transform of a given spectral function at fixed complex value z. This is the low-level routine, called from
   * other interfaces. The general strategy is to perform a direct integration of the defining integral for cases where z has a sufficiently
   * large imaginary part, otherwise the singularity is subtracted out and the calculation of the transformed integrand is performed using
   * an integration-variable substitution to better handle small values.
   *
   * @param rhor Real part of the spectral function
   * @param rhoi Imaginary part of the spectral function
   * @param B Half-bandwidth, i.e., the support of the spectral function is [-B:B]
   * @param z The complex value for which to evaluate the Hilbert transform
   * @param lim_direct value of y=Im(z) above which rho(E)/(x+Iy-E) is directly integrated, and below which the singularity is removed
   */
template <typename FNCR, typename FNCI> auto hilbert_transform(FNCR rhor, FNCI rhoi, const double B, const std::complex<double> z, const double lim_direct = 1e-3) {
  // Initialize GSL and set up the interpolation
  gsl_set_error_handler_off();
  integrator integr;
  const auto x = real(z);
  const auto y = imag(z);
  // Low-level Hilbert-transform routines. calcA routine handles the case with removed singularity and
  // perform the integration after a change of variables. calcB routine directly evaluates the defining
  // integral of the Hilbert transform. Real and imaginary parts are determined in separate steps.
  auto calcA = [&integr, x, y, B](auto f3p, auto f3m, auto d) -> double {
    const double W1 = (x - B) / abs(y); // Rescaled integration limits. Only the absolute value of y matters here.
    const double W2 = (B + x) / abs(y);
    assert(W2 >= W1);
    // Determine the integration limits depending on the values of (x,y).
    double lim1down = 1.0, lim1up = -1.0, lim2down = 1.0, lim2up = -1.0;
    bool inside;
    if (W1 < 0 && W2 > 0) {        // x within the band
      const double ln1016 = -36.8; // \approx log(10^-16)
      lim1down            = ln1016;
      lim1up              = log(-W1);
      lim2down            = ln1016;
      lim2up              = log(W2);
      inside              = true;
    } else if (W1 > 0 && W2 > 0) { // x above the band
      lim2down = log(W1);
      lim2up   = log(W2);
      inside   = false;
    } else if (W1 < 0 && W2 < 0) { // x below the band
      lim1down = log(-W2);
      lim1up   = log(-W1);
      inside   = false;
    } else { // special case: boundary points
      inside = true;
    }
    const auto result1 = (lim1down < lim1up ? integr(f3p, lim1down, lim1up) : 0.0);
    const auto result2 = (lim2down < lim2up ? integr(f3m, lim2down, lim2up) : 0.0);
    const auto result3 = (inside ? d : 0.0);
    return result1 + result2 + result3;
  };

  auto calcB = [&integr, B](auto f0) -> double { return integr(f0, -B, B); }; // direct integration
  auto calc  = [y, lim_direct, calcA, calcB](auto f3p, auto f3m, auto d, auto f0) { return (abs(y) < lim_direct ? calcA(f3p, f3m, d) : calcB(f0)); };

  // Re part of rho(omega)/(z-omega)
  auto ref0 = [x, y, &rhor, &rhoi](double omega) -> double { return (rhor(omega) * (x - omega) + rhoi(omega) * y) / (sqr(y) + sqr(x - omega)); };

  // Im part of rho(omega)/(z-omega)
  auto imf0 = [x,y,&rhor,&rhoi](double omega) -> double { return (rhor(omega)*(-y) + rhoi(omega)*(x-omega) )/(sqr(y)+sqr(x-omega)); };

  // Re part of rho(omega)/(z-omega) with the singularity subtracted out.
  auto ref1 = [x, y, &rhor, &rhoi](double omega) -> double {
    return ((rhor(omega) - rhor(x)) * (x - omega) + (rhoi(omega) - rhoi(x)) * (y)) / (sqr(y) + sqr(x - omega));
  };
  auto ref2  = [x, y, ref1](double W) -> double { return abs(y) * ref1(abs(y) * W + x); };
  auto ref3p = [ref2](double r) -> double { return ref2(exp(r)) * exp(r); };
  auto ref3m = [ref2](double r) -> double { return ref2(-exp(r)) * exp(r); };
  auto red   = rhor(x) * reQ(x, y, B) - rhoi(x) * imQ(x, y, B);

  // Im part of rho(omega)/(z-omega) with the singularity subtracted out.
  auto imf1 = [x, y, &rhor, &rhoi](double omega) -> double {
    return ((rhor(omega) - rhor(x)) * (-y) + (rhoi(omega) - rhoi(x)) * (x - omega)) / (sqr(y) + sqr(x - omega));
  };
  auto imf2  = [x, y, imf1](double W) -> double { return abs(y) * imf1(abs(y) * W + x); };
  auto imf3p = [imf2](double r) -> double { return imf2(exp(r)) * exp(r); };
  auto imf3m = [imf2](double r) -> double { return imf2(-exp(r)) * exp(r); };
  auto imd   = rhor(x) * imQ(x, y, B) + rhoi(x) * reQ(x, y, B);

  return std::complex(calc(ref3p, ref3m, red, ref0), calc(imf3p, imf3m, imd, imf0));
}

  /*
   * @param Xpts Frequency mesh
   * @param Rpts Real part of the spectral function
   * @param Ipts Imaginary part of the spectral function
   */ 
template <typename T> auto hilbert_transform(const T &Xpts, const T &Rpts, const T &Ipts, const std::complex<double> z, const double lim_direct = 1e-3) {
  interpolator rhor(Xpts, Rpts);
  interpolator rhoi(Xpts, Ipts);
  const double B = bandwidth(Xpts);
  return hilbert_transform(rhor, rhoi, B, z, lim_direct);
}

class Hilb {
  private:
  double scale = 1.0;         // scale factor
  double B     = 1.0 / scale; // half-bandwidth
  double Xmin = -B;
  double Xmax = +B;
  const int OUTPUT_PREC = 16;     // digits of output precision
  bool verbose     = false;
  bool G           = false; // G(z). Reports real and imaginary part.
  std::vector<double> Xpts, Ypts, Ipts;
  bool tabulated = false; // Use tabulated DOS. If false, use rho_Bethe().

  auto hilbert(const double x, const double y) {
    auto Bethe_fnc = [this](const auto w) { return abs(w*scale) < 1.0 ? 2.0 / M_PI * scale * sqrt(1 - sqr(w * scale)) : 0.0; };
    auto zero_fnc = []([[maybe_unused]] const auto w) { return 0.0; };
    const auto z = std::complex(x,y);
    return tabulated ? hilbert_transform(Xpts, Ypts, Ipts, z) : hilbert_transform(Bethe_fnc, zero_fnc, B, z);
  }

  void do_one(const double x, const double y, std::ostream &OUT) {
    if (verbose)std::cout << "z=" << std::complex(x, y) << std::endl;
    const auto res = hilbert(x, y);
    if (!G)
      OUT << res.imag() << std::endl;
    else
      OUT << res << std::endl;
  }

  void do_stream(std::istream &F, std::ostream &OUT) {
    while (F.good()) {
      double label, x, y;
      F >> label >> x >> y;
      if (!F.fail()) {
        const auto res = hilbert(x, y);
        if (!G)
          OUT << label << " " << res.imag() << std::endl;
        else
          OUT << label << " " << res.real() << " " << res.imag() << std::endl;      
      }
    }
  }

  void do_hilb(std::istream &Fr, std::istream &Fi, std::ostream &Or, std::ostream &Oi) {
    while (Fr.good()) {
      double label1, label2, x, y;
      Fr >> label1 >> x;
      Fi >> label2 >> y;
      if (!Fr.fail() && !Fi.fail()) {
        if (abs(label1 - label2) > 1e-6) throw std::runtime_error("Frequency mismatch in do_hilb()");
        // Ensure ImSigma is negative
        const double CLIPPING = 1e-8;
        y                     = std::min(y, -CLIPPING);
        // The argument to calcre/im() is actually omega-Sigma(omega)
        x              = label1 - x;
        y              = -y;
        const auto res = hilbert(x, y);
        double resre   = res.real();
        double resim   = res.imag();
        if (!G) {
          // We include the -1/pi factor here
          resre /= -M_PI;
          resim /= -M_PI;
        } // Otherwise we are saving the Green's function G itself and no factor is required!
        Or << label1 << " " << resre << std::endl;
        Oi << label2 << " " << resim << std::endl;
      }
    }
  }

  auto safe_open_rd(const std::string &filename) {
    std::ifstream F(filename);
    if (!F) throw std::runtime_error("Error opening file " + filename + " for reading.");
    return F;
  }

  auto safe_open_wr(const std::string &filename) {
    std::ofstream F(filename);
    if (!F) throw std::runtime_error("Error opening file " + filename + " for writing.");
    F << std::setprecision(OUTPUT_PREC);
    return F;
  }

  void load_dos(const std::string &filename) {
    if (verbose) { std::cout << "Density of states filename: " << filename << std::endl; }
    auto F = safe_open_rd(filename);
    while (F) {
      double x, y;
      F >> x >> y;
      if (!F.fail()) {
        assert(std::isfinite(x) && std::isfinite(y));
        Xpts.push_back(x);
        Ypts.push_back(y);
        Ipts.push_back(0);
      }
    }
  }

  void report_dos() {
    Xmin = Xpts.front();
    Xmax = Xpts.back();
    assert(Xmin < Xmax);
    B = std::max(std::abs(Xmin), std::abs(Xmax));
    interpolator rho(Xpts, Ypts);
    integrator integr;
    const auto sum = integr(rho, -B, B);
    if (!std::isfinite(sum)) throw std::runtime_error("Error: Integral is not a finite number.");
    if (verbose)std::cout << "Sum=" << sum << std::endl;
  }

  void info() {
    if (tabulated)
      std::cout << "Xmin=" << Xmin << " Xmax=" << Xmax << std::endl;
    else
      std::cout << "Semicircular DOS. scale=" << scale << std::endl;
    std::cout << "B=" << B << std::endl;
  }

  void about() {
    std::cout << "# hilb -- Hilbert transformer for arbitrary density of states." << std::endl;
    std::cout << "# Rok Zitko, rok.zitko@ijs.si, 2009-2020" << std::endl;
  }

  void usage() {
    std::cout << "Usage (1): hilb [options] <x> <y>" << std::endl;
    std::cout << "Usage (2): hilb [options] <inputfile>" << std::endl;
    std::cout << "Usage (3): hilb [options] <resigma.dat> <imsigma.dat> <reaw.dat> <imaw.dat>" << std::endl;
    std::cout << std::endl;
    std::cout << "Options:" << std::endl;
    std::cout << "-h        show help" << std::endl;
    std::cout << "-d <dos>  Load the density of state data from file 'dos'" << std::endl;
    std::cout << "          If this option is not used, the Bethe lattice DOS is assumed." << std::endl;
    std::cout << "-v        Increase verbosity" << std::endl;
    std::cout << "-s        Rescale factor 'scale' for the DOS." << std::endl;
    std::cout << "-B        Half-bandwidth 'B' of the Bethe lattice DOS." << std::endl;
    std::cout << "          Use either -s or -B. Default is scale=B=1." << std::endl;
    std::cout << "-G        Compute the Green's function. hilb then returns Re[G(z)] Im[G(z)]" << std::endl;
  }

  void parse_param_run(int argc, char *argv[]) {
    std::optional<std::ofstream> OUTFILE;
    char c;
    while (c = getopt(argc, argv, "hGd:vVs:B:o:"), c != -1) {
      switch (c) {
        case 'h': usage(); exit(EXIT_SUCCESS);
        case 'G': G = true; break;
        case 'd': 
          tabulated = true;
          load_dos(optarg);
          report_dos();
          break;
        case 'v': verbose = true; break;
        case 's':
          scale = atof(optarg);
          B     = 1 / scale;
          Xmin  = -B;
          Xmax  = +B;
          if (verbose) { std::cout << "scale=" << scale << " B=" << B << std::endl; }
          break;
        case 'B':
          B     = atof(optarg);
          scale = 1 / B;
          Xmin  = -B;
          Xmax  = +B;
          if (verbose) { std::cout << "scale=" << scale << " B=" << B << std::endl; }
          break;
        case 'o':
          OUTFILE = safe_open_wr(optarg);
          if (verbose) { std::cout << "Output file: " << optarg << std::endl; }
          break;
        default: abort();
      }
    }
    const auto remaining = argc - optind; // arguments left
    const std::vector<std::string> args(argv+optind, argv+argc); // NOLINT
    // Usage case 1: real (x,y) pairs from an input file.
    if (remaining == 1) {
      about();
      if (verbose) info();
      auto F = safe_open_rd(args[0]);
      do_stream(F, OUTFILE ? OUTFILE.value() : std::cout);
      return;
    }
    // Usage case 2: real a single (x,y) pair from the command line.
    if (remaining == 2) {
      const auto x = atof(args[0]);
      const auto y = atof(args[1]);
      do_one(x, y, OUTFILE ? OUTFILE.value() : std::cout);
      return;
    }
    // Usage case 3: convert self-energy to a spectral function
    if (remaining == 4) {
      about();
      if (verbose) info();
      auto Frs = safe_open_rd(args[0]); // Re[Sigma]
      auto Fis = safe_open_rd(args[1]); // Im[Sigma]
      auto Fra = safe_open_wr(args[2]); // Re[G] or Re[Aw] = -1/pi Re[G]
      auto Fia = safe_open_wr(args[3]); // Im[G] or Im[Aw] = -1/pi Im[G]
      do_hilb(Frs, Fis, Fra, Fia);
      return;
    }
    about();
    usage();
  }

  public:
  Hilb(int argc, char *argv[]) {
    std::cout << std::setprecision(OUTPUT_PREC);
    gsl_set_error_handler_off();
    parse_param_run(argc, argv);
  }
};

} // namespace

#endif
