// Discretization ODE solver for NRG
// ** LAMBDA class

#include <cmath>

// All things Lambda (stored to avoid recomputing).
class LAMBDA {
 private:
  double Lambda, logLambda{}, factorLambda{};
 public:
  LAMBDA() { Lambda = -1; };
  LAMBDA(double in_Lambda) {
    Lambda       = in_Lambda;
    logLambda    = log(Lambda);
    factorLambda = (1.0 - 1.0 / Lambda) / log(Lambda);
  }
  inline operator const double &() const { return Lambda; }
  inline double logL() const { return logLambda; }
  inline double factor() const { return factorLambda; }
  inline double power(double x) const { return pow(Lambda, x); }
};
