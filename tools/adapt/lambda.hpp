// Discretization ODE solver for NRG
// ** LAMBDA class

#ifndef _adapt_lambda_hpp_
#define _adapt_lambda_hpp_

#include <cmath>

namespace NRG::Adapt {

// All things Lambda (stored to avoid recomputing).
class LAMBDA {
 private:
   double Lambda, logLambda {}, factorLambda {};
 public:
   LAMBDA() : Lambda(-1) {}
   explicit LAMBDA(const double l) : Lambda(l), logLambda(log(l)), factorLambda((1.0-1.0/l)/log(l)) {}
   inline operator const double &() const { return Lambda; }
   inline double logL() const { return logLambda; }
   inline double factor() const { return factorLambda; }
   inline double power(double x) const { return pow(Lambda, x); }
};

} // namespace

#endif
