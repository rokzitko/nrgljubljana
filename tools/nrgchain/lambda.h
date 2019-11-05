// Discretization ODE solver for NRG
//
// ** LAMBDA class
//
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: lambda.h,v 1.1 2009/03/20 09:53:41 rok Exp $

// All things Lambda (stored to avoid recomputing).
class LAMBDA {
 private:
  double Lambda, logLambda, factorLambda;

 public:
  LAMBDA() { Lambda = -1; };
  LAMBDA(double in_Lambda) {
    Lambda = in_Lambda;
    logLambda = log(Lambda);
    factorLambda = (1.0-1.0/Lambda)/log(Lambda);
  }

  inline operator const double &() { return Lambda; }
  inline double logL() { return logLambda; }
  inline double factor() { return factorLambda; }
  inline double power(double x) { return pow(Lambda, x); }
};
