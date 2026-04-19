#ifndef _tools_common_lambda_hpp_
#define _tools_common_lambda_hpp_

#include <cmath>

namespace NRG::Tools {

class LambdaCache {
 private:
  double lambda_ = -1.0;
  double log_lambda_{};
  double factor_lambda_{};

 public:
  LambdaCache() = default;
  explicit LambdaCache(const double lambda)
    : lambda_(lambda), log_lambda_(std::log(lambda)), factor_lambda_((1.0 - 1.0 / lambda) / std::log(lambda)) {}

  inline operator const double &() const { return lambda_; }
  inline double logL() const { return log_lambda_; }
  inline double factor() const { return factor_lambda_; }
  inline double power(const double x) const { return std::pow(lambda_, x); }
};

} // namespace NRG::Tools

#endif
