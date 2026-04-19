#ifndef _tools_common_linint_hpp_
#define _tools_common_linint_hpp_

#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace NRG::Tools {

template<typename Vec, typename ErrorPolicy>
class LinIntBase {
 protected:
  Vec vec;
  int len{};
  int index{};
  double x0{}, x1{};
  double f0{}, f1{};
  double deriv{};
  bool newintegral_flag{};
  double xmin{}, xmax{};
  double fxmin{}, fxmax{};
  ErrorPolicy error_policy{};

 public:
  LinIntBase() = default;
  explicit LinIntBase(const Vec &in_vec) : vec(in_vec) {
    len = vec.size();
    if (len < 2)
      error_policy("At least two data points required for interpolation.");
    index            = -1;
    newintegral_flag = false;
    xmin             = vec.front().first;
    fxmin            = vec.front().second;
    xmax             = vec.back().first;
    fxmax            = vec.back().second;
  }

  void findindex(const double x) {
    if (index == -1) {
      for (int i = 0; i < len - 1; i++) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else if (x >= x1) {
      for (int i = index + 1; i < len - 1; i++) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else {
      for (int i = index - 1; i >= 0; i--) {
        x0 = vec[i].first;
        x1 = vec[i + 1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    }

    if (!(0 <= index && index < len && x0 <= x && x <= x1)) {
      std::ostringstream msg;
      msg << "findindex() error. x=" << x << " x0=" << x0 << " x1=" << x1 << " index=" << index << " len=" << len;
      error_policy(msg.str());
    }

    f0 = vec[index].second;
    f1 = vec[index + 1].second;
    const auto delta_y = f1 - f0;
    const auto delta_x = x1 - x0;
    deriv              = delta_y / delta_x;
  }

  double operator()(const double x) {
    if (x <= xmin) { return fxmin; }
    if (x >= xmax) { return fxmax; }
    if (index == -1 || !(x0 <= x && x < x1)) {
      newintegral_flag = true;
      findindex(x);
    }
    const auto dx = x - x0;
    return f0 + deriv * dx;
  }

  bool flag() const { return newintegral_flag; }
  void clear_flag() { newintegral_flag = false; }
};

template<typename Vec, typename ErrorPolicy>
class IntLinIntBase : public LinIntBase<Vec, ErrorPolicy> {
 protected:
  using Base = LinIntBase<Vec, ErrorPolicy>;
  Vec intvec;
  double intfxmin{}, intfxmax{};

 public:
  IntLinIntBase() = default;
  IntLinIntBase(const Vec &in_vec, const Vec &in_intvec) : Base(in_vec), intvec(in_intvec) {
    intfxmin = intvec.front().second;
    intfxmax = intvec.back().second;
  }

  double operator()(const double x) {
    if (x <= Base::xmin) { return intfxmin + Base::fxmin * (x - Base::xmin); }
    if (x >= Base::xmax) { return intfxmax + Base::fxmax * (x - Base::xmax); }
    if (Base::index == -1 || !(Base::x0 <= x && x < Base::x1)) { Base::findindex(x); }
    const auto dx = x - Base::x0;
    const auto int_from_0 = intvec[Base::index].second;
    return int_from_0 + dx * Base::f0 + dx * dx * Base::deriv / 2.0;
  }
};

} // namespace NRG::Tools

#endif
