// Discretization ODE solver for NRG
//
// ** Interpolation code
// 
// Rok Zitko, zitko@theorie.physik.uni-goettingen.de, Dec 2008
// $Id: linint.h,v 1.1 2009/03/20 09:53:41 rok Exp $

// Structures for storing tabulated data, such as rho(omega).
typedef pair<double, double> Pair;
typedef vector<Pair> Vec;

// Linear interpolation class
class LinInt
{
 protected:
  Vec vec; // tabulated data
  int len; // length of vec
  int index; // index of the interval where last x was found
  double x0, x1; // last x was in [x0:x1]
  double f0, f1; // f(x0), f(x1)
  double deriv; // (f1-f0)/(x1-x0)
  bool newintegral_flag; // set to true when we switch to a new interval
  double xmin, xmax; // lowest and highest x contained in vec
  double fxmin, fxmax; // f(xmin), f(xmax)

 public:
  LinInt() {};
  LinInt(Vec &in_vec) : vec(in_vec) {
    len = vec.size();
    if (len < 2) {
       cerr << "At least two data points required for interpolation." << endl;
       exit(1);
    }     
    index = -1;
    newintegral_flag = false;
    xmin = vec.front().first;
    fxmin = vec.front().second;
    xmax = vec.back().first;
    fxmax = vec.back().second;
  };
  void findindex(double x);
  double operator()(double x);
  bool flag() { return newintegral_flag; }
  void clear_flag() { newintegral_flag = false; }
};

// Integral of a linear interpolation class
class IntLinInt : public LinInt
{
 protected:
  Vec intvec;
  double intfxmin, intfxmax; // int f(x) @ xmin and @ xmax

 public:
  IntLinInt() {};
  IntLinInt(Vec &in_vec, Vec &in_intvec) :
    LinInt(in_vec), intvec(in_intvec) {
    intfxmin = intvec.front().second;
    intfxmax = intvec.back().second;
  };
  double operator()(double x);
};

// Serach for the interval so that x is contained in [x0:x1].
void LinInt::findindex(double x)
{
  if (index == -1) {
    // When interpolation class object is constructed, index is
    // initialized to -1. We have no prior knowledge, thus we search in
    // the full interval.
    for (int i = 0; i < len-1; i++) {
      x0 = vec[i].first;
      x1 = vec[i+1].first;
      if (x0 <= x && x <= x1) {
        index = i;
        break;
      }
    }
  } else {
    if (x >= x1) {
      for (int i = index+1; i < len-1; i++) {
        x0 = vec[i].first;
        x1 = vec[i+1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    } else {
      for (int i = index-1; i >= 0; i--) {
        x0 = vec[i].first;
        x1 = vec[i+1].first;
        if (x0 <= x && x <= x1) {
          index = i;
          break;
        }
      }
    }
  }

  if (!(0 <= index && index < len && x0 <= x && x <= x1)) {
      cerr << "findindex() error."
           << " x=" << x
           << " x0=" << x0
           << " x1=" << x1
           << " index=" << index
           << " len=" << len << endl;
      exit(1);
  }

  f0 = vec[index].second; // f(x0)
  f1 = vec[index+1].second; // f(x1)
  double Delta_y = f1-f0;
  double Delta_x = x1-x0;
  deriv = Delta_y/Delta_x;
}

// Return y(x) using linear interpolation between the tabulated values.
double LinInt::operator()(double x)
{
  // Extrapolate if necessary
  if (x <= xmin) {
    return fxmin;
  }
  if (x >= xmax) {
    return fxmax;
  }

  if (index == -1 || !(x0 <= x && x < x1)) {
    newintegral_flag = true;
    findindex(x);
  }

  double dx = x-x0;
  return f0 + deriv * dx;
}

// Return int(y(x),x) using quadratic formula in the integral [x0:x1]. The
// results of the integration is thus consistent with linear interpolation
// between the tabulated data.
double IntLinInt::operator()(double x)
{
  // Extrapolate if necessary. Consistent with constant function
  // extrapolations in LinInt.
  if (x <= xmin) {
    return intfxmin + fxmin * (x-xmin);
  }
  if (x >= xmax) {
    return intfxmax + fxmax * (x-xmax);
  }

  // As opposed to LinInt, here we don't set newintegral_flag! The
  // integral may be considered smooth as compared to rho itself, thus the
  // advantage of breaking the integration intervals into piecewise 2nd
  // order polynomials for the int[rho] brings little additional benefit.
  if (index == -1 || !(x0 <= x && x < x1)) {
    findindex(x);
  }

  double dx = x-x0;
  double int_from_0 = intvec[index].second;
  return int_from_0 + dx * f0 + dx*dx * deriv/2.0;
}
