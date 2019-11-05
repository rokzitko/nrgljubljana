// Computes \int rho(E)/(z-E) dE for given z and tabulated density
// of states rho(E). 
// 
// Mode 1: <Re z> <Im z> as imput. Returns Im part by default, or both Re and Im parts
// if the -G switch is used.
// Mode 2: read x and y from a file.
// Mode 3: convert imsigma/resigma.dat to imaw/reaw.dat files in the DMFT loop.
// 
// Rok Zitko, rok.zitko@ijs.si, 2009-2017

// CHANGE LOG
// 5.1.2017 - improved usage message
// 9.1.2017 - support for Im(z)<0
//          - integrand functions renamed for clarity
//          - support for calculating Re and Im parts.
// 22.11.2017 - support for imsigma/resigma -> imaw/reaw computation
// 13.12.2017 - Monotone cubic interpolation: gsl_interp_cspline -> gsl_interp_steffen

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <cassert>
#include <cfloat>
#include <utility>
#include <vector>
#include <string>
#include <map>
#include <ctime>
#include <unistd.h>
#include <complex>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_spline.h>

using namespace std;

typedef complex<double> cplx;

double x, y; // argument z=x+Iy. Global variable for simplicity.

double scale = 1.0; // scale factor
double B = 1.0/scale; // half-bandwidth

const size_t limit = 1000;
gsl_integration_workspace *w;

const double LIM_DIRECT = 1e-3; // value y where we switch over to direct
                                // integration of rho(E)/(z-E)
const double EPSABS = 1e-14; // numeric integration epsilon (absolute)
const double EPSREL = 1e-10; // (relative)
const double ln1016 = -36.8; // \approx log(10^-16)
const int OUTPUT_PREC = 16; // digits of output precision

bool verbose = false;
bool veryverbose = false;
bool G = false; // G(z). Reports real and imaginary part.

inline double sqr(double x) { return x*x; }

double rho_Bethe(double e)
{
   return  2.0/M_PI * scale * sqrt(1-sqr(e * scale));
}

typedef std::vector<double> DVEC;

DVEC Xpts, Ypts;
double Xmin, Xmax;

gsl_interp_accel *acc;
gsl_spline *spline;

bool tabulated = false; // Use tabulated DOS. If false, use rho_Bethe().

double rho(double e) 
{
   if (abs(e) >= B) {
      return 0.0;
   } else {
      if (tabulated) {
	 return gsl_spline_eval(spline, e, acc);
      } else {
	 return rho_Bethe(e);
      }      
   }   
}

// Im part of rho(omega)/(z-omega).
inline double imf0(double omega, void * params)
{
   return rho(omega) * (-y/(sqr(y)+sqr(x-omega)));
}

// As imf0, but with the singularity subtracted out.
// Called from imf2
inline double imf1(double omega)
{
   return (rho(omega)-rho(x)) * (-y/(sqr(y)+sqr(x-omega)));
}

// Called from imf3p, imf3m.
// W is the relative shift in units of abs(y).
inline double imf2(double W)
{
   return abs(y) * imf1(abs(y)*W+x);
}

inline double imf3p(double rho, void * params)
{
   return imf2(+exp(rho)) * exp(rho);
}

inline double imf3m(double rho, void * params)
{
   return imf2(-exp(rho)) * exp(rho);
}

inline double ref0(double omega, void * params)
{
   return rho(omega) * ((x-omega)/(sqr(y)+sqr(x-omega)));
}

inline double ref1(double omega)
{
   return (rho(omega)-rho(x)) * ((x-omega)/(sqr(y)+sqr(x-omega)));
}

inline double ref2(double W)
{
   return abs(y) * ref1(abs(y)*W+x);
}

inline double ref3p(double rho, void * params)
{
   return ref2(+exp(rho)) * exp(rho);
}

inline double ref3m(double rho, void * params)
{
   return ref2(-exp(rho)) * exp(rho);
}

// Integrate[(-y/(y^2 + (x - omega)^2)), {omega, -B, B}]
inline double atg(double x, double y)
{
   return atan((-B+x)/y)-atan((B+x)/y);
}

// Integrate[((x - omega)/(y^2 + (x - omega)^2)), {omega, -B, B}]
inline double logs(double x, double y)
{
   return (-log(sqr(B-x)+sqr(y))+log(sqr(B+x)+sqr(y)))/2.0;
}
   
void handle_qag(int status) 
{
   if (status && veryverbose) {
      cerr << "WARNING - qag error: " << status 
	<< " -- " << gsl_strerror(status) << endl;
   }
}

double calcim()
{
   gsl_function F;
   F.params = 0;
   
   if (abs(y) < LIM_DIRECT) {
      // Rescaled integration limits. Only the absolute value of y matters here.
      const double W1 = (x-B)/abs(y);
      const double W2 = (B+x)/abs(y);
      assert(W2 >= W1);

      // Default: no integration
      double lim1down = 1;
      double lim1up = -1;
      double lim2down = 1;
      double lim2up = -1;

      // x within the band
      if (W1 < 0 && W2 > 0) {
	 lim1down = ln1016;
	 lim1up = log(-W1);
	 lim2down = ln1016;
	 lim2up = log(W2);
      }
      
      // x above the band
      if (W1 > 0 && W2 > 0) {
	 lim2down = log(W1);
	 lim2up = log(W2);
      }
      
      // x below the band
      if (W1 < 0 && W2 < 0) {
	 lim1down = log(-W2);
	 lim1up = log(-W1);
      }
      
      if (veryverbose) {
	 cout << "W1=" << W1 << " W2=" << W2 << endl;
	 cout << "lim1 [" << lim1down << " " << lim1up << "]" << endl;
	 cout << "lim2 [" << lim2down << " " << lim2up << "]" << endl;
      }

      double result1, error1;
      F.function = &imf3p;

      if (lim1down < lim1up) {
	 int status = gsl_integration_qag(&F, lim1down, lim1up,
					  EPSABS, EPSREL, limit,
					  GSL_INTEG_GAUSS15, w, 
					  &result1, &error1);
	 handle_qag(status);
	 
      } else {
	 result1 = error1 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result1=" << result1 << " error1=" << error1 << endl;
      }

      double result2, error2;
      F.function = &imf3m;

      if (lim2down < lim2up) {
	 int status = gsl_integration_qag(&F, lim2down, lim2up,
					  EPSABS, EPSREL, limit,
					  GSL_INTEG_GAUSS15, w, 
					  &result2, &error2);
	 
	 handle_qag(status);
	 
      } else {
	 result2 = error2 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result2=" << result2 << " error2=" << error2 << endl;
      }
      
      // Contribution from the direct term
      double result3;
      
      if (W1 < 0 && W2 > 0) {
	 result3 = rho(x) * atg(x, y);
      } else {
	 result3 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result3=" << result3 << endl;
      }
      
      return result1 + result2 + result3;

   } else {
      // Case 2: abs(y) > LIM_DIRECT
      // In this case we may integrate the function directly, since it is
      // sufficiently smooth.
      
      double result, error;
      F.function = &imf0;
      
      int status = gsl_integration_qag(&F, -B, B,
				       EPSABS, EPSREL, limit,
				       GSL_INTEG_GAUSS15, w, 
				       &result, &error);
      
      handle_qag(status);
      
      if (veryverbose) {
	 cout << "Error=" << error << endl;
      }
      
      return result;
   }

   abort();
}

double calcre()
{
   gsl_function F;
   F.params = 0;

   if (abs(y) < LIM_DIRECT) {
      // Rescaled integration limits. Only the absolute value of y matters here.
      const double W1 = (x-B)/abs(y);
      const double W2 = (B+x)/abs(y);
      assert(W2 >= W1);

      // Default: no integration
      double lim1down = 1;
      double lim1up = -1;
      double lim2down = 1;
      double lim2up = -1;

      // x within the band
      if (W1 < 0 && W2 > 0) {
	 lim1down = ln1016;
	 lim1up = log(-W1);
	 lim2down = ln1016;
	 lim2up = log(W2);
      }
      
      // x above the band
      if (W1 > 0 && W2 > 0) {
	 lim2down = log(W1);
	 lim2up = log(W2);
      }
      
      // x below the band
      if (W1 < 0 && W2 < 0) {
	 lim1down = log(-W2);
	 lim1up = log(-W1);
      }
      
      if (veryverbose) {
	 cout << "W1=" << W1 << " W2=" << W2 << endl;
	 cout << "lim1 [" << lim1down << " " << lim1up << "]" << endl;
	 cout << "lim2 [" << lim2down << " " << lim2up << "]" << endl;
      }

      double result1, error1;
      F.function = &ref3p;

      if (lim1down < lim1up) {
	 int status = gsl_integration_qag(&F, lim1down, lim1up,
					  EPSABS, EPSREL, limit,
					  GSL_INTEG_GAUSS15, w, 
					  &result1, &error1);
	 handle_qag(status);
	 
      } else {
	 result1 = error1 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result1=" << result1 << " error1=" << error1 << endl;
      }

      double result2, error2;
      F.function = &ref3m;

      if (lim2down < lim2up) {
	 int status = gsl_integration_qag(&F, lim2down, lim2up,
					  EPSABS, EPSREL, limit,
					  GSL_INTEG_GAUSS15, w, 
					  &result2, &error2);
	 
	 handle_qag(status);
	 
      } else {
	 result2 = error2 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result2=" << result2 << " error2=" << error2 << endl;
      }
      
      // Contribution from the direct term
      double result3;
      
      if (W1 < 0 && W2 > 0) {
	 result3 = rho(x) * logs(x, y);
      } else {
	 result3 = 0.0;
      }
      
      if (veryverbose) {
	 cout << "result3=" << result3 << endl;
      }
      
      return result1 + result2 + result3;
   } else {
      // Case 2: abs(y) > LIM_DIRECT
      // In this case we may integrate the function directly, since it is
      // sufficiently smooth.
      
      double result, error;
      F.function = &ref0;
      
      int status = gsl_integration_qag(&F, -B, B,
				       EPSABS, EPSREL, limit,
				       GSL_INTEG_GAUSS15, w, 
				       &result, &error);
      
      handle_qag(status);
      
      if (veryverbose) {
	 cout << "Error=" << error << endl;
      }
      
      return result;
   }

   abort();
}


void about()
{
   cout << "# hilb -- Hilbert transformer for arbitrary density of states." << endl;
   cout << "# Rok Zitko, rok.zitko@ijs.si, 2009-2017" << endl;
}

void usage()
{
   cout << "Usage (1): hilb [options] <x> <y>" << endl;
   cout << "Usage (2): hilb [options] <inputfile>" << endl;
   cout << "Usage (3): hilb [options] <resigma.dat> <imsigma.dat> <reaw.dat> <imaw.dat>" << endl;
   cout << endl;
   cout << "Options:" << endl;
   cout << "-d <dos>  Load the density of state data from file 'dos'" << endl;
   cout << "          If this option is not used, the Bethe lattice DOS is assumed." << endl;
   cout << "-v        Increase verbosity" << endl;
   cout << "-V        Increase verbosity even further" << endl;
   cout << "-s        Rescale factor 'scale' for the DOS." << endl;
   cout << "-B        Half-bandwidth 'B' of the Bethe lattice DOS." << endl;
   cout << "          Use either -s or -B. Default is scale=B=1." << endl;
   cout << "-G        Compute the Green's function. hilb then returns Re[G(z)] Im[G(z)] (mode 1 and 2)" << endl;
}

// Note that this one uses C++ convention for streaming out complex numbers.
void do_one(ostream &OUT = cout)
{
   if (verbose) {
      const cplx z = cplx(x,y);
      cout << "z=" << z << endl;
   }
   
   double resim = calcim();
   
   double resre = 0;
   if (G) {
      resre = calcre();
   }
   
   if (!G) {
      OUT << resim << endl;
   } else {
      const cplx res = cplx(resre, resim);
      OUT << res << endl;
   }
}


void do_stream(istream &F, ostream &OUT = cout)
{
  while (F.good()) {
    double label;
    F >> label >> x >> y;
    if (!F.fail()) {
      double resim = calcim();

      double resre;
      if (G) {
	resre = calcre();
      }

      if (!G) {
	OUT << label << " " << resim << endl;
      } else {
	OUT << label << " " << resre << " " << resim << endl;
      }
    }
  }
}

void do_hilb(istream &Fr, istream &Fi, ostream &Or, ostream &Oi)
{
   while (Fr.good()) {
      double label1, label2;
      Fr >> label1 >> x;
      Fi >> label2 >> y;
      if (!Fr.fail() && !Fi.fail()) {
	if (abs(label1-label2) > 1e-6) {
	  cerr << "Mismatch in do_hilb(). Exiting." << endl;
	  exit(1);
	}

	// Ensure ImSigma is negative
	const double CLIPPING = 1e-8;
	y = min(y, -CLIPPING);

	// The argument to calcre/im() is actually omega-Sigma(omega)
	x = label1 - x;
	y = -y;

	double resre = calcre();
	double resim = calcim();
	 
	// We include the conventional -1/pi factor here
	resre /= -M_PI;
	resim /= -M_PI;
	
	Or << label1 << " " << resre << endl;
	Oi << label2 << " " << resim << endl;
      }
   }
}

void load_dos(char *dosfilename)
{
   if (verbose) {
      cout << "Density of states filename: " << dosfilename << endl;
   }
   
   ifstream F(dosfilename);
   if (!F) {
      cerr << "Can't open " << dosfilename << " for reading." << endl;
      exit(1);
   }
   
   while (F) {
      F >> x >> y;
      
      if (!F.fail()) {
	 assert(isfinite(x) && isfinite(y));
	 
	 Xpts.push_back(x);
	 Ypts.push_back(y);
      }
   }

   F.close();

   int len = Xpts.size();
   
   acc = gsl_interp_accel_alloc();
   spline = gsl_spline_alloc(gsl_interp_steffen, len);
   gsl_spline_init(spline, &Xpts[0], &Ypts[0], len);
   
   Xmin = Xpts[0];
   Xmax = Xpts[len-1];
   assert(Xmin < Xmax);
   B = max(Xmin, Xmax);
   
   double sum = gsl_spline_eval_integ(spline, Xmin, Xmax, acc);
   
   if (!isfinite(sum)) {
      cerr << "Error: Integral is not a finite number." << endl;
      exit(1);
   }
   
   if (verbose)
     cout << "Sum=" << sum << endl;
   
   tabulated = true;
}

// Open with error testing
void safeopen(ifstream &F, char *filename)
{
   F.open(filename);
   if (!F) {
      cerr << "Error opening file " << filename << " for reading." << endl;
      exit(1);
   }
}

void info()
{
   if (!verbose)
      return;
   
   if (tabulated) {
      cout << "Xmin=" << Xmin << " Xmax=" << Xmax << endl;
   } else {
      cout << "Semicircular DOS. scale=" << scale << endl;
   }
   cout << "B=" << B << endl;
}

void safeopen(ofstream &F, char *filename)
{
   F.open(filename);
   if (!F) {
      cerr << "Error opening file " << filename << " for writing." << endl;
      exit(1);
   }
}

void parse_param_run(int argc, char *argv[])
{
   char c;
   
   ostream * OUT = &cout;
   ofstream OUTFILE;

   while ((c = getopt(argc, argv, "Gd:vVs:B:o:")) != -1) {
      switch (c) {
       case 'G':
	 G = true;
	 break;
	 
       case 'd':
	 load_dos(optarg);
	 break;
	 
       case 'v':
	 verbose = true;
	 break;
	 
       case 'V':
	 veryverbose = true;
	 break;
	 
       case 's':
	 scale = atof(optarg);
	 B = 1/scale;
	 if (verbose) {
	    cout << "scale=" << scale << " B=" << B << endl;
	 }
	 break;
	 
       case 'B':
	 B = atof(optarg);
	 scale = 1/B;
	 if (verbose) {
	    cout << "scale=" << scale << " B=" << B << endl;
	 }
	 break;
	 
       case 'o':
	   {
	      OUTFILE.open(optarg);
	      if (!OUTFILE) {
		 cerr << "Can't open " << optarg << " for writing." << endl;
		 exit(1);
	      }
	      OUT = &OUTFILE;
	      OUTFILE << setprecision(OUTPUT_PREC);
	      if (verbose) {
		 cout << "Output file: " << optarg << endl;
	      }
	   }
	 break;
	 
       default:
	 abort();
      }
   }
   
   int remaining = argc-optind; // arguments left

   // Usage case 1: real (x,y) pairs from an input file.
   if (remaining == 1) {
      about();
      info();
      
      char *filename = argv[optind];
      
      ifstream F;
      safeopen(F, filename);
      do_stream(F, *OUT);
      F.close();
      return;
   }
   
   // Usage case 2: real a single (x,y) pair from the command line.
   if (remaining == 2) {
      x = atof(argv[optind]);
      y = atof(argv[optind+1]);
      
      do_one(*OUT);
      return;
   }
   
   if (remaining == 4) {
      about();
      info();
      
      char *fnresigma = argv[optind];
      char *fnimsigma = argv[optind+1];
      char *fnreaw = argv[optind+2];
      char *fnimaw = argv[optind+3];
      
      ifstream Frs, Fis;
      ofstream Fra, Fia;
      safeopen(Frs, fnresigma);
      safeopen(Fis, fnimsigma);
      safeopen(Fra, fnreaw);
      safeopen(Fia, fnimaw);
      Fra << setprecision(OUTPUT_PREC);
      Fia << setprecision(OUTPUT_PREC);
      
      do_hilb(Frs, Fis, Fra, Fia);
      return;
   }
   
   about();
   usage();
}

void init()
{
   cout << setprecision(OUTPUT_PREC);
   gsl_set_error_handler_off();
   w = gsl_integration_workspace_alloc(limit);
}

void done()
{
   gsl_integration_workspace_free(w);
}

int main(int argc, char *argv[])
{
   init();
   parse_param_run(argc, argv);
   done();
}
