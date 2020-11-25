// Resample a tabulated function on a new X grid using smooth interpolation functions (GSL).
// Rok Zitko, rok.zitko@ijs.si, May 2014

// The input file must consist of a table of space-separated (energy,
// value) pairs. Gauss-Kronrod quadrature rules are used.

// CHANGE LOG
// 22.5.2014 - first version

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <utility>
#include <cassert>
#include <string>
#include <cstring>
#include <algorithm>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <getopt.h>

#include "misc.hpp"
#include "basicio.hpp"

using namespace NRG;

namespace NRG::Resample{

template <typename T>
class Resample
{
    private:
    
        // Dump additional information to stdout?
        bool verbose = false; // enable with -v
        std::string inputfn;  // Filename for input data
        std::string gridfn;   // Filename for a new X grid
        std::string outputfn; // Filename for resampled data. May be the same as gridfn.
        int output_precision = 16; // number of digits of precision in the output
        // const size_t limit = 1000;
        // const T EPSABS = 1e-12; // numeric integration epsilon (absolute)
        // const T EPSREL = 1e-8;  // numeric integration epsilon (relative)

        int len;           // number of data points
        T Xmin, Xmax; // the interval boundaries
        std::vector<T>  Xpts, Ypts;
        
        gsl_interp_accel *acc;
        gsl_spline *spline;
        // gsl_integration_workspace *w;



        void usage() 
        {
            std::cout << "\nUsage: resample [-h] [-v] <input> <grid> <output>" << std::endl;
            std::cout << "-v: toggle verbose messages (now=" << verbose << ")" << std::endl;
        }

        void parse_param(int argc, char *argv[]) 
        {
            char c;
            while (c = getopt(argc, argv, "hvp:"), c != -1) 
            {
                switch (c) 
                {
                case 'h':
                    usage();
                    exit(EXIT_SUCCESS);
                case 'v': verbose = true; break;
                case 'p': output_precision = atoi(optarg); break;
                default: throw std::runtime_error("Unknown argument " + c);
                }
            }
            int remaining = argc - optind;
            if (remaining != 3) {
                about();
                usage();
                exit(1);
            }
            inputfn  = std::string(argv[optind++]);
            gridfn   = std::string(argv[optind++]);
            outputfn = std::string(argv[optind++]);
        }

        void about() 
        {
            std::cout << "Resampling tool" << std::endl;
            #ifdef __TIMESTAMP__
            std::cout << "Timestamp: " << __TIMESTAMP__ << std::endl;
            #endif
            std::cout << "Compiled on " << __DATE__ << " at " << __TIME__ << std::endl;
        }
             

    public:
        Resample(int argv, char *argc[])
        {
            parse_param(argv, argc);
            if (verbose) about();
        }

        Resample(std::string inputfn, std::string gridfn, std::string outputfn, bool verbose = false, int output_precision = 16){
            this->inputfn = inputfn;
            this->gridfn = gridfn;
            this->outputfn = outputfn;
            this->verbose = verbose;
            this->output_precision = output_precision;
            if (verbose) about();
        }
        
        ~Resample() 
        {
            gsl_spline_free(spline);
            gsl_interp_accel_free(acc);
        }

        void run()
        {
            auto Fin = safe_open_for_reading(inputfn);
            std::vector<std::pair<double, double>> f;
            readtable(Fin, f, verbose);
            init(f);

            auto Fgrid = safe_open_for_reading(gridfn);
            std::vector<std::pair<double, double>> grid;
            readtable(Fgrid, grid, verbose);
            resample(grid);

            auto Fout = safe_open(outputfn);
            writetable(grid, Fout, output_precision);
        }

        void init(std::vector<std::pair<T, T>> &im) 
        {
            std::sort(im.begin(), im.end());
            len  = im.size();
            Xmin = im[0].first;
            Xmax = im[len - 1].first;
            if (verbose) std::cout << "Range: [" << Xmin << " ; " << Xmax << "]" << std::endl;
            // Xpts are increasing
            Xpts = std::vector<T> (len);
            Ypts = std::vector<T> (len);
            for (int i = 0; i < len; i++) 
            {
                Xpts[i] = im[i].first;
                Ypts[i] = im[i].second;
            }
            acc = gsl_interp_accel_alloc();
            const gsl_interp_type *Interp_type = gsl_interp_akima;
            spline = gsl_spline_alloc(Interp_type, len);
            gsl_spline_init(spline, &Xpts[0], &Ypts[0], len);
            // w = gsl_integration_workspace_alloc(limit);
            gsl_set_error_handler_off();
        }

        void resample(std::vector<std::pair<T, T>> &grid) 
        {
            for (auto & i : grid) i.second = gsl_spline_eval(spline, i.first, acc);
        }  


};

} //namespace