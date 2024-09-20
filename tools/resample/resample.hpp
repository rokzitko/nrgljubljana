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
using namespace std::string_literals;
#include <cstring>
#include <algorithm>
#include <optional>
#include <memory>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>

#include <unistd.h>
#include <getopt.h>

#include "misc.hpp"
#include "basicio.hpp"

using namespace NRG;

namespace NRG::Resample{

struct gsl_acc_del{
    void operator()(gsl_interp_accel* acc) {
        if (acc) gsl_interp_accel_free(acc);
    }
};

struct gsl_spline_del{
    void operator()(gsl_spline* spline){
        if(spline) gsl_spline_free(spline);
    }
};

template <typename T>
class Resample
{
    private:

        // Dump additional information to stdout?
        std::string inputfn;  // Filename for input data
        std::string gridfn;   // Filename for a new X grid
        std::optional<std::string> outputfn; // Filename for resampled data. May be the same as gridfn.
        bool verbose = false; // enable with -v
        int output_precision = 16; // number of digits of precision in the output

        std::unique_ptr<gsl_interp_accel, gsl_acc_del> acc;
        std::unique_ptr<gsl_spline, gsl_spline_del> spline;

        std::vector<std::pair<T, T>> grid;

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
                default: throw std::runtime_error("Unknown argument "s + c);
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
            std::vector<std::pair<T, T>> f = readtable<T,T>(inputfn, verbose);
            grid = readtable<T,T>(gridfn, verbose);
            init(f);
        }

        Resample(std::string inputfn, std::string gridfn, std::optional<std::string> outputfn = std::nullopt, bool verbose = false, int output_precision = 16):
        inputfn(inputfn), gridfn(gridfn), outputfn(outputfn), verbose(verbose), output_precision(output_precision)
        {
            std::vector<std::pair<T, T>> f = readtable<T,T>(inputfn, verbose);
            grid = readtable<T,T>(gridfn, verbose);
            init(f);
        }

        Resample(std::vector<std::pair<T, T>> f, std::vector<std::pair<T, T>> grid, std::optional<std::string> outputfn = std::nullopt, bool verbose = false, int output_precision = 16):
        grid(grid), outputfn(outputfn), verbose(verbose), output_precision(output_precision)
        {
            init(f);
        }

        std::optional<std::vector<std::pair<T, T>>> run()
        {
            resample(grid);

            if (outputfn)
            {
                writetable(grid, *outputfn, output_precision);
                return std::nullopt;
            }
            else return grid;
        }

        void init(std::vector<std::pair<T, T>> &im)
        {
            int len;      // number of data points
            T Xmin, Xmax; // the interval boundaries
            std::vector<T>  Xpts, Ypts;

            std::sort(im.begin(), im.end());
            len  = im.size();
            Xmin = im.front().first;
            Xmax = im.back().first;
            if (verbose) std::cout << "Range: [" << Xmin << " ; " << Xmax << "]" << std::endl;

            std::transform(im.begin(), im.end(), std::back_inserter(Xpts), [] (const auto& pair){return pair.first;});
            std::transform(im.begin(), im.end(), std::back_inserter(Ypts), [] (const auto& pair){return pair.second;});

            acc.reset(gsl_interp_accel_alloc());
            const gsl_interp_type *Interp_type = gsl_interp_akima;
            spline.reset(gsl_spline_alloc(Interp_type, len));
            gsl_spline_init(spline.get(), Xpts.data(), Ypts.data(), len);
            gsl_set_error_handler_off();
        }

        void resample(std::vector<std::pair<T, T>> &grid)
        {
            for (auto & i : grid) i.second = gsl_spline_eval(spline.get(), i.first, acc.get());
        }
};

} //namespace
