// Resample a tabulated function on a new X grid using smooth interpolation functions (GSL).
// Rok Zitko, rok.zitko@ijs.si, May 2014

// The input file must consist of a table of space-separated (energy,
// value) pairs. Gauss-Kronrod quadrature rules are used.

// CHANGE LOG
// 22.5.2014 - first version

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <iterator>
#include <ostream>
#include <stdexcept>
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
#include "../common/gsl_config.hpp"

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
        NRG::Tools::InterpolationMethod interpolation_method = NRG::Tools::InterpolationMethod::akima;
        bool extrapolate = false; // enable constant extrapolation with -e
        std::optional<T> extrapolation_below;
        std::optional<T> extrapolation_above;
        T Xmin{};
        T Xmax{};

        std::unique_ptr<gsl_interp_accel, gsl_acc_del> acc;
        std::unique_ptr<gsl_spline, gsl_spline_del> spline;

        std::vector<std::pair<T, T>> grid;

        void usage() 
        {
            std::cout << "\nUsage: resample [-h] [-v] [-i method] [-p precision] [-e [-a A] [-b B]] <input> <grid> <output>" << std::endl;
            std::cout << "-h, --help: show this help" << std::endl;
            std::cout << "-i, --interpolation METHOD: linear, cspline, akima, or steffen (default: akima)" << std::endl;
            std::cout << "-v: toggle verbose messages (now=" << verbose << ")" << std::endl;
            std::cout << "-e: enable constant extrapolation outside the input x range" << std::endl;
            std::cout << "-a A: use A for x < x_min (requires -e; default: y at x_min)" << std::endl;
            std::cout << "-b B: use B for x > x_max (requires -e; default: y at x_max)" << std::endl;
        }

        void parse_param(int argc, char *argv[]) 
        {
            const option long_options[] = {
                {"help", no_argument, nullptr, 'h'},
                {"interpolation", required_argument, nullptr, 'i'},
                {nullptr, 0, nullptr, 0}
            };
            int c;
            while ((c = getopt_long(argc, argv, "hvei:p:a:b:", long_options, nullptr)) != -1)
            {
                switch (c) 
                {
                case 'h':
                    usage();
                    std::exit(EXIT_SUCCESS);
                case 'i': interpolation_method = NRG::Tools::parse_interpolation_method(optarg); break;
                case 'v': verbose = true; break;
                case 'e': extrapolate = true; break;
                case 'p': output_precision = atoi(optarg); break;
                case 'a': extrapolation_below = static_cast<T>(std::stod(optarg)); break;
                case 'b': extrapolation_above = static_cast<T>(std::stod(optarg)); break;
                default: throw std::runtime_error("Unknown argument "s + std::string(1, static_cast<char>(c)));
                }
            }
            if (!extrapolate && (extrapolation_below || extrapolation_above))
              throw std::invalid_argument("Options -a and -b require -e.");
            int remaining = argc - optind;
            if (remaining != 3) {
                about();
                usage();
                std::exit(1);
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

        Resample(std::string inputfn_, std::string gridfn_, std::optional<std::string> outputfn_ = std::nullopt, bool verbose_ = false,
                 int output_precision_ = 16,
                 NRG::Tools::InterpolationMethod interpolation_method_ = NRG::Tools::InterpolationMethod::akima):
        inputfn(inputfn_), gridfn(gridfn_), outputfn(outputfn_), verbose(verbose_), output_precision(output_precision_),
        interpolation_method(interpolation_method_)
        {
            std::vector<std::pair<T, T>> f = readtable<T,T>(inputfn, verbose);
            grid = readtable<T,T>(gridfn, verbose);
            init(f);
        }

        Resample(std::vector<std::pair<T, T>> f, std::vector<std::pair<T, T>> grid_,
                 std::optional<std::string> outputfn_ = std::nullopt, bool verbose_ = false, int output_precision_ = 16,
                 NRG::Tools::InterpolationMethod interpolation_method_ = NRG::Tools::InterpolationMethod::akima):
        outputfn(outputfn_), verbose(verbose_), output_precision(output_precision_), interpolation_method(interpolation_method_), grid(grid_)
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
            if (im.empty()) throw std::runtime_error("No input data points available for resampling.");
            int len;      // number of data points
            std::vector<T>  Xpts, Ypts;

            std::sort(im.begin(), im.end());
            len  = im.size();
            const auto minimum_size = NRG::Tools::interpolation_minimum_size(interpolation_method);
            if (im.size() < minimum_size)
              throw std::runtime_error("Interpolation method " + std::string(NRG::Tools::interpolation_method_name(interpolation_method))
                                       + " requires at least " + std::to_string(minimum_size) + " input points.");
            Xmin = im.front().first;
            Xmax = im.back().first;
            if (!extrapolation_below) extrapolation_below = im.front().second;
            if (!extrapolation_above) extrapolation_above = im.back().second;
            if (verbose) std::cout << "Range: [" << Xmin << " ; " << Xmax << "]" << std::endl;

            std::transform(im.begin(), im.end(), std::back_inserter(Xpts), [] (const auto& pair){return pair.first;});
            std::transform(im.begin(), im.end(), std::back_inserter(Ypts), [] (const auto& pair){return pair.second;});

            gsl_set_error_handler_off();
            acc.reset(gsl_interp_accel_alloc());
            if (!acc) throw std::runtime_error("Failed to allocate GSL interpolation accelerator.");
            const gsl_interp_type *Interp_type = NRG::Tools::gsl_interpolation_type(interpolation_method);
            spline.reset(gsl_spline_alloc(Interp_type, len));
            if (!spline) throw std::runtime_error("Failed to allocate GSL spline.");
            if (const auto status = gsl_spline_init(spline.get(), Xpts.data(), Ypts.data(), len); status != 0)
              throw std::runtime_error("Failed to initialize GSL spline: "s + gsl_strerror(status));
        }

        void resample(std::vector<std::pair<T, T>> &grid_)
        {
            for (auto & i : grid_) {
              if (extrapolate && i.first < Xmin) i.second = *extrapolation_below;
              else if (extrapolate && i.first > Xmax) i.second = *extrapolation_above;
              else i.second = gsl_spline_eval(spline.get(), i.first, acc.get());
            }
        }
};

} //namespace
