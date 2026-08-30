#ifndef _h5write_h5write_hpp_
#define _h5write_h5write_hpp_

#include <cstdlib>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <traits.hpp>
#include <basicio.hpp>
#include <io.hpp>
#include <misc.hpp>
#include <h5.hpp>

#include "../common/diagnostics.hpp"

namespace NRG::H5Write {

using namespace NRG;

void h5_copy_matrix_from_file(H5Easy::File &file, const std::string &path, const std::string &input_filename) {
  const auto M = read_matrix_text(generate_matrix<double>, input_filename);
  h5_dump_matrix(file, path, M);
}

class H5Write {
 private:
    bool truncate = false;
    bool scalar = false;
    bool verbose = false;
    bool veryverbose = false;
    double scalar_value = 0.0;
   std::string filename_h5;
   std::string ds_path;
   std::string input;
   std::unique_ptr<H5Easy::File> h5;

   void usage() {
      std::cout << "\nUsage: h5write [-h] [-v|-vv] [-s] [-t] <hdf5 file> <path> <input_file|scalar>\n"
                << "  -v              print resolved configuration to stderr\n"
                << "  -vv             enable very verbose diagnostics\n"
                << "  -V, --version   print version and exit\n";
   }

   void parse_cmd_line(int argc, char *argv[]) {
     int c;
      while ((c = getopt(argc, argv, "hvst")) != -1) {
       switch (c) {
        case 'h': usage(); std::exit(EXIT_SUCCESS);
        case 'v':
          if (verbose) veryverbose = true;
          verbose = true;
          break;
       case 's': scalar = true; break;
       case 't': truncate = true; break;
       default: 
         usage();
         throw std::runtime_error("Invalid input switch.");
       }
     }
     const auto remaining = argc-optind;
     if (remaining != 3) {
       usage();
       throw std::runtime_error("Invalid number of parameters.");
     }
     const std::vector<std::string> args(argv+optind, argv+argc);
     filename_h5 = args[0];
      ds_path = args[1];
      input = args[2];
      if (scalar) scalar_value = atof(input);
    }

    void report_configuration() const {
      NRG::Tools::ConfigurationReport report("h5write");
      report.value("verbosity", veryverbose ? 2 : 1);
      report.value("hdf5_file", filename_h5);
      report.value("dataset_path", ds_path);
      report.value("dataset_mode", scalar ? "scalar" : "matrix");
      report.value("open_mode", truncate ? "truncate" : "read_write");
      if (scalar) {
        report.value("scalar_input", input);
        report.value("scalar_value", scalar_value);
      } else {
        report.value("matrix_input_file", input);
        report.value("matrix_input_mode", "text");
      }
      report.write(std::cerr);
    }
 public:
    H5Write(int argc, char *argv[]) {
      parse_cmd_line(argc, argv);
      if (verbose) report_configuration();
   }

   void run() {
     h5 = std::make_unique<H5Easy::File>(filename_h5, truncate ? H5Easy::File::Truncate : H5Easy::File::ReadWrite);
      if (scalar)
        h5_dump_scalar(*h5, ds_path, scalar_value);
     else
       h5_copy_matrix_from_file(*h5, ds_path, input);
   }
};

} // namespace

#endif
