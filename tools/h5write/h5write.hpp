#ifndef _h5write_h5write_hpp_
#define _h5write_h5write_hpp_

#include <iostream>

#include <traits.hpp>
#include <basicio.hpp>
#include <misc.hpp>
#include <h5.hpp>

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
   std::string filename_h5;
   std::string ds_path;
   std::string input;
   std::unique_ptr<H5Easy::File> h5;

   void usage() {
     std::cout << "\nUsage: h5write [-h] [-s] [-t] <hdf5 file> <path> <input_file|scalar>\n";
   }

   void parse_cmd_line(int argc, char *argv[]) {
     char c;
     while (c = getopt(argc, argv, "hst"), c != -1) {
       switch (c) {
       case 'h': usage(); exit(EXIT_SUCCESS);
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
   }
 public:
   H5Write(int argc, char *argv[]) {
     parse_cmd_line(argc, argv);
   }

   void run() {
     h5 = std::make_unique<H5Easy::File>(filename_h5, truncate ? H5Easy::File::Truncate : H5Easy::File::ReadWrite);
     if (scalar)
       h5_dump_scalar(*h5, ds_path, atof(input));
     else
       h5_copy_matrix_from_file(*h5, ds_path, input);
   }
};

} // namespace

#endif

