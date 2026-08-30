#include <common/version.hpp>

#include "h5write.hpp"

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "h5write")) return EXIT_SUCCESS;
  NRG::H5Write::H5Write h5write(argc, argv);
  h5write.run();
}
