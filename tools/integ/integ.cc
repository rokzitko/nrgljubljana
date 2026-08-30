#include <cstdlib>
#include <exception>
#include <iostream>

#include <common/version.hpp>

#include "integ.hpp"

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "integ")) return EXIT_SUCCESS;
  try {
    NRG::Integ::run(argc, argv);
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "integ: error: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
