#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>

#include <common/version.hpp>

#include "hilb.hpp"

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "hilb")) return EXIT_SUCCESS;
  try {
    NRG::Hilb::Hilb hilb(argc, argv);
    return EXIT_SUCCESS;
  } catch (const std::exception &error) {
    std::cerr << "hilb: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
