#include <exception>
#include <iostream>
#include <ostream>

#include <common/version.hpp>

#include "kk.hpp"

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "kk")) return EXIT_SUCCESS;
  try {
    NRG::KK::KK kk(argc, argv);
    return EXIT_SUCCESS;
  } catch (const std::exception &e) {
    std::cerr << "kk: error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "kk: error: unknown exception" << std::endl;
    return 1;
  }
}
