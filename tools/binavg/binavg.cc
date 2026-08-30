#include <exception>
#include <iostream>
#include <iomanip>
#include <ostream>
#include <cstdlib>

#include <common/version.hpp>

#include "binavg.hpp"

const int cout_PREC = 18; // Precision for verbose reporting on console

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "binavg")) return EXIT_SUCCESS;
  try {
    std::cout << "binavg - binned binary data averaging tool" << std::endl;
    std::cout << std::setprecision(cout_PREC);
    std::cerr << std::setprecision(cout_PREC);
    NRG::BinAvg::BinAvg binavg(argc, argv);
    binavg.calc();
  } catch (const std::exception &e) {
    std::cerr << "binavg: error: " << e.what() << std::endl;
    return EXIT_FAILURE;
  } catch (...) {
    std::cerr << "binavg: error: unknown exception" << std::endl;
    return EXIT_FAILURE;
  }
}
