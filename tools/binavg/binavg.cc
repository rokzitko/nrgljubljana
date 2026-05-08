#include <iostream>
#include <iomanip>
#include <string_view>
#include <cstdlib>
#include "binavg.hpp"

constexpr std::string_view VERSION = "0.0.2";
const int cout_PREC = 18; // Precision for verbose reporting on console

int main(int argc, char *argv[]) {
  try {
    std::cout << "binvag - binned binary data averaging tool - " << VERSION << std::endl;
    std::cout << "Rok Zitko, rok.zitko@ijs.si, 2013" << std::endl;
    std::cout << std::setprecision(cout_PREC);
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
