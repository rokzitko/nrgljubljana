#include <iostream>
#include "broaden.hpp"

int main(int argc, char *argv[]) {
  try {
    std::cout << "broaden - finite-temperature broadening tool" << std::endl;
    NRG::Broaden::Broaden broaden(argc, argv);
    broaden.calc();
  } catch (const std::exception &e) {
    std::cerr << "broaden: error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "broaden: error: unknown exception" << std::endl;
    return 1;
  }
}
