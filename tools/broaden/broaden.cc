#include <iostream>
#include "broaden.hpp"

int main(int argc, char *argv[]) {
  std::cout << "broaden - finite-temperature broadening tool" << std::endl;
  NRG::Broaden::Broaden broaden(argc, argv);
  broaden.calc();
}
