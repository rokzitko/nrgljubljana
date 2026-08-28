#include <cstdlib>
#include <exception>
#include <iostream>
#include <ostream>

#include "hilb.hpp"

int main(int argc, char *argv[]) {
  try {
    NRG::Hilb::Hilb hilb(argc, argv);
  } catch (const std::exception &error) {
    std::cerr << "hilb: " << error.what() << '\n';
    return EXIT_FAILURE;
  }
}
