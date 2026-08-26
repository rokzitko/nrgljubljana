#include "hilb.hpp"

int main(int argc, char *argv[]) {
  try {
    NRG::Hilb::Hilb hilb(argc, argv);
  } catch (const std::exception &error) {
    std::cerr << "hilb: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }
}
