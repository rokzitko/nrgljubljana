#include "kk.hpp"

int main(int argc, char *argv[]) {
  try {
    NRG::KK::KK kk(argc, argv);
  } catch (const std::exception &e) {
    std::cerr << "kk: error: " << e.what() << std::endl;
    return 1;
  } catch (...) {
    std::cerr << "kk: error: unknown exception" << std::endl;
    return 1;
  }
}
