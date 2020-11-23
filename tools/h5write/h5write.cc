#include "h5write.hpp"

int main(int argc, char *argv[]) {
  NRG::H5Write::H5Write h5write(argc, argv);
  h5write.run();
}
