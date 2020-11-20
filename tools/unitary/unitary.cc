#include "unitary.hpp"

int main(int argc, char *argv[]) {
  using namespace NRG::Unitary;
  Unitary unit(argc, argv);
  unit.run(argc, argv);
}
