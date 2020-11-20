#include "unitary.hpp"

int main(int argc, char *argv[]) {
  using namespace NRG::Unitary;
  parse_param(argc, argv);
  about();
  run(argc, argv);
}
