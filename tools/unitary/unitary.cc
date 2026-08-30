#include "unitary.hpp"

#include <common/version.hpp>

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "unitary")) return EXIT_SUCCESS;
  using namespace NRG::Unitary;
  Unitary unit(argc, argv);
  unit.run(argc, argv);
}
