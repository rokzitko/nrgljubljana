#include "resample.hpp"

#include <common/version.hpp>

int main(int argc, char *argv[]) {
  if (NRG::Tools::report_version_if_requested(argc, argv, "resample")) return EXIT_SUCCESS;

  NRG::Resample::Resample<double> resample(argc, argv);
  resample.run();
}
