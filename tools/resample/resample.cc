#include "resample.hpp"

int main(int argv, char *argc[]) {

  NRG::Resample::Resample<double> resample(argv, argc);
  resample.run();
}