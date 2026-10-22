#include <nrgljubljana/numerics.hpp>

#include <cmath>

int main() {
  const NRG::Numerics::Table dos{{-1.0, 0.0, 1.0}, {0.5, 0.5, 0.5}};
  const NRG::Numerics::Table sigma{{0.0}, {0.0}};
  const auto result = NRG::Numerics::lattice(dos, sigma, sigma, 0.0, 0.1);
  return std::abs(NRG::Numerics::integrate(dos) - 1.0) < 1e-14
           && std::abs(result.real.y[0]) < 1e-14 && result.spectral.y[0] > 0.0 ? 0 : 1;
}
