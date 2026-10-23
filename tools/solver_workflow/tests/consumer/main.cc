#include <nrgljubljana/solver_workflow.hpp>
#include <cmath>
#include <iostream>
#include <stdexcept>

int main(int argc, char **argv) {
  try {
    if (argc != 3 && argc != 5)
      throw std::invalid_argument("Expected calculation and absolute scratch directories [--prepare-only nmax]");
    const bool prepare_only = argc == 5 && std::string(argv[3]) == "--prepare-only";
    if (argc == 5 && !prepare_only) throw std::invalid_argument("Expected --prepare-only nmax");
    const auto nmax = prepare_only ? static_cast<unsigned>(std::stoul(argv[4])) : 6U;
    using namespace NRG::Workflow;
    discretize({argv[1]});
    discretize({argv[1], "param", Branch::negative});
    const auto result = instantiate({argv[1], "param", "template", nmax});
    if (!std::isfinite(result.ground_energy)) throw std::runtime_error("Invalid ground energy");
    if (!prepare_only) solve(argv[1], argv[2]);
    std::cout << result.data_file << '\n';
    return 0;
  } catch (const std::exception &error) {
    std::cerr << error.what() << '\n';
    return 1;
  }
}
