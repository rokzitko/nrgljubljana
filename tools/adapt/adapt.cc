// Discretization ODE solver for NRG
// Adaptable discretization mesh code
//
// Linear interpolation of rho(omega), integration using trapezoidal method. 4th-order Runge-Kutta ODE solver. Secant
// method for refinement of parameter A.

#include <iostream>
#include <fstream>
#include <utility>
#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cctype>

using namespace std::string_literals;

#include "adapt.hpp"
using namespace NRG::Adapt;

void about(std::ostream &F = std::cout) {
  F << "# Discretization ODE solver" << std::endl;
  F << "# Rok Zitko, rok.zitko@ijs.si, 2008-2020" << std::endl;
}

void help(int argc, char **argv, const std::string &help_message)
{
  std::vector<std::string> args(argv+1, argv+argc); // NOLINT
  if (args.size() >= 1 && args[0] == "-h") {
    std::cout << help_message << std::endl;
    exit(EXIT_SUCCESS);
  }
}

const auto usage = "Usage: adapt [-h] [P|N] [param_filename]"s;

std::pair<Sign, std::string> cmd_line(int argc, char *argv[]) {
  Sign sign = Sign::POS;
  if (argc >= 2) {
    const auto first_char = toupper(argv[1][0]);
    switch (first_char) {
    case 'P': sign = Sign::POS; break;
    case 'N': sign = Sign::NEG; break;
    default: std::cerr << usage << std::endl; exit(1);
    }
  }
  const std::string param_fn = argc == 3 ? std::string(argv[2]) : "param";
  std::cout << "# ++ " << (sign == Sign::POS ? "POSITIVE" : "NEGATIVE") << std::endl;
  return {sign, param_fn};
}

int main(int argc, char *argv[]) {
  const clock_t start_clock = clock();
  about();
  help(argc, argv, usage);
  const auto [sign, param_fn] = cmd_line(argc, argv);
  Params P(param_fn);
  Adapt calc(P, sign);
  calc.run();
  const clock_t end_clock = clock();
  std::cout << "# Elapsed " << double(end_clock - start_clock) / CLOCKS_PER_SEC << " s" << std::endl;
}
