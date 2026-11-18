// Solve calculations through the workflow from a host that initializes MPI itself. With one rank, every calculation
// is solved in turn; with one rank per calculation, the ranks solve them at the same time.
//
// Usage: workflow_solve_mpi <calculation directory>... <absolute scratch directory>

#include <nrgljubljana/solver_workflow.hpp>

#include <mpi.h>

#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);
  int rank = 0;
  int size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  int status = 0;
  try {
    if (argc < 3) throw std::invalid_argument("Usage: workflow_solve_mpi <calculation directory>... <absolute scratch directory>");
    const std::vector<std::string> calculations(argv + 1, argv + argc - 1);
    const std::string scratch = argv[argc - 1];
    if (size != 1 && static_cast<std::size_t>(size) != calculations.size())
      throw std::invalid_argument("Run with one rank or one rank per calculation");
    for (std::size_t i = 0; i < calculations.size(); ++i)
      if (size == 1 || static_cast<int>(i) == rank) NRG::Workflow::solve(calculations[i], scratch, "solve.log");
  } catch (const std::exception &error) {
    std::cerr << "Rank " << rank << ": " << error.what() << '\n';
    status = 1;
  }
  MPI_Finalize();
  return status;
}
