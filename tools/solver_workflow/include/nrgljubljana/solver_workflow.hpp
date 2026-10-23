#ifndef NRGLJUBLJANA_SOLVER_WORKFLOW_HPP
#define NRGLJUBLJANA_SOLVER_WORKFLOW_HPP

#include <string>

namespace NRG::Workflow {

enum class Branch { positive, negative };

struct DiscretizationOptions {
  std::string calculation_directory;
  std::string parameter_file = "param";
  Branch branch = Branch::positive;
  bool integral = true;
};

// Existing calculation directory required. Relative input paths are resolved there.
// integral=true overrides f_method; false retains the parameter/CLI default.
void discretize(const DiscretizationOptions &);

struct InitializerOptions {
  std::string calculation_directory;
  std::string parameter_file = "param";
  std::string template_directory = "template";
  unsigned nmax = 0;
};

struct InitializerResult {
  std::string data_file;
  double ground_energy;
};

// nmax is required, already resolved by the caller, in [1,998]. The solver's
// parameter file is never modified. data_file is the absolute published path.
InitializerResult instantiate(const InitializerOptions &);

void solve(const std::string &calculation_directory, const std::string &absolute_scratch_directory);

// Calls are serialized, not generally thread-safe. The host must not overlap
// other solver/GSL operations, changes to GMP defaults or stream formatting, or
// any CWD-dependent I/O. CWD, GMP defaults and cout/cerr formatting are restored
// on success and exceptions; the caller's directory must remain accessible.
// Diagnostics may be emitted. See solver-workflow.md for supported modes.

} // namespace NRG::Workflow

#endif
