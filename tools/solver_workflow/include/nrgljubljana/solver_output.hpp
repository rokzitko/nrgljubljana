#ifndef NRGLJUBLJANA_SOLVER_OUTPUT_HPP
#define NRGLJUBLJANA_SOLVER_OUTPUT_HPP

#include <map>
#include <string>
#include <vector>

namespace NRG::Workflow {

// Readers for the text the solver writes, without I/O: the caller passes the contents of the file. Malformed text
// throws std::runtime_error naming the line. See output-formats.md.

// custom, td, customfdm or tdfdm: the column names of the last '#' line before the data (custom has a line of column
// numbers before it), and one row of values per further line.
struct StepTable {
  std::vector<std::string> columns;
  std::vector<std::vector<double>> rows; // one per NRG step (custom, td) or the final values (customfdm, tdfdm)
};
StepTable read_step_table(const std::string &text);

// customfdm or tdfdm: a table with exactly one row, as column name -> value
std::map<std::string, double> read_final_values(const std::string &text);

// A group of degenerate levels of one NRG step in annotated.dat
struct EnergyLevel {
  double energy = 0;
  std::vector<std::vector<int>> subspaces; // the quantum numbers of each invariant subspace in the group
  long degeneracy = 0;                     // the number of states, multiplicity included
  bool operator==(const EnergyLevel &) const = default;
};

// annotated.dat as written with dumpgroups=true (the default): per line "energy (qn ...) (qn ...) [degeneracy]", an
// empty line after each NRG step. The levels of each step; a step without levels is empty. A file written with
// dumpgroups=false (one level per line, no degeneracy) is refused.
std::vector<std::vector<EnergyLevel>> read_energy_flow(const std::string &text);

// The file of a spectrum without its extension (.dat broadened, .bin the bins with savebins=true), as the solver
// names it: list is the parameter that requests it (specd, specs, spect, specot, specq, specgt, speci1t, speci2t,
// specchit), correlator is "op1-op2", algorithm is FT, DMNRG, CFS, FDM, FDMmats, ..., and spin is +1 or -1 for the
// spectra of a doublet list (specd, specgt, speci1t, speci2t) under a symmetry type with a field (TemplateInfo::field),
// otherwise 0. Example: ("specd", "A_d-A_d", "FDM", 1) -> "spec_FDM_dens_A_d-A_d-u".
std::string spectrum_file_name(const std::string &list, const std::string &correlator, const std::string &algorithm,
                               int spin = 0);

} // namespace NRG::Workflow

#endif
