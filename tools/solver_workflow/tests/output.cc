// The readers of solver_output.hpp on the stored outputs of the regression tests.
//
// Usage: workflow_output <test/c++ directory>
// Every ref/ directory there: td, custom, tdfdm and customfdm read as tables, the last two also as final values;
// annotated.dat read as an energy flow with one step per empty line; every spectrum file name rebuilt by
// spectrum_file_name. Then the values of mpitest1_qs, and malformed text.

#include <nrgljubljana/solver_output.hpp>

#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

namespace fs = std::filesystem;
using namespace NRG::Workflow;
unsigned assertions = 0;

void check(bool condition, const std::string &message) {
  ++assertions;
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &file) {
  std::ifstream in(file);
  check(bool(in), "Cannot read " + file.string());
  std::ostringstream text;
  text << in.rdbuf();
  return text.str();
}

template <typename Function> void rejects(Function function, const std::string &detail) {
  try {
    function();
  } catch (const std::runtime_error &error) {
    check(std::string(error.what()).find(detail) != std::string::npos,
          "Expected an error containing '" + detail + "', got '" + error.what() + "'");
    return;
  }
  throw std::runtime_error("Expected an error containing: " + detail);
}

std::size_t empty_lines(const std::string &text) {
  std::size_t count = 0;
  std::istringstream lines(text);
  for (std::string line; std::getline(lines, line);)
    if (line.find_first_not_of(" \t\r") == std::string::npos) count++;
  return count;
}

void stored_outputs(const fs::path &tests) {
  const std::map<std::string, std::string> lists{{"corr", "specs"},   {"spin", "spect"}, {"orbspin", "specot"},
                                                 {"spec", "specd"},   {"specq", "specq"}, {"gt", "specgt"},
                                                 {"i1t", "speci1t"}, {"i2t", "speci2t"}, {"chit", "specchit"}};
  std::size_t tables = 0, flows = 0, spectra = 0;
  for (const auto &test : fs::directory_iterator(tests)) {
    const auto ref = test.path() / "ref";
    if (!fs::is_directory(ref)) continue;
    for (const auto &entry : fs::directory_iterator(ref)) {
      const auto name = entry.path().filename().string();
      const auto where = entry.path().string();
      try {
        if (name == "td" || name == "custom" || name == "tdfdm" || name == "customfdm") {
          const auto table = read_step_table(read(entry.path()));
          check(!table.columns.empty() && !table.rows.empty() && table.columns.front() == "T", where + ": columns");
          if (name.ends_with("fdm")) check(read_final_values(read(entry.path())).contains("T"), where + ": final values");
          tables++;
        } else if (name == "annotated.dat") {
          const auto text = read(entry.path());
          const auto flow = read_energy_flow(text);
          check(flow.size() == empty_lines(text), where + ": one step per empty line");
          flows++;
        } else if (const auto dens = name.find("_dens_"); dens != std::string::npos && (name.ends_with(".dat") || name.ends_with(".bin"))) {
          // <prefix>_<algorithm>_dens_<correlator>[-u|-d].<extension>
          const auto head      = name.substr(0, dens);
          const auto underline = head.find('_');
          if (underline == std::string::npos || !lists.contains(head.substr(0, underline))) continue;
          const auto list      = lists.at(head.substr(0, underline));
          auto correlator      = name.substr(dens + 6, name.size() - dens - 6 - 4);
          int spin             = 0;
          const bool doublet   = list == "specd" || list == "specgt" || list == "speci1t" || list == "speci2t";
          if (doublet && (correlator.ends_with("-u") || correlator.ends_with("-d"))) {
            spin       = correlator.ends_with("-u") ? 1 : -1;
            correlator = correlator.substr(0, correlator.size() - 2);
          }
          check(spectrum_file_name(list, correlator, head.substr(underline + 1), spin) + name.substr(name.size() - 4) == name,
                where + ": spectrum file name");
          spectra++;
        }
      } catch (const std::runtime_error &error) { throw std::runtime_error(where + ": " + error.what()); }
    }
  }
  check(tables > 100 && flows > 50 && spectra > 100, "too few stored outputs found");
  std::cout << "  stored outputs: " << tables << " tables, " << flows << " energy flows, " << spectra << " spectrum names\n";
}

void mpitest1_qs(const fs::path &tests) {
  const auto ref = tests / "mpitest1_qs" / "ref";
  const auto td  = read_step_table(read(ref / "td"));
  check(td.columns.size() == 9 && td.columns[1] == "<Sz^2>" && td.rows.size() == 5, "td shape");
  check(td.rows[0][0] == 1.020139447 && td.rows[0][2] == -0.03678911966, "td values");
  const auto custom = read_step_table(read(ref / "custom"));
  check((custom.columns == std::vector<std::string>{"T", "n_d"}) && custom.rows.size() == 5, "custom: the header after the column numbers");
  const auto tdfdm = read_final_values(read(ref / "tdfdm"));
  check(tdfdm.at("T") == 0.1 && tdfdm.at("E_fdm") == -2.041644463, "tdfdm values");
  check(read_final_values(read(ref / "customfdm")).at("n_d") == 0.9624330014, "customfdm values");
  const auto flow = read_energy_flow(read(ref / "annotated.dat"));
  check(flow.size() == 5 && flow[0].size() == 10, "annotated.dat: steps and levels");
  check(flow[0][0] == EnergyLevel{0.0, {{0, 1}}, 1} && flow[0][1] == EnergyLevel{0.19078253, {{-1, 2}}, 2}, "annotated.dat: levels");
  check(spectrum_file_name("specd", "A_d-A_d", "FDM") == "spec_FDM_dens_A_d-A_d", "spectrum name");
  check(spectrum_file_name("specd", "A_d-A_d", "CFS", 1) == "spec_CFS_dens_A_d-A_d-u", "spectrum name with a field");
  std::cout << "  mpitest1_qs: ok\n";
}

void malformed() {
  rejects([] { read_step_table("1 2\n"); }, "line 1: data before a header line");
  rejects([] { read_step_table("# T x\n1 2\n3\n"); }, "line 3: expected 2 values, found 1");
  rejects([] { read_step_table("# T x\n1 two\n"); }, "line 2: not a number: \"two\"");
  rejects([] { read_step_table("# T x\n1 2\n# T y\n"); }, "line 3: a header line after the data");
  rejects([] { read_step_table(""); }, "no header line");
  rejects([] { read_final_values("# T x\n1 2\n3 4\n"); }, "exactly one row, found 2");
  rejects([] { read_energy_flow("0.5 (0 1)\n"); }, "line 1: expected [degeneracy]");
  rejects([] { read_energy_flow("0.5 0 1\n"); }, "dumpgroups=false");
  rejects([] { read_energy_flow("x (0 1) [1]\n"); }, "expected an energy");
  rejects([] { read_energy_flow("0.5 (0 1 [1]\n"); }, "unclosed");
  rejects([] { read_energy_flow("0.5 (0 a) [1]\n"); }, "not an integer quantum number");
  rejects([] { read_energy_flow("0.5 (0 1) [0]\n"); }, "expected [degeneracy]");
  rejects([] { spectrum_file_name("specx", "A-A", "FDM"); }, "unknown spectrum list");
  rejects([] { spectrum_file_name("specs", "n-n", "FDM", 1); }, "only the doublet lists");
  rejects([] { spectrum_file_name("specd", "A-A", "FDM", 2); }, "only the doublet lists");
  check(read_energy_flow("\n0 (0 1) [1]\n\n").size() == 2, "a step without levels");
  std::cout << "  malformed text: ok\n";
}

void floating_point_fields() {
  const auto table = read_step_table("# decimal exponent zero subnormal infinity nan\n.5 -1.25e+2 -0 5e-324 inf nan\n");
  const auto &row  = table.rows.at(0);
  check(row.at(0) == 0.5 && row.at(1) == -125.0, "decimal and scientific notation");
  check(row.at(2) == 0.0 && std::signbit(row.at(2)), "signed zero");
  check(row.at(3) == std::numeric_limits<double>::denorm_min(), "representable subnormal");
  check(std::isinf(row.at(4)) && std::isnan(row.at(5)), "infinity and NaN");
  for (const std::string token : {"+1", "-0X1p0", "1,25", "1.25x", "1e400", "1e-400"}) {
    rejects([&] { read_step_table("# T\n" + token + "\n"); }, "line 2: not a number: \"" + token + "\"");
  }
}

} // namespace

int main(int argc, char **argv) {
  try {
    check(argc == 2, "Usage: workflow_output <test/c++ directory>");
    stored_outputs(argv[1]);
    mpitest1_qs(argv[1]);
    malformed();
    floating_point_fields();
    std::cout << "Output readers: tests passed (" << assertions << " checks).\n";
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "Output readers: tests failed: " << error.what() << '\n';
    return 1;
  }
}
