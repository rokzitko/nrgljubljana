#include <nrgljubljana/solver_workflow.hpp>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <gmp.h>
#include <gsl/gsl_errno.h>
#include <invar.hpp>
#include "../../common/parser.hpp"

namespace fs = std::filesystem;
using namespace NRG::Workflow;

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &path) {
  std::ifstream file(path);
  require(bool(file), "Cannot read " + path.string());
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void gsl_handler(const char *, const char *, int, int) {}

struct StreamState {
  std::ostream &stream;
  const std::ios_base::fmtflags flags = stream.flags();
  const std::streamsize precision = stream.precision(), width = stream.width();
  const char fill = stream.fill();
  void check() const {
    require(stream.flags() == flags && stream.precision() == precision && stream.width() == width && stream.fill() == fill,
            "Stream formatting leaked");
  }
};

class FailingParameterBuffer : public std::stringbuf {
  const bool fail_seek;
  int_type underflow() override { throw std::ios_base::failure("Injected parameter read failure"); }
  pos_type seekoff(off_type offset, std::ios_base::seekdir direction, std::ios_base::openmode mode) override {
    return fail_seek ? pos_type(off_type(-1)) : std::stringbuf::seekoff(offset, direction, mode);
  }
 public:
  FailingParameterBuffer(const std::string &prefix, bool fail_seek_) : std::stringbuf(prefix), fail_seek(fail_seek_) {}
};

int main(int argc, char **argv) {
  try {
    require(argc == 2, "Expected isolated test root");
    const auto root = fs::canonical(argv[1]);
    const auto cwd = fs::current_path();
    std::cout << std::scientific << std::setprecision(5) << std::setfill('@');
    std::cerr << std::fixed << std::setprecision(3) << std::setfill('!');
    std::cout.width(7);
    std::cerr.width(9);
    const StreamState out{std::cout}, err{std::cerr};
    mpf_set_default_prec(192);
    const auto precision = mpf_get_default_prec();
    auto *old_handler = gsl_set_error_handler(gsl_handler);
    // Deliberately unlike QS: prove that both preparation and solve restore the host's metadata.
    NRG::Invar::invdim = 3;
    NRG::Invar::qntype = {NRG::additive, NRG::multiplicative, NRG::mod3};
    NRG::Invar::names = {{"host_charge", 0}, {"host_parity", 1}, {"host_mod3", 2}};
    const auto invdim = NRG::Invar::invdim;
    const auto qntype = NRG::Invar::qntype;
    const auto names = NRG::Invar::names;
    auto check_state = [&] {
      require(fs::current_path() == cwd, "CWD leaked");
      require(mpf_get_default_prec() == precision, "GMP precision leaked");
      require(gsl_set_error_handler(gsl_handler) == gsl_handler, "GSL handler leaked");
      require(NRG::Invar::invdim == invdim && NRG::Invar::qntype == qntype && NRG::Invar::names == names,
              "Invariant metadata leaked");
      out.check();
      err.check();
    };
    auto fails = [&](auto function) {
      bool failed = false;
      std::string message;
      try { function(); } catch (const std::exception &error) { failed = true; message = error.what(); }
      require(failed, "Expected failure did not occur");
      check_state();
      return message;
    };
    // Inject failures after an actual seek, during the block search and after reading a key.
    for (const auto &[prefix, seek_failure] : {std::pair{"", true}, {"# comment\n", false}, {"[param]\nLambda=2\n", false}}) {
      FailingParameterBuffer buffer(prefix, seek_failure);
      std::ifstream input;
      static_cast<std::istream &>(input).rdbuf(&buffer);
      const auto message = fails([&] {
        std::map<std::string, std::string> parameters;
        if (NRG::Tools::find_block(input, "param")) NRG::Tools::parse_key_value_block(input, parameters);
      });
      require(message.find(seek_failure ? "Error seeking parameter input" : "Error reading parameter input") != std::string::npos,
              "Wrong parameter I/O failure");
    }
    // Missing blocks, empty files, EOF without a newline, and stopping at the next section retain their grammar.
    for (const auto &text : {"", "[other]\nLambda=4", " [param]\nLambda=4"}) {
      std::stringbuf buffer(text);
      std::ifstream input;
      static_cast<std::istream &>(input).rdbuf(&buffer);
      require(!NRG::Tools::find_block(input, "param"), "Missing-block grammar changed");
      require(!NRG::Tools::find_block(input, "param"), "EOF rewind failed");
    }
    for (const auto &text : {"[param]\nLambda=4", "[param]\nLambda=4\n[other]\nLambda=8\n"}) {
      std::stringbuf buffer(text);
      std::ifstream input;
      static_cast<std::istream &>(input).rdbuf(&buffer);
      std::map<std::string, std::string> parameters;
      require(NRG::Tools::find_block(input, "param"), "Parameter block not found");
      NRG::Tools::parse_key_value_block(input, parameters);
      require(parameters.size() == 1 && parameters.at("Lambda") == "4", "Parameter-block grammar changed");
    }
    auto run = [&](const std::string &name, unsigned nmax) {
      const auto directory = root / name;
      // Exercise calculation-relative custom paths, not paths relative to the parameter file.
      for (auto branch : {Branch::positive, Branch::negative}) {
        discretize({fs::relative(directory, cwd).string(), "inputs/solver.param", branch});
        check_state();
      }
      auto result = instantiate({directory.string(), "inputs/solver.param", "template", nmax});
      check_state();
      require(fs::path(result.data_file) == directory / "data", "Incorrect data_file result");
      require(std::isfinite(result.ground_energy), "Nonfinite ground energy");
      return result;
    };

    const auto parameter_a = read(root / "A/inputs/solver.param");
    const auto a = run("A", 6);
    const auto a_data = read(a.data_file);
    const auto a_fsol = read(root / "A/FSOL.dat");
    const auto a_negative = read(root / "A/FSOLNEG.dat");
    const auto b = run("B", 5);
    require(read(b.data_file) != a_data, "A/B inputs did not exercise distinct calculations");
    auto again = run("A", 6);
    require(read(again.data_file) == a_data && again.ground_energy == a.ground_energy, "A/B/A initializer contamination");
    require(read(root / "A/FSOL.dat") == a_fsol && read(root / "A/FSOLNEG.dat") == a_negative,
            "A/B/A discretization contamination");

    for (const unsigned nmax : {0U, 999U})
      fails([&] { instantiate({(root / "A").string(), "inputs/solver.param", "template", nmax}); });
    fails([&] { instantiate({(root / "A").string(), "inputs/solver.param", "broken-template", 6}); });
    require(read(a.data_file) == a_data, "Failure replaced valid data");
    fails([&] { discretize({(root / "A").string(), "bad-adapt", Branch::positive}); });
    require(fails([&] { discretize({(root / "A").string(), "."}); }).find("regular file") != std::string::npos,
            "Directory parameter input was not rejected");
#if defined(__linux__)
    // Linux reports this as a regular file, but reading unmapped address zero fails with EIO.
    require(fs::is_regular_file("/proc/self/mem"), "Read-error fixture is not a regular file");
    require(fails([&] { discretize({(root / "A").string(), "/proc/self/mem"}); }).find("Error reading parameter input") != std::string::npos,
            "File read failure was interpreted as defaults");
#endif
    require(read(root / "A/FSOL.dat") == a_fsol, "Failure replaced valid FSOL");
    require(read(root / "A/FSOLNEG.dat") == a_negative, "Failure replaced valid negative FSOL");
    fails([&] { discretize({(root / "missing").string()}); });
    fails([&] { discretize({(root / "A").string(), "inputs/solver.param", static_cast<Branch>(7)}); });
    fails([&] { instantiate({(root / "A").string(), "unsupported", "template", 6}); });
    fails([&] { instantiate({(root / "A").string(), "inputs/solver.param", "complex-template", 6}); });
    require(!fs::exists(root / "missing"), "Missing calculation directory was created");
    again = run("A", 6);
    require(read(again.data_file) == a_data && again.ground_energy == a.ground_energy, "Failure contaminated next call");
    require(read(root / "A/inputs/solver.param") == parameter_a, "Initializer edited solver parameters");

    // Explicit false retains f_method in the file, including integral selection.
    discretize({(root / "A").string(), "integral-param", Branch::positive, false});
    check_state();
    require(read(root / "A/FSOL.dat") == a_fsol, "integral=false ignored the parameter method");
    discretize({(root / "A").string(), "ode-param", Branch::positive, false});
    check_state();
    fs::rename(root / "A/FSOL.dat", root / "A/FSOL-ode.dat");
    discretize({(root / "A").string(), "inputs/solver.param"});
    check_state();

    fails([&] { solve((root / "A").string(), "relative-scratch"); });
    fails([&] { solve((root / "empty").string(), (root / "scratch").string()); });
    solve((root / "A").string(), (root / "scratch").string());
    check_state();
    require(fs::exists(root / "A/DONE"), "Solver did not produce DONE");
    require(fs::is_empty(root / "scratch"), "Solver scratch leaked");
    {
      std::ofstream data(a.data_file, std::ios::app);
      data << "\n? invalid-workflow-block\n";
    }
    const auto late_failure = fails([&] { solve((root / "A").string(), (root / "scratch").string()); });
    require(late_failure.find("Unknown block ? in data file") != std::string::npos,
            "Solver failure did not reach the block reader after symmetry initialization");
    require(!fs::exists(root / "A/DONE"), "Late failure left a stale DONE marker");
    require(fs::is_empty(root / "scratch"), "Late solver failure leaked scratch");
    again = run("A", 6);
    require(read(again.data_file) == a_data, "Solver contaminated next preparation");
    std::ofstream result(root / "A/result");
    result << std::setprecision(17) << a.ground_energy << '\n';
    gsl_set_error_handler(old_handler);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "workflow contract: " << error.what() << '\n';
    return 1;
  }
}
