#include <nrgljubljana/solver_workflow.hpp>

#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace fs = std::filesystem;
using namespace NRG::Workflow;
static_assert(workflow_api_version >= 2);

namespace {

void require(bool condition, const std::string &message) {
  if (!condition) throw std::runtime_error(message);
}

std::string read(const fs::path &path) {
  std::ifstream file(path);
  require(bool(file), "Cannot read " + path.string());
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void flush() {
  std::cout.flush();
  std::cerr.flush();
  require(std::fflush(nullptr) == 0, "Cannot flush test streams");
}

void emit(const std::string &label) {
  std::fprintf(stdout, "%s: C stdout\n", label.c_str());
  std::fprintf(stderr, "%s: C stderr\n", label.c_str());
  const auto out = label + ": C++ stdout\n", err = label + ": C++ stderr\n";
  std::cout.write(out.data(), out.size());
  std::cerr.write(err.data(), err.size());
  // Leave a buffered tail: capture must flush it before restoring the host fds.
  std::fputs("C tail\n", stdout);
  std::cout.write("C++ tail\n", 9);
}

void contains_all(const std::string &text, const std::string &label) {
  for (const auto *stream : {"C stdout", "C stderr", "C++ stdout", "C++ stderr"})
    require(text.find(label + ": " + stream + '\n') != std::string::npos, "Missing captured " + std::string(stream));
  require(text.find("C tail\n") != std::string::npos && text.find("C++ tail\n") != std::string::npos, "Unflushed diagnostic tail");
}

template <typename Function> void fails(Function function, const std::string &message) {
  try { function(); }
  catch (const std::exception &error) {
    require(std::string(error.what()).find(message) != std::string::npos, "Wrong failure: " + std::string(error.what()));
    return;
  }
  throw std::runtime_error("Expected failure: " + message);
}

struct TemporaryRoot {
  fs::path path;
  explicit TemporaryRoot(const fs::path &parent) {
    auto pattern = (fs::canonical(parent) / "workflow-diagnostics-XXXXXX").string();
    require(::mkdtemp(pattern.data()) != nullptr, "Cannot create test directory");
    path = pattern;
  }
  ~TemporaryRoot() { std::error_code error; fs::remove_all(path, error); }
};

struct HostStreams {
  const int out = ::dup(STDOUT_FILENO), err = ::dup(STDERR_FILENO);
  const int out_flags = ::fcntl(STDOUT_FILENO, F_GETFD), err_flags = ::fcntl(STDERR_FILENO, F_GETFD);
  std::unique_ptr<std::FILE, decltype(&std::fclose)> output{std::tmpfile(), &std::fclose}, error{std::tmpfile(), &std::fclose};
  HostStreams() {
    require(out != -1 && err != -1 && output && error, "Cannot save host streams");
    flush();
    require(::dup2(::fileno(output.get()), STDOUT_FILENO) != -1 && ::dup2(::fileno(error.get()), STDERR_FILENO) != -1,
            "Cannot redirect test streams");
    require(::fcntl(STDOUT_FILENO, F_SETFD, FD_CLOEXEC) == 0 && ::fcntl(STDERR_FILENO, F_SETFD, FD_CLOEXEC) == 0,
            "Cannot set host descriptor flags");
  }
  ~HostStreams() {
    try { flush(); } catch (...) {}
    ::dup2(out, STDOUT_FILENO);
    ::dup2(err, STDERR_FILENO);
    ::fcntl(STDOUT_FILENO, F_SETFD, out_flags);
    ::fcntl(STDERR_FILENO, F_SETFD, err_flags);
    ::close(out);
    ::close(err);
  }
};

struct StreamState {
  std::ostream &stream;
  const int fd;
  const std::ios_base::fmtflags flags = stream.flags();
  const std::streamsize precision = stream.precision(), width = stream.width();
  const char fill = stream.fill();
  std::streambuf *const buffer = stream.rdbuf();
  const std::ios_base::iostate state = stream.rdstate(), exceptions = stream.exceptions();
  const int descriptor_flags = ::fcntl(fd, F_GETFD), status_flags = ::fcntl(fd, F_GETFL);
  const off_t offset = ::lseek(fd, 0, SEEK_CUR);
  struct stat status {};
  StreamState(std::ostream &stream_, int fd_) : stream(stream_), fd(fd_) {
    require(::fstat(fd, &status) == 0, "Cannot inspect stream descriptor");
  }
  void check(bool silent) const {
    require(stream.flags() == flags && stream.precision() == precision && stream.width() == width && stream.fill() == fill &&
            stream.rdbuf() == buffer && stream.rdstate() == state && stream.exceptions() == exceptions, "C++ stream state leaked");
    struct stat after {};
    require(::fstat(fd, &after) == 0 && after.st_dev == status.st_dev && after.st_ino == status.st_ino &&
            ::fcntl(fd, F_GETFD) == descriptor_flags && ::fcntl(fd, F_GETFL) == status_flags, "Stream descriptor changed");
    if (silent) require(::lseek(fd, 0, SEEK_CUR) == offset, "Diagnostics leaked to host stream");
  }
};

} // namespace

// These two backends are deliberately test-only. Routing and discretization are
// production sources; no initializer, solver, MPI initialization or diagonalization runs.
namespace NRG::Workflow::detail {
double instantiate(const std::string &parameter, const std::string &, unsigned) {
  std::cout.precision(2);
  std::cerr.fill('?');
  emit("instantiate");
  if (parameter == "fail") throw std::runtime_error("initializer fixture failure");
  return -0.25;
}
} // namespace NRG::Workflow::detail

namespace NRG {
void run_nrg(const std::string &directory, const std::string &scratch) {
  const auto original = fs::current_path();
  try {
    fs::current_path(directory);
    emit("solve");
    if (scratch == "fail") throw std::runtime_error("solver fixture failure");
    fs::current_path(original);
  } catch (...) {
    fs::current_path(original);
    throw;
  }
}
} // namespace NRG

int main(int argc, char **argv) {
  try {
    require(argc == 2 || (argc == 3 && std::string(argv[2]) == "unsynced"), "Expected temporary parent [unsynced]");
    if (argc == 3) std::ios::sync_with_stdio(false);
    const TemporaryRoot root(argv[1]);
    const auto original = fs::current_path();
    const auto calculation = root.path / "calculation";
    fs::create_directory(calculation);
    fs::create_directory(calculation / "logs");
    const auto relative = fs::relative(calculation, original).string();
    std::ofstream(calculation / "param") << "[param]\nLambda=2\ndos=missing-density.dat\n";
    const HostStreams host;
    std::cout.precision(5);
    std::cerr.precision(3);
    std::cout.width(7);
    std::cerr.width(9);
    std::cout.fill('@');
    std::cerr.fill('!');
    auto checked = [&](auto function, bool silent = true) {
      flush();
      const StreamState out(std::cout, STDOUT_FILENO), err(std::cerr, STDERR_FILENO);
      function();
      flush();
      require(fs::current_path() == original, "CWD leaked");
      out.check(silent);
      err.check(silent);
      require(!std::ferror(stdout) && !std::ferror(stderr), "C stream error leaked");
    };

    // Real discretization reports its parameters before failing to open the DOS.
    // A relative calculation directory catches accidental resolution after chdir.
    checked([&] {
      emit("before");
      fails([&] { discretize({relative}, "logs/discretize.log"); }, "missing-density.dat");
    }, false);
    const auto diagnostic = read(calculation / "logs/discretize.log");
    require(diagnostic.find("# Lambda=") != std::string::npos, "Discretization produced no diagnostics");
    require(diagnostic.find("before") == std::string::npos, "Pre-call output entered diagnostic log");
    require(diagnostic.find("tail\n") == std::string::npos, "Pre-call buffered output entered diagnostic log");
    require(!fs::exists(calculation / "FSOL.dat"), "Failed discretization published output");

    checked([&] { fails([&] { discretize({relative}, "logs/discretize.log"); }, "Opening workflow diagnostic log"); });
    require(read(calculation / "logs/discretize.log") == diagnostic, "Existing log was truncated");
    for (const auto *target : {"discretize.log", "absent.log"}) {
      const auto link = calculation / "logs/link.log";
      fs::create_symlink(target, link);
      checked([&] { fails([&] { discretize({relative}, "logs/link.log"); }, "Opening workflow diagnostic log"); });
      require(fs::is_symlink(link), "Log symlink was replaced");
      fs::remove(link);
    }
    require(!fs::exists(calculation / "logs/absent.log"), "Dangling log symlink was followed");
    checked([&] { fails([&] { discretize({relative}, "absent/log"); }, "Opening workflow diagnostic log"); });
    require(!fs::exists(calculation / "absent"), "Missing log parent was created");
    checked([&] {
      fails([&] { discretize({(root.path / "absent").string()}, (root.path / "invalid-directory.log").string()); },
            "existing directory");
    });
    require(!fs::exists(root.path / "invalid-directory.log"), "Invalid preparation created a diagnostic log");

    checked([&] {
      const auto result = instantiate({relative}, "logs/instantiate.log");
      require(result.ground_energy == -0.25 && result.data_file == (calculation / "data").string(), "Initializer result changed");
    });
    contains_all(read(calculation / "logs/instantiate.log"), "instantiate");
    checked([&] { fails([&] { instantiate({relative, "fail"}, "logs/instantiate-fail.log"); }, "initializer fixture failure"); });
    contains_all(read(calculation / "logs/instantiate-fail.log"), "instantiate");
    checked([&] { solve(relative, "unused", "logs/solve.log"); });
    contains_all(read(calculation / "logs/solve.log"), "solve");
    checked([&] { fails([&] { solve(relative, "fail", "logs/solve-fail.log"); }, "solver fixture failure"); });
    contains_all(read(calculation / "logs/solve-fail.log"), "solve");
    const auto absolute_log = root.path / "absolute.log";
    checked([&] { solve(relative, "unused", absolute_log.string()); });
    contains_all(read(absolute_log), "solve");
    require((fs::status(absolute_log).permissions() & (fs::perms::group_all | fs::perms::others_all)) == fs::perms::none,
            "Diagnostic log is not private");

    // Empty paths and original ABI entry points must retain separate host streams.
    checked([&] {
      fails([&] { discretize({relative}); }, "missing-density.dat");
      fails([&] { discretize({relative}, ""); }, "missing-density.dat");
      instantiate({relative});
      instantiate({relative}, "");
      solve(relative, "unused");
      solve(relative, "unused", "");
    }, false);
    emit("after");
    flush();
    require(read(calculation / "logs/discretize.log") == diagnostic, "Post-call output entered diagnostic log");
    auto host_text = [](std::FILE *file) {
      require(std::fseek(file, 0, SEEK_SET) == 0, "Cannot rewind host stream");
      std::string text;
      for (int ch; (ch = std::fgetc(file)) != EOF;) text += static_cast<char>(ch);
      require(!std::ferror(file), "Cannot read host stream");
      return text;
    };
    const auto out = host_text(host.output.get()), err = host_text(host.error.get());
    for (const auto *label : {"before", "after", "instantiate", "solve"}) {
      for (const auto *language : {"C", "C++"}) {
        require(out.find(std::string(label) + ": " + language + " stdout\n") != std::string::npos, "Host stdout not restored");
        require(err.find(std::string(label) + ": " + language + " stderr\n") != std::string::npos, "Host stderr not restored");
      }
    }
    require(out.find("stderr") == std::string::npos && err.find("stdout") == std::string::npos, "Host streams were combined");
    require(out.find("# Lambda=") != std::string::npos, "Uncaptured discretization was silenced");
    for (const auto &entry : fs::directory_iterator(calculation))
      require(!entry.is_directory() || entry.path().filename() == "logs", "Discretization staging leaked");
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "workflow diagnostics: " << error.what() << '\n';
    return 1;
  }
}
