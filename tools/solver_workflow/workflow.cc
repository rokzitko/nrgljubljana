#include <nrgljubljana/solver_workflow.hpp>

#include <cerrno>
#include <cstdio>
#include <exception>
#include <fcntl.h>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <system_error>
#include <type_traits>
#include <unistd.h>
#include <boost/io/ios_state.hpp>
#include <gmp.h>

#include <invar.hpp>
#include <nrg-lib.hpp>
#include "../instantiate/workflow.hpp"

namespace NRG::Workflow {
namespace detail {
void discretize(const DiscretizationOptions &);
}
namespace {

std::mutex invocation_mutex;

struct SavedState {
  const mp_bitcnt_t precision = mpf_get_default_prec();
  boost::io::ios_base_all_saver cout_format{std::cout}, cerr_format{std::cerr};
  boost::io::ios_fill_saver cout_fill{std::cout}, cerr_fill{std::cerr};
  const size_t invdim = Invar::invdim;
  decltype(Invar::qntype) qntype = Invar::qntype;
  decltype(Invar::names) names = Invar::names;
  ~SavedState() {
    mpf_set_default_prec(precision);
    Invar::invdim = invdim;
    Invar::qntype.swap(qntype);
    Invar::names.swap(names);
  }
};

class DiagnosticLog {
  int original_out = -1, original_error = -1;
  int out_flags = 0, error_flags = 0;

  static int duplicate_to(int source, int destination) noexcept {
    int result;
    do { result = ::dup2(source, destination); } while (result == -1 && errno == EINTR);
    return result;
  }

  void restore() noexcept {
    if (original_out != -1) {
      duplicate_to(original_out, STDOUT_FILENO);
      ::fcntl(STDOUT_FILENO, F_SETFD, out_flags);
      ::close(original_out);
      original_out = -1;
    }
    if (original_error != -1) {
      duplicate_to(original_error, STDERR_FILENO);
      ::fcntl(STDERR_FILENO, F_SETFD, error_flags);
      ::close(original_error);
      original_error = -1;
    }
  }

 public:
  DiagnosticLog(const std::string &directory, const std::string &file) {
    if (file.empty()) return;
    // Constructed under invocation_mutex, before preparation/run_nrg changes CWD.
    const auto path = std::filesystem::absolute(std::filesystem::path(directory) / file);
    std::cout.flush();
    std::cerr.flush();
    if (std::fflush(nullptr)) throw std::runtime_error("Cannot flush workflow diagnostic streams");
    if ((out_flags = ::fcntl(STDOUT_FILENO, F_GETFD)) == -1 || (error_flags = ::fcntl(STDERR_FILENO, F_GETFD)) == -1)
      throw std::system_error(errno, std::generic_category(), "Reading workflow diagnostic descriptors");
    const int destination = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
    if (destination == -1) throw std::system_error(errno, std::generic_category(), "Opening workflow diagnostic log " + path.string());
    if ((original_out = ::fcntl(STDOUT_FILENO, F_DUPFD_CLOEXEC, 3)) == -1 ||
        (original_error = ::fcntl(STDERR_FILENO, F_DUPFD_CLOEXEC, 3)) == -1 || duplicate_to(destination, STDOUT_FILENO) == -1 ||
        duplicate_to(destination, STDERR_FILENO) == -1) {
      const int error = errno;
      ::close(destination);
      restore();
      throw std::system_error(error, std::generic_category(), "Redirecting workflow diagnostics");
    }
    ::close(destination);
  }
  DiagnosticLog(const DiagnosticLog &) = delete;
  DiagnosticLog &operator=(const DiagnosticLog &) = delete;
  ~DiagnosticLog() {
    if (original_out == -1) return;
    try { std::cout.flush(); } catch (...) {}
    try { std::cerr.flush(); } catch (...) {}
    std::fflush(nullptr);
    restore();
  }
};

template <typename Function>
auto prepare(const std::string &directory, const std::string &diagnostic_file, Function function) {
  const std::lock_guard lock(invocation_mutex);
  if (directory.empty() || !std::filesystem::is_directory(directory))
    throw std::invalid_argument("Workflow calculation_directory must be an existing directory.");
  const SavedState saved_state;
  const auto original = std::filesystem::current_path();
  const DiagnosticLog diagnostics(directory, diagnostic_file);
  try {
    std::filesystem::current_path(directory);
    if constexpr (std::is_void_v<decltype(function())>) {
      function();
      std::filesystem::current_path(original);
    } else {
      auto result = function();
      std::filesystem::current_path(original);
      return result;
    }
  } catch (...) {
    std::error_code error;
    std::filesystem::current_path(original, error);
    if (error)
      std::throw_with_nested(std::filesystem::filesystem_error("Can't restore workflow caller's directory", original, error));
    throw;
  }
}

} // namespace

void discretize(const DiscretizationOptions &options) {
  discretize(options, {});
}

void discretize(const DiscretizationOptions &options, const std::string &diagnostic_file) {
  prepare(options.calculation_directory, diagnostic_file, [&] { detail::discretize(options); });
}

InitializerResult instantiate(const InitializerOptions &options) {
  return instantiate(options, {});
}

InitializerResult instantiate(const InitializerOptions &options, const std::string &diagnostic_file) {
  return prepare(options.calculation_directory, diagnostic_file, [&] {
    const auto energy = detail::instantiate(options.parameter_file, options.template_directory, options.nmax);
    return InitializerResult{std::filesystem::absolute("data").string(), energy};
  });
}

void solve(const std::string &calculation_directory, const std::string &absolute_scratch_directory) {
  solve(calculation_directory, absolute_scratch_directory, {});
}

void solve(const std::string &calculation_directory, const std::string &absolute_scratch_directory, const std::string &diagnostic_file) {
  const std::lock_guard lock(invocation_mutex);
  const SavedState saved_state;
  const DiagnosticLog diagnostics(calculation_directory, diagnostic_file);
  NRG::run_nrg(calculation_directory, absolute_scratch_directory);
}

} // namespace NRG::Workflow
