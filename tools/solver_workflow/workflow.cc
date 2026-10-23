#include <nrgljubljana/solver_workflow.hpp>

#include <exception>
#include <filesystem>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <type_traits>
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

template <typename Function>
auto prepare(const std::string &directory, Function function) {
  const std::lock_guard lock(invocation_mutex);
  if (directory.empty() || !std::filesystem::is_directory(directory))
    throw std::invalid_argument("Workflow calculation_directory must be an existing directory.");
  const SavedState saved_state;
  const auto original = std::filesystem::current_path();
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
  prepare(options.calculation_directory, [&] { detail::discretize(options); });
}

InitializerResult instantiate(const InitializerOptions &options) {
  return prepare(options.calculation_directory, [&] {
    const auto energy = detail::instantiate(options.parameter_file, options.template_directory, options.nmax);
    return InitializerResult{std::filesystem::absolute("data").string(), energy};
  });
}

void solve(const std::string &calculation_directory, const std::string &absolute_scratch_directory) {
  const std::lock_guard lock(invocation_mutex);
  const SavedState saved_state;
  NRG::run_nrg(calculation_directory, absolute_scratch_directory);
}

} // namespace NRG::Workflow
