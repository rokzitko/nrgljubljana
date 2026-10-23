#include <nrgljubljana/solver_workflow.hpp>
#include "../adapt/adapt.hpp"
#include <workdir.hpp>

namespace NRG::Workflow::detail {

void discretize(const DiscretizationOptions &options) {
  NRG::Adapt::Sign sign;
  switch (options.branch) {
    case Branch::positive: sign = NRG::Adapt::Sign::POS; break;
    case Branch::negative: sign = NRG::Adapt::Sign::NEG; break;
    default: throw std::invalid_argument("Unknown discretization branch.");
  }
  const auto directory = std::filesystem::current_path();
  NRG::Workdir staging(directory.string(), true);
  NRG::Adapt::Params parameters(options.parameter_file);
  NRG::Adapt::Adapt calculation(parameters, sign, std::nullopt, options.integral, {}, staging.get());
  calculation.run();
  // Publish only completed files, replacing destination symlinks rather than following them.
  for (const auto &name : {calculation.g_fn(sign), calculation.f_fn(sign)}) {
    const auto staged = std::filesystem::path(staging.get()) / name;
    if (std::filesystem::exists(staged)) std::filesystem::rename(staged, directory / name);
  }
}

} // namespace NRG::Workflow::detail
