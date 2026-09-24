#include <gtest/gtest.h>

#include <cstddef>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <cerrno>
#include <csignal>
#include <filesystem>
#include <system_error>
#include <sys/resource.h>
#include <sys/wait.h>
#include <unistd.h>

#include <adapt/adapt.hpp>

using namespace NRG::Adapt;

namespace {

void write_file(const std::string &filename, const std::string &contents) {
  std::ofstream file(filename);
  file << contents;
}

std::string read_file(const std::string &filename) {
  std::ifstream file(filename);
  return {std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>()};
}

void expect_no_temporary_output(const std::string &filename) {
  for (const auto &entry : std::filesystem::directory_iterator("."))
    EXPECT_FALSE(entry.path().filename().string().starts_with(filename + ".tmp.")) << entry.path();
}

} // namespace

TEST(Adapt, parser_skips_blank_lines) { // NOLINT
  const auto filename = "adapt_blank_lines.param";
  write_file(filename, "[param]\n\nLambda=3\n\n# comment\n\nxmax=7\n");

  Params P(filename);
  EXPECT_EQ(P.P("Lambda", 0.0), 3.0);
  EXPECT_EQ(P.P("xmax", 0.0), 7.0);

  std::remove(filename);
}

TEST(Adapt, parser_accepts_extended_bool_values_and_rejects_invalid_values) { // NOLINT
  const auto filename = "adapt_bool.param";
  write_file(filename,
             "[param]\n"
             "yes = YeS\n"
             "true = tRuE\n"
             "one = 1\n"
             "no = nO\n"
             "false = FaLsE\n"
             "zero = 0\n"
             "invalid = maybe\n");

  Params P(filename);
  EXPECT_TRUE(P.Pbool("yes", false));
  EXPECT_TRUE(P.Pbool("true", false));
  EXPECT_TRUE(P.Pbool("one", false));
  EXPECT_FALSE(P.Pbool("no", true));
  EXPECT_FALSE(P.Pbool("false", true));
  EXPECT_FALSE(P.Pbool("zero", true));
  EXPECT_TRUE(P.Pbool("missing", true));
  EXPECT_THROW(P.Pbool("invalid", false), std::runtime_error);

  std::remove(filename);
}

TEST(Adapt, parser_rejects_partial_and_nonfinite_numbers) { // NOLINT
  const auto filename = "adapt_invalid_numbers.param";
  write_file(filename,
             "[param]\n"
             "partial=2junk\n"
             "infinite=inf\n"
             "integer=12junk\n");

  Params params(filename);
  EXPECT_THROW(params.P("partial", 0.0), std::invalid_argument);
  EXPECT_THROW(params.P("infinite", 0.0), std::invalid_argument);
  EXPECT_THROW(params.Pint("integer", 0), std::invalid_argument);

  std::remove(filename);
}

TEST(Adapt, linint_requires_two_points) { // NOLINT
  EXPECT_THROW(LinInt(Vec{}), std::runtime_error);
  EXPECT_THROW(LinInt(Vec{{1.0, 2.0}}), std::runtime_error);
}

TEST(Adapt, int_with_to_throws_when_step_limit_is_exceeded) { // NOLINT
  const auto filename = "adapt_int_with_to.param";
  write_file(filename, "[param]\nLambda=2\nadapt=false\nxmax=2\nxfine=2\n");

  Params P(filename);
  Adapt calc(P, Sign::POS);
  calc.rho = NRG::Tools::TabulatedDensity({{0.0, 1.0}, {1.0, 1.0}},
                                          NRG::Tools::InterpolationMethod::linear);
  calc.x = 0.0;
  calc.y = 0.0;

  EXPECT_THROW(calc.int_with_to(0.0, 1.0, []([[maybe_unused]] const auto x, [[maybe_unused]] const auto y) { return 0.0; }, false, 1e-10, 0), std::runtime_error);

  std::remove(filename);
}

TEST(Adapt, split_string_keeps_high_bit_bytes) { // NOLINT
  const std::string input = std::string(1, static_cast<char>(0x80)) + "1 2";
  const auto columns = split_string(input, 2);

  ASSERT_EQ(columns.size(), 2U);
  EXPECT_EQ(columns[0], std::string(1, static_cast<char>(0x80)) + "1");
  EXPECT_EQ(columns[1], "2");
}

TEST(Adapt, integral_method_matches_flat_band) { // NOLINT
  const auto filename = "adapt_integral_flat.param";
  write_file(filename,
             "[param]\n"
             "Lambda=10\n"
             "adapt=false\n"
             "f_method=integral\n"
             "xmax=30\n"
             "outputstep=0.25\n"
             "allowed_error=1e-10\n");

  Params P(filename);
  Adapt calc(P, Sign::POS, 0.01);
  calc.run();

  const auto values = load_g("FSOL.dat");
  const double lambda = 10.0;
  const double log_lambda = std::log(lambda);
  for (const auto &[x, f] : values) {
    const double expected = x <= 2.0
                              ? 2.0 - x + (1.0 - std::pow(lambda, 1.0 - x)) / log_lambda
                              : std::pow(lambda, 2.0 - x) * (1.0 - 1.0 / lambda) / log_lambda;
    const double expected_f = expected / std::pow(lambda, 2.0 - x);
    EXPECT_NEAR(f, expected_f, 1e-10) << "x=" << x;
  }

  std::remove("FSOL.dat");
  std::remove(filename);
}

TEST(Adapt, cumulative_inverse_uses_upper_plateau_edge) { // NOLINT
  const auto param_filename = "adapt_integral_plateau.param";
  const auto dos_filename = "adapt_integral_plateau.dat";
  write_file(param_filename,
             "[param]\n"
             "dos=adapt_integral_plateau.dat\n"
             "f_method=integral\n"
             "xmax=2\n");
  write_file(dos_filename,
             "0.1 0\n"
             "0.25 0\n"
             "0.4 1\n"
             "0.5 0\n"
             "0.7 0\n"
             "0.8 1\n"
             "1.0 1\n");

  Params P(param_filename);
  Adapt calc(P, Sign::POS);
  calc.load_init_rho();
  calc.init_cumulative();

  EXPECT_NEAR(calc.inverse_normalized_cumulative(0.0), 0.25, 1e-14);
  const auto plateau_weight = calc.normalized_cumulative(0.6);
  EXPECT_NEAR(calc.inverse_normalized_cumulative(plateau_weight), 0.7, 1e-14);
  EXPECT_LT(calc.inverse_normalized_cumulative(std::nextafter(plateau_weight, 0.0)), 0.5);
  EXPECT_GT(calc.inverse_normalized_cumulative(std::nextafter(plateau_weight, 1.0)), 0.7);

  std::remove(dos_filename);
  std::remove(param_filename);
}

TEST(Adapt, hardgap_integral_matches_flat_band_on_both_rescaled_branches) { // NOLINT
  const auto filename = "adapt_hardgap_flat.param";
  write_file(filename, "[param]\nLambda=2\nadapt=false\nhardgap=true\nboundary=0.25\n"
                       "f_method=integral\nxmax=6\noutputstep=0.25\n");
  Params params(filename);
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    for (const double bandwidth : {0.5, 2.0}) {
      SCOPED_TRACE(sign == Sign::POS ? "positive" : "negative");
      SCOPED_TRACE(bandwidth);
      params["bandrescale"] = std::to_string(bandwidth);
      Adapt calc(params, sign, 0.01);
      calc.run();
      const auto values = load_g(calc.f_fn(sign));
      ASSERT_EQ(values.size(), 21U);
      EXPECT_DOUBLE_EQ(values.front().first, 1.0);
      EXPECT_DOUBLE_EQ(values.front().second, 0.5);
      EXPECT_DOUBLE_EQ(values.back().first, 6.0);
      for (const auto &[x, f] : values) {
        SCOPED_TRACE(x);
        const double scale = std::pow(2.0, 2.0 - x);
        const double ungapped = x < 2.0 ? 2.0 - x + (1.0 - std::pow(2.0, 1.0 - x)) / std::log(2.0)
                                       : scale * 0.5 / std::log(2.0);
        const double expected = 0.25 + 0.75 * ungapped;
        EXPECT_NEAR(f, expected / scale, 1e-10);
        const double energy = calc.mesh.Eps(x, f);
        EXPECT_NEAR(energy * bandwidth, expected * bandwidth, 1e-10);
        if (x > 1.0) {
          const double lower = calc.eps(x + 1.0);
          const double upper = calc.eps(x);
          EXPECT_NEAR(lower, 0.25 + 0.75 * std::pow(2.0, 1.0 - x), 1e-14);
          EXPECT_NEAR(upper, x <= 2.0 ? 1.0 : 0.25 + 0.75 * scale, 1e-14);
          EXPECT_GT(lower, 0.25);
          EXPECT_LT(lower, upper);
          EXPECT_LE(upper, 1.0);
          EXPECT_GT(energy, lower);
          EXPECT_LT(energy, upper);
        }
      }
      std::remove(calc.f_fn(sign).c_str());
    }
  }
  std::remove(filename);
}

TEST(Adapt, hardgap_rejects_unsupported_modes_after_cli_override_without_touching_outputs) { // NOLINT
  const auto filename = "adapt_hardgap_modes.param";
  write_file(filename, "[param]\nhardgap=true\nboundary=0.25\n");
  Params params(filename);
  const std::string original = "existing table\n";
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    const auto f_filename = sign == Sign::POS ? "FSOL.dat" : "FSOLNEG.dat";
    const auto g_filename = sign == Sign::POS ? "GSOL.dat" : "GSOLNEG.dat";
    write_file(f_filename, original);
    write_file(g_filename, original);
    for (const bool adaptive : {false, true}) {
      params["adapt"] = adaptive ? "true" : "false";
      for (const auto method : {"ode", "integral"}) {
        params["f_method"] = method;
        for (const bool force_integral : {false, true}) {
          const bool supported = !adaptive && (force_integral || params["f_method"] == "integral");
          try {
            Adapt calc(params, sign, 0.01, force_integral);
            EXPECT_TRUE(supported);
            EXPECT_EQ(calc.f_method, FMethod::INTEGRAL);
          } catch (const std::invalid_argument &error) {
            EXPECT_FALSE(supported);
            const std::string message = error.what();
            EXPECT_NE(message.find("adapt=false"), std::string::npos);
            EXPECT_NE(message.find("--integral"), std::string::npos);
          }
          EXPECT_EQ(read_file(f_filename), original);
          EXPECT_EQ(read_file(g_filename), original);
        }
      }
    }
    std::remove(f_filename);
    std::remove(g_filename);
  }
  std::remove(filename);
}

TEST(Adapt, hardgap_max_abs_truncation_fails_and_preserves_previous_table) { // NOLINT
  const auto filename = "adapt_hardgap_max_abs.param";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\n"
                       "xmax=6\noutputstep=1\nmax_abs=1\n");
  Params params(filename);
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    Adapt calc(params, sign, 0.01);
    const auto output = calc.f_fn(sign);
    const std::string original = "previous valid table\n";
    write_file(output, original);
    try {
      calc.run();
      FAIL() << "Expected a hardgap truncation error";
    } catch (const std::runtime_error &error) {
      const std::string message = error.what();
      for (const auto text : {"x_last=3", "xmax=6", "max_abs=1", "increase max_abs", "reduce xmax"})
        EXPECT_NE(message.find(text), std::string::npos) << message;
    }
    EXPECT_EQ(read_file(output), original);
    calc.max_abs = 100.0;
    EXPECT_NO_THROW(calc.run());
    EXPECT_DOUBLE_EQ(load_g(output).back().first, 6.0);
    // Exceeding the coefficient bound at the requested endpoint does not truncate the table.
    calc.max_abs = 1.0;
    calc.xmax = 3.0;
    EXPECT_NO_THROW(calc.run());
    EXPECT_DOUBLE_EQ(load_g(output).back().first, 3.0);
    std::remove(output.c_str());
  }
  std::remove(filename);
}

TEST(Adapt, hardgap_rejects_collapsed_intervals_and_unresolvable_positive_weight_energies) { // NOLINT
  const auto filename = "adapt_hardgap_collapsed.param";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\n"
                       "xmax=3\noutputstep=1\nmax_abs=1e100\n");
  Params params(filename);
  Adapt calc(params, Sign::POS, 1.0);
  calc.load_init_rho();
  calc.init_cumulative();
  calc.max_error = 0.0;
  std::unique_ptr<gsl_integration_cquad_workspace, GslWorkspaceDeleter> workspace(gsl_integration_cquad_workspace_alloc(1000));
  ASSERT_TRUE(workspace);
  EXPECT_DOUBLE_EQ(calc.Eps_integral(1.0, workspace.get()), 1.0);
  EXPECT_GT(calc.eps(55.0), calc.mesh.boundary);
  EXPECT_LT(calc.eps(55.0), calc.eps(54.0));
  EXPECT_GT(calc.rho.integral(calc.eps(55.0), calc.eps(54.0)), 0.0);
  EXPECT_THROW(calc.Eps_integral(54.0, workspace.get()), std::runtime_error);
  EXPECT_DOUBLE_EQ(calc.eps(60.0), calc.mesh.boundary);
  EXPECT_THROW(calc.Eps_integral(60.0, workspace.get()), std::runtime_error);

  // A very large Lambda makes the first exported interval collapse onto the gap.
  params["Lambda"] = "1e20";
  Adapt collapsed(params, Sign::POS, 1.0);
  const std::string original = "previous valid table\n";
  write_file("FSOL.dat", original);
  try {
    collapsed.run();
    FAIL() << "Expected a collapsed hardgap interval error";
  } catch (const std::runtime_error &error) {
    EXPECT_NE(std::string(error.what()).find("interval collapsed"), std::string::npos);
  }
  EXPECT_EQ(read_file("FSOL.dat"), original);
  std::remove("FSOL.dat");
  std::remove(filename);
}

TEST(Adapt, hardgap_publication_replaces_symlink_without_touching_target) { // NOLINT
  const auto filename = "adapt_hardgap_symlink.param";
  const auto target = "adapt_hardgap_symlink_target.dat";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\nxmax=3\noutputstep=1\n");
  Params params(filename);
  const std::string original = "existing target\n";
  write_file(target, original);
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    Adapt calc(params, sign, 0.01);
    const auto output = calc.f_fn(sign);
    std::filesystem::create_symlink(target, output);
    ASSERT_NO_THROW(calc.run());
    EXPECT_FALSE(std::filesystem::is_symlink(output));
    EXPECT_TRUE(std::filesystem::is_regular_file(output));
    EXPECT_DOUBLE_EQ(load_g(output).back().first, 3.0);
    EXPECT_EQ(read_file(target), original);
    expect_no_temporary_output(output);
    std::remove(output.c_str());
  }
  std::remove(target);
  std::remove(filename);
}

TEST(Adapt, hardgap_publication_rename_failure_preserves_destination_and_removes_temporary) { // NOLINT
  const auto filename = "adapt_hardgap_rename.param";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\nxmax=3\noutputstep=1\n");
  Params params(filename);
  const std::string original = "existing directory contents\n";
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    Adapt calc(params, sign, 0.01);
    const auto output = calc.f_fn(sign);
    ASSERT_TRUE(std::filesystem::create_directory(output));
    const auto marker = output + "/previous";
    write_file(marker, original);
    try {
      calc.run();
      FAIL() << "Expected a publication rename error";
    } catch (const std::system_error &error) {
      EXPECT_NE(std::string(error.what()).find("Failed to rename temporary output to " + output), std::string::npos);
    }
    EXPECT_EQ(read_file(marker), original);
    expect_no_temporary_output(output);
    std::remove(marker.c_str());
    std::filesystem::remove(output);
  }
  std::remove(filename);
}

TEST(Adapt, hardgap_publication_write_and_close_failures_preserve_existing_table) { // NOLINT
  const auto filename = "adapt_hardgap_io_failure.param";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\nxmax=3\n");
  Params params(filename);
  const std::string original = "previous valid table\n";
  for (const auto sign : {Sign::POS, Sign::NEG}) {
    // Small output fails when fclose flushes it; output larger than the stdio buffer fails in fwrite.
    for (const auto step : {"1", "0.001953125"}) {
      SCOPED_TRACE(step);
      params["outputstep"] = step;
      Adapt calc(params, sign, 0.01);
      const auto output = calc.f_fn(sign);
      write_file(output, original);
      const std::string operation = params["outputstep"] == "1" ? "close" : "write";
      const auto child = ::fork();
      ASSERT_GE(child, 0);
      if (child == 0) {
        // Never alter the test runner's resource limits or signal handling.
        struct rlimit limit;
        if (::getrlimit(RLIMIT_FSIZE, &limit) != 0) ::_exit(2);
        limit.rlim_cur = 0;
        if (::signal(SIGXFSZ, SIG_IGN) == SIG_ERR || ::setrlimit(RLIMIT_FSIZE, &limit) != 0) ::_exit(2);
        std::cout.setstate(std::ios_base::badbit);
        try {
          calc.run();
        } catch (const std::system_error &error) {
          const bool expected = error.code() == std::errc::file_too_large
                                && std::string(error.what()).find("Failed to " + operation + " temporary output for " + output)
                                     != std::string::npos;
          ::_exit(expected ? 0 : 3);
        } catch (...) {
          ::_exit(4);
        }
        ::_exit(1);
      }
      int status = 0;
      pid_t waited;
      do {
        waited = ::waitpid(child, &status, 0);
      } while (waited == -1 && errno == EINTR);
      ASSERT_EQ(waited, child);
      ASSERT_TRUE(WIFEXITED(status));
      EXPECT_EQ(WEXITSTATUS(status), 0) << "1=unexpected success, 2=setup failure, 3=wrong I/O error, 4=non-I/O exception";
      EXPECT_EQ(read_file(output), original);
      expect_no_temporary_output(output);
      std::remove(output.c_str());
    }
  }
  std::remove(filename);
}

TEST(Adapt, hardgap_zero_weight_allows_an_edge_but_not_an_outside_plateau_inverse) { // NOLINT
  const auto filename = "adapt_hardgap_plateau.param";
  write_file(filename, "[param]\nLambda=2\nhardgap=true\nboundary=0.25\nf_method=integral\n");
  Params params(filename);
  Adapt calc(params, Sign::POS);
  calc.vecrho = {{0.0, 0.0}, {0.625, 0.0}, {1.0, 1.0}};
  calc.rho = NRG::Tools::TabulatedDensity(calc.vecrho);
  calc.init_cumulative();
  calc.max_error = 0.0;
  std::unique_ptr<gsl_integration_cquad_workspace, GslWorkspaceDeleter> workspace(gsl_integration_cquad_workspace_alloc(1000));
  ASSERT_TRUE(workspace);
  EXPECT_DOUBLE_EQ(calc.rho.integral(calc.eps(4.0), calc.eps(3.0)), 0.0);
  EXPECT_DOUBLE_EQ(calc.Eps_integral(3.0, workspace.get()), 0.625);
  EXPECT_DOUBLE_EQ(calc.rho.integral(calc.eps(4.5), calc.eps(3.5)), 0.0);
  EXPECT_THROW(calc.Eps_integral(3.5, workspace.get()), std::runtime_error);
  std::remove(filename);
}

TEST(Adapt, zero_or_disabled_hardgap_keeps_ungapped_modes_and_truncation) { // NOLINT
  const auto filename = "adapt_hardgap_ungapped.param";
  write_file(filename, "[param]\nLambda=2\nxmax=6\noutputstep=0.25\nmax_abs=0.1\n");
  Params params(filename);
  for (const bool hardgap : {false, true}) {
    params["hardgap"] = hardgap ? "true" : "false";
    params["boundary"] = hardgap ? "0" : "0.25";
    for (const bool adaptive : {false, true}) {
      params["adapt"] = adaptive ? "true" : "false";
      for (const auto method : {"ode", "integral"}) {
        params["f_method"] = method;
        EXPECT_NO_THROW(Adapt(params, Sign::POS, 0.01));
      }
    }
    params["adapt"] = "false";
    Adapt calc(params, Sign::POS, 0.01);
    EXPECT_NO_THROW(calc.run());
    const auto values = load_g("FSOL.dat");
    ASSERT_EQ(values.size(), 2U);
    EXPECT_DOUBLE_EQ(values.back().first, 1.25);
    const double expected = 0.75 + (1.0 - std::pow(2.0, -0.25)) / std::log(2.0);
    EXPECT_NEAR(values.back().second * std::pow(2.0, 0.75), expected, 1e-10);
    std::remove("FSOL.dat");
  }
  std::remove(filename);
}

TEST(Adapt, cumulative_inverse_extends_terminal_plateau_to_band_edge) { // NOLINT
  const auto param_filename = "adapt_integral_terminal_plateau.param";
  const auto dos_filename = "adapt_integral_terminal_plateau.dat";
  write_file(param_filename,
             "[param]\n"
             "dos=adapt_integral_terminal_plateau.dat\n"
             "f_method=integral\n"
             "xmax=2\n");
  write_file(dos_filename,
             "0.1 1\n"
             "0.5 1\n"
             "0.8 0\n"
             "0.9 0\n");

  Params P(param_filename);
  Adapt calc(P, Sign::POS);
  calc.load_init_rho();
  calc.init_cumulative();

  EXPECT_DOUBLE_EQ(calc.inverse_normalized_cumulative(1.0), 1.0);

  std::remove(dos_filename);
  std::remove(param_filename);
}

TEST(Adapt, integral_method_matches_linear_density_oracle) { // NOLINT
  const auto filename = "adapt_integral_linear.param";
  write_file(filename,
             "[param]\n"
             "Lambda=2\n"
             "adapt=false\n"
             "f_method=integral\n"
             "xmax=40\n"
             "allowed_error=1e-10\n");

  Params P(filename);
  Adapt calc(P, Sign::POS);
  calc.vecrho = {{0.0, 0.0}, {1.0, 2.0}};
  calc.rho = NRG::Tools::TabulatedDensity(calc.vecrho, NRG::Tools::InterpolationMethod::linear);
  calc.init_cumulative();
  calc.max_error = 0.0;

  constexpr std::size_t limit = 1000;
  std::unique_ptr<gsl_integration_cquad_workspace, GslWorkspaceDeleter> workspace(
    gsl_integration_cquad_workspace_alloc(limit));
  ASSERT_TRUE(workspace);
  const double log_lambda = std::log(2.0);
  for (const double x : {1.25, 2.0, 10.0, 20.0, 30.0, 40.0}) {
    const double expected_weight = x < 2.0
                                     ? 2.0 - x + (1.0 - std::pow(2.0, 2.0 - 2.0 * x)) / (2.0 * log_lambda)
                                     : std::pow(2.0, 4.0 - 2.0 * x) * (1.0 - 0.25) / (2.0 * log_lambda);
    const double expected = std::sqrt(expected_weight);
    EXPECT_NEAR(calc.Eps_integral(x, workspace.get()) / expected, 1.0, 1e-9) << "x=" << x;
  }

  std::remove(filename);
}

TEST(Adapt, rejects_unknown_f_method) { // NOLINT
  const auto filename = "adapt_invalid_f_method.param";
  write_file(filename, "[param]\nf_method=unknown\n");

  Params P(filename);
  EXPECT_THROW(Adapt(P, Sign::POS), std::invalid_argument);

  std::remove(filename);
}

TEST(Adapt, selects_only_supported_density_interpolation) { // NOLINT
  const auto steffen_filename = "adapt_steffen.param";
  write_file(steffen_filename, "[param]\ndensity_interpolation=steffen\n");
  Params steffen_params(steffen_filename);
  Adapt steffen(steffen_params, Sign::POS);
  EXPECT_EQ(steffen.density_interpolation, NRG::Tools::InterpolationMethod::steffen);

  const auto invalid_filename = "adapt_invalid_density_interpolation.param";
  write_file(invalid_filename, "[param]\ndensity_interpolation=cspline\n");
  Params invalid_params(invalid_filename);
  EXPECT_THROW(Adapt(invalid_params, Sign::POS), std::invalid_argument);

  std::remove(invalid_filename);
  std::remove(steffen_filename);
}

TEST(Adapt, validates_bandrescale_and_uses_exact_density_weight) { // NOLINT
  const auto invalid_filename = "adapt_invalid_bandrescale.param";
  write_file(invalid_filename, "[param]\nbandrescale=0\n");
  Params invalid_params(invalid_filename);
  EXPECT_THROW(Adapt(invalid_params, Sign::POS), std::invalid_argument);
  std::remove(invalid_filename);

  const auto filename = "adapt_exact_density_weight.param";
  write_file(filename, "[param]\n");
  Params params(filename);
  Adapt calc(params, Sign::POS);
  calc.vecrho = {{0.25, 2.0}, {0.5, 2.0}};
  calc.rho = NRG::Tools::TabulatedDensity(calc.vecrho, NRG::Tools::InterpolationMethod::linear);
  calc.init_A();
  EXPECT_DOUBLE_EQ(calc.intA, 2.0);
  EXPECT_DOUBLE_EQ(calc.A, 2.0);
  std::remove(filename);
}

TEST(Adapt, cquad_error_policies_handle_failed_results) { // NOLINT
  const auto filename = "adapt_cquad_policy.param";
  write_file(filename, "[param]\nf_method=integral\n");
  Params params(filename);

  Adapt default_fail(params, Sign::POS);
  EXPECT_THROW(default_fail.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0), std::runtime_error);

  CquadOptions ignore_options;
  ignore_options.gsl_error_policy = NRG::Tools::GslErrorPolicy::ignore;
  Adapt ignored(params, Sign::POS, std::nullopt, false, ignore_options);
  EXPECT_NO_THROW(ignored.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0));

  CquadOptions warn_options;
  warn_options.gsl_error_policy = NRG::Tools::GslErrorPolicy::warn;
  Adapt warned(params, Sign::POS, std::nullopt, false, warn_options);
  testing::internal::CaptureStderr();
  EXPECT_NO_THROW(warned.handle_cquad_result(GSL_EMAXITER, 1.0, 0.1, 2.0));
  const auto warning = testing::internal::GetCapturedStderr();
  EXPECT_NE(warning.find("adapt: warning: Integral method failed at x=2.000000"), std::string::npos);

  std::remove(filename);
}

TEST(Adapt, cquad_override_does_not_bypass_allowed_error_validation) { // NOLINT
  const auto filename = "adapt_invalid_allowed_error.param";
  write_file(filename, "[param]\nadapt=true\nf_method=integral\nallowed_error=-1\n");
  Params params(filename);
  CquadOptions options;
  options.epsrel = 1e-8;

  EXPECT_THROW(Adapt(params, Sign::POS, std::nullopt, false, options), std::invalid_argument);
  std::remove(filename);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
