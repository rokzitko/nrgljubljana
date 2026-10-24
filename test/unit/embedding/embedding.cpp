#include <gtest/gtest.h>
#include <nrg-lib.hpp>
#include <portabil.hpp>
#include <boost/io/ios_state.hpp>
#include <gmp.h>
#include <array>
#include <complex>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <sstream>
#include <thread>

using namespace NRG;
namespace fs = std::filesystem;

TEST(EmbeddingAssertions, Finiteness) {
  EXPECT_EQ(assert_isfinite(1.25), 1.25);
  EXPECT_EQ(assert_isfinite(std::complex<double>(1.0, 2.0)), std::complex<double>(1.0, 2.0));
  const auto infinity = std::numeric_limits<double>::infinity();
  const auto nan = std::numeric_limits<double>::quiet_NaN();
  for (const auto value : {infinity, -infinity, nan}) {
    try {
      finite_test_fnc(value, "embedding-input", 42);
      FAIL() << "Expected a finiteness exception";
    } catch (const std::runtime_error &error) {
      EXPECT_EQ(std::string(error.what()),
                "#### EXITING DUE TO FAILED ASSERTION.\nFile embedding-input, line 42.\nFiniteness assertion\n");
    }
    EXPECT_THROW(assert_isfinite(std::complex<double>(value, 0.0)), std::runtime_error);
    EXPECT_THROW(assert_isfinite(std::complex<double>(0.0, value)), std::runtime_error);
  }
}

TEST(EmbeddingAssertions, Equality) {
  EXPECT_NO_THROW(my_assert_equal(2, 2));
  const int left = 2, right = 3;
  try {
    my_assert_equal(left, right);
    FAIL() << "Expected an equality exception";
  } catch (const std::runtime_error &error) {
    const std::string message = error.what();
    EXPECT_NE(message.find("#### EXITING DUE TO FAILED ASSERTION.\nFile " + std::string(__FILE__) + ", line "), std::string::npos);
    EXPECT_NE(message.find("left=2\nright=3\nExiting.\n"), std::string::npos);
  }
}

TEST(EmbeddingAssertions, NotReached) {
  try {
    my_assert_not_reached();
  } catch (const std::runtime_error &error) {
    const std::string message = error.what();
    EXPECT_NE(message.find("#### EXITING DUE TO FAILED ASSERTION.\nFile " + std::string(__FILE__) + ", line "), std::string::npos);
    EXPECT_NE(message.find("Should never be reached.\n"), std::string::npos);
    return;
  }
  FAIL() << "Expected a not-reached exception";
}

class Embedding : public testing::Test {
protected:
  Workdir root{fs::temp_directory_path().string(), true};
  const fs::path calculation = fs::path(root.get()) / "calculation";
  const fs::path scratch = fs::path(root.get()) / "scratch";

  void SetUp() override {
    fs::create_directory(calculation);
    fs::create_directory(scratch);
  }
};

TEST_F(Embedding, RejectsInvalidDirectories) {
  const auto original_directory = fs::current_path();
  EXPECT_THROW(run_nrg(calculation.string(), ""), std::invalid_argument);
  EXPECT_THROW(run_nrg(calculation.string(), "."), std::invalid_argument);
  EXPECT_THROW(run_nrg(calculation.string(), (scratch / "missing").string()), std::invalid_argument);
  std::ofstream(scratch / "file") << "not a directory";
  EXPECT_THROW(run_nrg(calculation.string(), (scratch / "file").string()), std::invalid_argument);
  EXPECT_THROW(run_nrg("", scratch.string()), std::invalid_argument);
  EXPECT_THROW(run_nrg((calculation / "missing").string(), scratch.string()), fs::filesystem_error);
  EXPECT_THROW(run_nrg((scratch / "file").string(), scratch.string()), fs::filesystem_error);
  EXPECT_EQ(fs::current_path(), original_directory);
}

TEST_F(Embedding, CleansWorkdirWhenLoggingThrows) {
  const auto original_directory = fs::current_path();
  const auto checkpoint = fs::path(root.get()) / "checkpoint";
  fs::create_directory(checkpoint);
  std::ofstream(checkpoint / "keep") << "host checkpoint file";
  struct ThrowingOutput : std::streambuf {
    fs::path observed_directory;
    std::streamsize xsputn(const char *, std::streamsize) override {
      observed_directory = fs::current_path();
      throw std::runtime_error("injected output failure");
    }
  } output;

  for (const auto mode : {WorkdirMode::unique_temporary, WorkdirMode::persistent_exact}) {
    std::string message;
    {
      boost::io::ios_all_saver cout_state(std::cout);
      std::cout.rdbuf(&output);
      std::cout.width(0);
      std::cout.exceptions(std::ios::badbit | std::ios::failbit);
      try {
        if (mode == WorkdirMode::unique_temporary)
          run_nrg(calculation.string(), scratch.string());
        else
          Workdir workdir(checkpoint.string(), mode);
      } catch (const std::runtime_error &error) {
        message = error.what();
      }
      std::cout.clear();
    }
    EXPECT_EQ(message, "injected output failure");
    EXPECT_TRUE(fs::equivalent(output.observed_directory, mode == WorkdirMode::unique_temporary ? calculation : original_directory));
    EXPECT_EQ(fs::current_path(), original_directory);
    EXPECT_TRUE(fs::is_empty(scratch));
    EXPECT_TRUE(fs::exists(checkpoint / "keep"));
    EXPECT_NO_THROW(Workdir(checkpoint.string(), WorkdirMode::persistent_exact, true));
  }
}

TEST_F(Embedding, RestoresStateAndCleansScratchAfterFailure) {
  const auto original_directory = fs::current_path();
  struct PrecisionSaver {
    const mp_bitcnt_t precision = mpf_get_default_prec();
    ~PrecisionSaver() { mpf_set_default_prec(precision); }
  } precision_saver;
  boost::io::ios_base_all_saver cout_format(std::cout), cerr_format(std::cerr);
  boost::io::ios_fill_saver cout_fill(std::cout), cerr_fill(std::cerr);

  // Mutate host state at the workdir log flush, before missing data stops the run.
  // This exercises restoration without any numerical solver work.
  struct MutatingOutput : std::stringbuf {
    bool mutated = false;
    fs::path observed_directory;
    int sync() override {
      if (!mutated) {
        mutated = true;
        observed_directory = fs::current_path();
        mpf_set_default_prec(1024);
        for (auto *stream : {&std::cout, &std::cerr}) {
          stream->flags(std::ios::dec | std::ios::scientific);
          stream->precision(17);
          stream->width(0);
          stream->fill('!');
        }
      }
      return 0;
    }
  } output;
  boost::io::ios_rdbuf_saver cout_buffer(std::cout, &output);
  mpf_set_default_prec(192);
  const auto host_precision = mpf_get_default_prec();
  std::cout.flags(std::ios::hex | std::ios::fixed);
  std::cerr.flags(std::ios::oct | std::ios::showbase);
  std::cout.precision(5);
  std::cerr.precision(7);
  std::cout.width(11);
  std::cerr.width(13);
  std::cout.fill('@');
  std::cerr.fill('#');

  std::ofstream(fs::path(root.get()) / "DONE") << "outside calculation";
  std::ofstream(scratch / "keep") << "host scratch file";
  const auto alias = fs::path(root.get()) / "calculation-alias";
  fs::create_directory_symlink(calculation, alias);
  for (const auto &directory : {calculation, fs::relative(calculation, original_directory), alias}) {
    output.mutated = false;
    std::ofstream(calculation / "DONE") << "stale";
    try {
      run_nrg(directory.string(), scratch.string());
      FAIL() << "Expected missing data to stop the run";
    } catch (const std::runtime_error &error) {
      EXPECT_EQ(std::string(error.what()), "Can't load initial data.");
    }
    EXPECT_TRUE(output.mutated);
    EXPECT_TRUE(fs::equivalent(output.observed_directory, calculation));
    EXPECT_EQ(fs::current_path(), original_directory);
    EXPECT_EQ(mpf_get_default_prec(), host_precision);
    EXPECT_EQ(std::cout.flags(), std::ios::hex | std::ios::fixed);
    EXPECT_EQ(std::cerr.flags(), std::ios::oct | std::ios::showbase);
    EXPECT_EQ(std::cout.precision(), 5);
    EXPECT_EQ(std::cerr.precision(), 7);
    EXPECT_EQ(std::cout.width(), 11);
    EXPECT_EQ(std::cerr.width(), 13);
    EXPECT_EQ(std::cout.fill(), '@');
    EXPECT_EQ(std::cerr.fill(), '#');
    EXPECT_FALSE(fs::exists(calculation / "DONE"));
    EXPECT_TRUE(fs::exists(fs::path(root.get()) / "DONE"));
    EXPECT_TRUE(fs::exists(scratch / "keep"));
    EXPECT_EQ(std::distance(fs::directory_iterator(scratch), fs::directory_iterator{}), 1);
  }
}

TEST_F(Embedding, SerializesConcurrentCalls) {
  const auto original_directory = fs::current_path();
  constexpr size_t count = 4;
  std::array<fs::path, count> directories;
  std::array<std::string, count> errors;
  std::array<std::jthread, count> threads;
  for (size_t i = 0; i < count; ++i) {
    directories[i] = calculation / std::to_string(i);
    fs::create_directory(directories[i]);
    std::ofstream(directories[i] / "DONE") << "stale";
  }
  for (size_t i = 0; i < count; ++i) {
    threads[i] = std::jthread([&, i] {
      try {
        for (size_t repeat = 0; repeat < 3; ++repeat) {
          try {
            run_nrg(directories[i].string(), scratch.string());
            errors[i] = "Unexpected success";
          } catch (const std::runtime_error &error) {
            if (std::string(error.what()) != "Can't load initial data.") errors[i] = error.what();
          }
        }
      } catch (const std::exception &error) {
        errors[i] = error.what();
      }
    });
  }
  for (auto &thread : threads) thread.join();
  for (size_t i = 0; i < count; ++i) {
    EXPECT_TRUE(errors[i].empty()) << errors[i];
    EXPECT_FALSE(fs::exists(directories[i] / "DONE"));
  }
  EXPECT_EQ(fs::current_path(), original_directory);
  EXPECT_TRUE(fs::is_empty(scratch));
}
