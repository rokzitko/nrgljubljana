#include <string>
#include <sstream>
#include <exception>
#include <filesystem>
#include <fstream>
#include <gtest/gtest.h>
#include <ios>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdexcept>

#include "test_common.hpp"
#include "test_clean.hpp"
#include <operators.hpp>

using namespace NRG;

TEST(Operators, MatrixElements) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag<double>(P, Sym);
  {
    std::string str =
      "1\n"
      "0 1 0 1\n"
      "1 2\n 3 4\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "2\n"
      "0 1 0 1\n"
      "1 2\n 3 4\n"
      "1 2 1 2\n"
      "1 2 3\n 4 5 6\n 7 8 9\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "1\n"
      "0 1 1 2\n"
      "1 2 3\n 4 5 6\n";
    std::istringstream ss(str);
    auto me = MatrixElements(ss, diag);
    std::cout << me << std::endl;
  }
  {
    std::string str =
      "1\n"
      "0 1 0 99\n"
      "1 2\n 3 4\n";
    std::istringstream ss(str);
    EXPECT_THROW(MatrixElements(ss, diag), std::runtime_error);
  }
}

TEST(Operators, Opch_empty) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag<double>(P, Sym);
  {
    auto o = Opch<double>(P);
    o.dump();
    EXPECT_EQ(o.size(), 1); // channels
    EXPECT_EQ(o[0].size(), 1); // perchannel
    EXPECT_EQ(o[0][0].size(), 0);
  }
}

TEST(Operators, Opch_clean) { // NOLINT
  Params P;
  auto SymSP = setup_Sym<double>(P);
  auto Sym = SymSP.get();
  auto diag = setup_diag_clean<double>(P, Sym);
  diag.states_report(Sym->multfnc());
  std::string str =
    "f 0 0\n"
    "2\n"
    "1 1 0 2\n"
    "1.4142135623730951\n"
    "0 2 -1 1\n"
    "1.\n";
  std::istringstream ss(str);
  auto o = Opch<double>(ss, diag, P);
  o.dump();
}

TEST(Operators, DensMatSaveUsesTemporaryFile) { // NOLINT
  Params P;
  DensMatElements<double> rho;
  rho[Invar()] = zero_matrix<double>(1);

  const auto fn = P.workdir->rhofn(0, fn_rho);
  rho.save(0, P, fn_rho);

  EXPECT_TRUE(std::filesystem::exists(fn));
  EXPECT_FALSE(std::filesystem::exists(fn + ".tmp"));
}

TEST(Operators, DensMatLoadPreservesCorruptFileOnFailure) { // NOLINT
  Params P;
  const auto fn = P.workdir->rhofn(0, fn_rho);
  std::ofstream(fn, std::ios::binary | std::ios::out) << "bad";

  DensMatElements<double> rho;
  EXPECT_THROW(rho.load(0, P, fn_rho, true), std::exception);
  EXPECT_TRUE(std::filesystem::exists(fn));
}

TEST(Operators, DensMatLoadDoesNotMutateExistingDataOnFailure) { // NOLINT
  Params P;
  const auto fn = P.workdir->rhofn(0, fn_rho);
  std::ofstream(fn, std::ios::binary | std::ios::out) << "bad";

  DensMatElements<double> rho;
  rho[Invar()] = zero_matrix<double>(1);
  rho[Invar()](0,0) = 7.0;

  EXPECT_THROW(rho.load(0, P, fn_rho, false), std::exception);
  ASSERT_EQ(rho.size(), 1U);
  EXPECT_DOUBLE_EQ(rho.at(Invar())(0,0), 7.0);
  EXPECT_TRUE(std::filesystem::exists(fn));
}

TEST(Operators, DensMatLoadRemovesFileAfterSuccessWhenRequested) { // NOLINT
  Params P;
  DensMatElements<double> rho_out;
  rho_out[Invar()] = zero_matrix<double>(1);
  rho_out[Invar()](0,0) = 5.0;
  const auto fn = P.workdir->rhofn(0, fn_rho);
  rho_out.save(0, P, fn_rho);

  DensMatElements<double> rho_in;
  rho_in.load(0, P, fn_rho, true);

  EXPECT_FALSE(std::filesystem::exists(fn));
  ASSERT_EQ(rho_in.size(), 1U);
  EXPECT_DOUBLE_EQ(rho_in.at(Invar())(0,0), 5.0);
}

TEST(Operators, DensMatRejectsIncompatibleCheckpoint) { // NOLINT
  Params writer;
  writer.Nmax = 2;
  writer.checkpoint_data_fingerprint = 123;
  DensMatElements<double> rho_out;
  rho_out[Invar()] = zero_matrix<double>(1);
  rho_out.save(0, writer, fn_rho);

  Params reader;
  reader.Nmax = 2;
  reader.workdir = std::make_unique<Workdir>(writer.workdir->get(), WorkdirMode::persistent_exact, true);

  EXPECT_FALSE(DensMatElements<double>::compatible_file(0, reader, fn_rho));
  EXPECT_THROW(DensMatElements<double>().load(0, reader, fn_rho, false), std::runtime_error);
  EXPECT_TRUE(std::filesystem::exists(writer.workdir->rhofn(0, fn_rho)));
}

TEST(Operators, DensMatRejectsDifferentUnitaryChain) { // NOLINT
  Params writer;
  writer.Nmax = 2;
  writer.checkpoint_unitary_fingerprint = 123;
  DensMatElements<double> rho_out;
  rho_out[Invar()] = zero_matrix<double>(1);
  rho_out.save(0, writer, fn_rho);

  Params reader;
  reader.Nmax = 2;
  reader.workdir = std::make_unique<Workdir>(writer.workdir->get(), WorkdirMode::persistent_exact, true);

  EXPECT_FALSE(DensMatElements<double>::compatible_file(0, reader, fn_rho));
}

TEST(Operators, DensMatCompatibilityReadsEntirePayload) { // NOLINT
  Params P;
  DensMatElements<double> rho_out;
  rho_out[Invar()] = zero_matrix<double>(2);
  rho_out.save(0, P, fn_rho);
  const auto fn = P.workdir->rhofn(0, fn_rho);
  const auto size = std::filesystem::file_size(fn);
  ASSERT_GT(size, 1U);
  std::filesystem::resize_file(fn, size - 1);

  EXPECT_FALSE(DensMatElements<double>::compatible_file(0, P, fn_rho));
  EXPECT_THROW(DensMatElements<double>().load(0, P, fn_rho, false), std::exception);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
