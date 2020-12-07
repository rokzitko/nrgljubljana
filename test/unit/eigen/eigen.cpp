#include <string>
#include <sstream>
#include <gtest/gtest.h>

#include "test_common.hpp"
#include <traits.hpp>
#include <eigen.hpp>

using namespace NRG;

TEST(RawEigen, constructor) { // NOLINT
  {
    NRG::RawEigen<double> e;
    EXPECT_EQ(e.getnrcomputed(), 0);
    EXPECT_EQ(e.getdim(), 0);
  }
  {
    NRG::RawEigen<double> e(2,3);
    EXPECT_EQ(e.getnrcomputed(), 2);
    EXPECT_EQ(e.getdim(), 3);
  }
}

TEST(Eigen, get) { // NOLINT
  NRG::Eigen<double> e(5,8);
  e.truncate_prepare(2);
  EXPECT_EQ(e.getnrcomputed(), 5);
  EXPECT_EQ(e.getdim(), 8);
  EXPECT_EQ(e.getnrall(), 5);
  EXPECT_EQ(e.getnrkept(), 2);
  EXPECT_EQ(e.getnrdiscarded(), 3);
  EXPECT_EQ(e.all(), range0(5));
  EXPECT_EQ(e.kept(), range0(2));
  EXPECT_EQ(e.discarded(), boost::irange(2,5) );
  EXPECT_EQ(e.getnrstored(), 5);
  EXPECT_EQ(e.stored(), range0(5));
}

using EVEC = evec_traits<double>;

TEST(Eigen, diagonal) { // NOLINT
  NRG::Eigen<double> e(3,3);
  EVEC v = { 1.0, 2.0, 3.0 };
  e.diagonal(v);
  EXPECT_EQ(e.values.rel(0), 1.0);
  EXPECT_EQ(e.value_zero[0], 1.0);
  EXPECT_EQ(e.matrix(0,0), 1.0); // identity matrix
  EXPECT_EQ(e.matrix(0,1), 0.0);
  EXPECT_EQ(e.matrix(1,1), 1.0);
  EXPECT_EQ(e.matrix(2,2), 1.0);
  e.subtract_Egs(1.0);
  EXPECT_EQ(e.value_zero[0], 0.0);
}

TEST(Eigen, io) { // NOLINT
  NRG::Eigen<double> e(2,2);
  EVEC v1 = { 1.0, 2.0 };
  EVEC v2 = { 3.0, 4.0 };
  e.diagonal(v1);
  EXPECT_EQ(e.values.rel(0), 1.0);
  std::ostringstream oss;
  boost::archive::binary_oarchive oa(oss);
  e.save(oa);
  e.diagonal(v2);
  EXPECT_EQ(e.values.rel(0), 3.0);
  std::istringstream iss(oss.str());
  boost::archive::binary_iarchive ia(iss);
  e.load(ia);
  EXPECT_EQ(e.values.rel(0), 1.0); // original value
}

TEST(Eigen, hdf5io) { // NOLINT
  NRG::Eigen<double> e(2,2);
  EVEC v1 = { 1.0, 2.0 };
  e.diagonal(v1);
  auto h5 = H5Easy::File("Eigen.h5", H5Easy::File::Overwrite);
  e.h5save(h5, "test", false);
}

template<typename U, typename V> // XXX: concept vector: .size, [] accessor
void VECTOR_EQ(const U &A, const V &B)
{
  EXPECT_EQ(A.size(), B.size());
  for(int i = 0; i < A.size(); i++) EXPECT_EQ(A[i], B[i]);
}

template<typename U, typename V> // XXX: concept matrix: size1(), size2(), (i,j) accessor
void MATRIX_EQ(const U &A, const V &B)
{
  EXPECT_EQ(A.size1(), B.size1());
  EXPECT_EQ(A.size2(), B.size2());
  for(int i = 0; i < A.size1(); i++) 
    for(int j = 0; j < A.size2(); j++) 
      EXPECT_EQ(A(i,j), B(i,j));
}

TEST(io, read_vector) { // NOLINT
  std::string data = "5 1 2 3 4 5";
  std::istringstream ss(data);
  auto vec = read_vector<double>(ss);
  EVEC ref(5);
  ref[0] = 1; ref[1] = 2; ref[2] = 3; ref[3] = 4; ref[4] = 5;
  VECTOR_EQ(vec, ref);
}

TEST(io, read_matrix) { // NOLINT
  std::string data = "1 2 3\n 4 5 6\n";
  std::istringstream ss(data);
  auto mat = read_matrix<double>(ss, 2, 3);
  ublas::matrix<double> ref(2,3);
  ref(0,0) = 1; ref(0,1) = 2; ref(0,2) = 3; ref(1,0) = 4; ref(1,1) = 5; ref(1,2) = 6;
  MATRIX_EQ(mat, ref);
}

template<typename T>
auto range_size(T t)
{
  size_t ctr = 0;
  for(const auto &x: t) ctr++;
  return ctr;
}

TEST(Diag, constructor) { // NOLINT
  std::string data =
    "0 1\n"
    "2 1 2\n"
    "1 2\n"
    "3 4 5 6\n";
  std::istringstream ss(data);
  Params P;
  auto Sym = setup_Sym<double>(P); // need working Invar
  P.absolute.setvalue(true); // disable any energy rescaling
  DiagInfo<double> diag(ss, 2, P);
  EXPECT_EQ(diag.size(), 2);
  EXPECT_EQ(range_size(diag.subspaces()), 2);
  EXPECT_EQ(range_size(diag.eigs()), 2);
  EXPECT_EQ(diag.find_groundstate(), 1.0);
  std::vector ref_energies = { 1.0, 2.0, 4.0, 5.0, 6.0 };
  EXPECT_EQ(diag.sorted_energies_rel(), ref_energies);
  EXPECT_EQ(diag.size_subspace(Invar(0,1)), 2);
  EXPECT_EQ(diag.size_subspace(Invar(1,2)), 3);
  EXPECT_EQ(diag.size_subspace(Invar(0,0)), 0);
  EXPECT_EQ(diag.dims(Invar(0,1),Invar(1,2)), std::make_pair(2ul,3ul));
  EXPECT_EQ(diag.count_states(Sym->multfnc()), 2*1+3*2); // multiplicity!
  EXPECT_EQ(diag.trace([](const auto x){ return 1; }, 1.0, Sym->multfnc()), 
                       exp(-1.0)+exp(-2.0)+2*exp(-4.0)+2*exp(-5.0)+2*exp(-6.0));
//  diag.save(3, P);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS(); // NOLINT
}
