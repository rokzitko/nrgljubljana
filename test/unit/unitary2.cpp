#include <gtest/gtest.h>
#include "unitary/unitary.hpp"
#include <string>
#include <cmake_configure.hpp>

#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <array>
#include <vector>

using namespace NRG::Unitary;
using namespace NRG;
using namespace std::string_literals;

void compare_matrices(const ublas::matrix<double> &m1, const ublas::matrix<double> &m2) {
    ASSERT_EQ(m1.size1(), m2.size1());
    ASSERT_EQ(m1.size2(), m2.size2());
    for(int i = 0; i < m1.size1(); i++)
        for(int j = 0; j < m1.size2(); j++)
            EXPECT_EQ(m1(i, j), m2(i, j));
}

template<typename T>
void compare_matrices(const std::vector<std::vector<T>> m1, const ublas::matrix<double> &m2){
    ASSERT_EQ(m2.size1(), m1.size());
    ASSERT_EQ(m2.size2(), m1[0].size());
    for(int i = 0; i < m2.size1(); i++)
        for(int j = 0; j < m2.size2(); j++)
            EXPECT_EQ(m2(i, j), m1[i][j]);
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) 
       throw std::runtime_error("popen() failed!");
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) 
       result += buffer.data();
    return result;
}

TEST(unitary, unitary_help){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -h");
    std::string  expected = "Usage: unitary [-h] [-b | -B] [-qvV] [-tl] [-s scale] [-o output_fn] [-c chop_tol] <A> <B> <C>\n";
    EXPECT_EQ(out, expected);
}

class unitaryProdTest : public ::testing::Test {
    protected:
        void SetUp() override{
            A = read_matrix("txt/matrix_A.txt");
            B = read_matrix("txt/matrix_B.txt");
            C = read_matrix("txt/matrix_C.txt");
        }

        void TearDown() override{
            std::remove("txt/temp_result_matrix.txt");
        }

    void Compare() {
        ublas::matrix<double> D = ublas::prod(B,C);
        my_result = ublas::prod(A, D);
        func_result = read_matrix("txt/temp_result_matrix.txt");
        compare_matrices(my_result, func_result);
    }
    ublas::matrix<double> A;
    ublas::matrix<double> B;
    ublas::matrix<double> C;
    ublas::matrix<double> my_result;
    ublas::matrix<double> func_result;
};

TEST_F(unitaryProdTest, noTrans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt");
    std::cout << out << std::endl;
    Compare();
}

TEST_F(unitaryProdTest, Trans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -t -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt");
    std::cout << out << std::endl;
    A = ublas::trans(A);
    Compare();
}

TEST_F(unitaryProdTest, Trans_Scale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -s 2 -t -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt");
    std::cout << out << std::endl;
    A = 2 * ublas::trans(A);
    Compare();
}
