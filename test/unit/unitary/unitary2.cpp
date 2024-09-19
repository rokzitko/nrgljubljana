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

#include <traits.hpp>
#include <numerics.hpp>

#include "compare.hpp"
using namespace NRG::Unitary;
using namespace NRG;
using namespace std::string_literals;


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
        EigenMatrix<double> D = B*C;
        my_result = A*D;
        func_result = read_matrix("txt/temp_result_matrix.txt");
        compare(my_result, func_result);
    }
    EigenMatrix<double> A;
    EigenMatrix<double> B;
    EigenMatrix<double> C;
    EigenMatrix<double> my_result;
    EigenMatrix<double> func_result;
};



TEST_F(unitaryProdTest, noTrans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt");
    std::cout << out << std::endl;
    Compare();
}

TEST_F(unitaryProdTest, Trans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -t -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt"); // t = transpose first (A)
    std::cout << out << std::endl;
    A.transposeInPlace();
    Compare();
}

TEST_F(unitaryProdTest, Trans_Scale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -s 2 -t -q -o txt/temp_result_matrix.txt txt/matrix_A.txt txt/matrix_B.txt txt/matrix_C.txt");
    std::cout << out << std::endl;
    A.transposeInPlace();
    A = 2.0*A;
    Compare();
}
