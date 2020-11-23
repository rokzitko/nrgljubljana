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

TEST(unitary, count_words_in_string) {
    EXPECT_EQ(count_words_in_string("one two three four a"s), 5);
    EXPECT_EQ(count_words_in_string(""s), 0);
}

TEST(unitary, get_dims) {
    auto file = safe_open_for_reading("matrix.txt");
    auto const [dim1, dim2] = get_dims(file);
    EXPECT_EQ(dim1, 3);
    EXPECT_EQ(dim2, 4);
    auto file_err = safe_open_for_reading("matrix_err.txt");
    EXPECT_THROW(get_dims(file_err), std::runtime_error);
}

using MAT = ublas::matrix<double>;

void compare_matrices(const MAT &m1, const MAT &m2) {
    ASSERT_EQ(m1.size1(), m2.size1());
    ASSERT_EQ(m1.size2(), m2.size2());
    for(int i = 0; i < m1.size1(); i++)
        for(int j = 0; j < m1.size2(); j++)
            EXPECT_EQ(m1(i, j), m2(i, j));
}

template<typename T>
void compare_matrices(const std::vector<std::vector<T>> m1, const MAT &m2){
    ASSERT_EQ(m2.size1(), m1.size());
    ASSERT_EQ(m2.size2(), m1[0].size());
    for(int i = 0; i < m2.size1(); i++)
        for(int j = 0; j < m2.size2(); j++)
            EXPECT_EQ(m2(i, j), m1[i][j]);
}

TEST(unitary, read_matrix_text){
    std::vector<std::vector<int>> const ref_matrix = {{31,41,53,46},{12,5,1,41},{5,2,4,7}};
    auto matrix = read_matrix("matrix.txt");
    compare_matrices(ref_matrix, matrix);
}

TEST(unitary, save_matrix){
    auto matrix = read_matrix("matrix.txt");
    save_matrix("matrix_temp.txt", matrix);
    auto matrix_temp = read_matrix("matrix_temp.txt");
    compare_matrices(matrix, matrix_temp);
    std::remove("matrix_temp.txt");
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
            A = read_matrix("matrix_A.txt");
            B = read_matrix("matrix_B.txt");
            C = read_matrix("matrix_C.txt");
        }

        void TearDown() override{
            std::remove("temp_result_matrix.txt");
        }

    void Compare() {
        MAT D = ublas::prod(B,C);
        my_result = ublas::prod(A, D);
        func_result = read_matrix("temp_result_matrix.txt");
        compare_matrices(my_result, func_result);
    }
    MAT A;
    MAT B;
    MAT C;
    MAT my_result;
    MAT func_result;
};

TEST_F(unitaryProdTest, noTrans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -q -o temp_result_matrix.txt matrix_A.txt matrix_B.txt matrix_C.txt");
    std::cout << out << std::endl;
    Compare();
}

TEST_F(unitaryProdTest, Trans_noScale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -t -q -o temp_result_matrix.txt matrix_A.txt matrix_B.txt matrix_C.txt");
    std::cout << out << std::endl;
    A = ublas::trans(A);
    Compare();
}

TEST_F(unitaryProdTest, Trans_Scale){
    auto out = exec(PROJECT_BINARY_DIR "/tools/unitary -s 2 -t -q -o temp_result_matrix.txt matrix_A.txt matrix_B.txt matrix_C.txt");
    std::cout << out << std::endl;
    A = 2 * ublas::trans(A);
    Compare();
}
