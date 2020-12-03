#include <gtest/gtest.h>
#include <complex>
#include <list>
#include <map>
#include <sstream>
#include <string>
using namespace std::string_literals;
#include <fstream>
#include <exception>
#include <misc.hpp>
#include "compare.hpp"

using namespace NRG;

TEST(misc, containers) {
  {
    std::list l = {1, 2, 3, 4};
    auto b = get_back(l);
    EXPECT_EQ(b, 4);
    auto f = get_front(l);
    EXPECT_EQ(f, 1);
  }
}

TEST(misc, strings) {
  {
    auto str = "123  "s;
    auto res = strip_trailing_whitespace(str);
    EXPECT_EQ(res, "123"s);
  }
}

TEST(misc, tokenizer) {
  {
    auto str = "1 2 3 4"s;
    string_token st(str);
    EXPECT_EQ(st.find("1"), true);
    EXPECT_EQ(st.find("2"), true);
    EXPECT_EQ(st.find("3"), true);
    EXPECT_EQ(st.find("4"), true);
    EXPECT_EQ(st.find("5"), false);
  }
  {
    auto str = "ab cd ef";
    string_token st(str);
    EXPECT_EQ(st.find("ab"), true);
    EXPECT_EQ(st.find("cd"), true);
    EXPECT_EQ(st.find("ef"), true);
    EXPECT_EQ(st.find("gh"), false);
  }
}

TEST(misc, get_back) { 
  std::list a = {31,44,55,66,71,82};
  const auto b = get_back(a);   
  EXPECT_EQ(std::size(a), 5); 
  EXPECT_EQ(b, 82);          
  EXPECT_EQ(a.back(), 71);
  EXPECT_EQ(a.front(), 31);

	std::list<int> empty_list;
 	EXPECT_THROW(get_back(empty_list), std::runtime_error);

}

TEST(misc, get_front) {
  std::list a = {31,44,55,66,71,82};
  const auto b = get_front(a);
  EXPECT_EQ(std::size(a), 5);
  EXPECT_EQ(b, 31);
  EXPECT_EQ(a.back(), 82);
  EXPECT_EQ(a.front(), 44);
}

TEST(misc, switch3) {
	EXPECT_EQ(switch3(2,3,10,4,15,2,88),88);
  EXPECT_EQ(switch3(3,3,10,4,15,2,88),10);
  EXPECT_EQ(switch3(4,3,10,4,15,2,88),15);
  EXPECT_THROW(switch3(1,3,10,4,15,2,88), std::runtime_error);
}

TEST(misc, nextline) {
	auto file = safe_open_for_reading("txt/nextline.txt");
  const std::vector<std::string> result = {"zdravo", "nekaj", "adijo"};
	for(int i = 0; i < 3; i++){
		EXPECT_EQ(nextline(file),result[i]);
	}
	EXPECT_EQ(nextline(file), std::nullopt);	

	auto empty_file = safe_open_for_reading("txt/empty.txt");
	EXPECT_EQ(nextline(empty_file), std::nullopt);
}

TEST(misc, strip_trailing_whitespace) {
	const auto a = " test  \t \n  ";
	EXPECT_EQ(strip_trailing_whitespace(a), " test");

	EXPECT_EQ(strip_trailing_whitespace("  \t   \n  "s), ""s);
}



TEST(misc, block) {
	const auto nekaj = parser("txt/block.txt", "nekaj");
  const auto nekaj_drugega = parser("txt/block.txt", "nekaj drugega");
  const auto nekaj_tretjega = parser("txt/block.txt", "nekaj tretjega");
	
  const std::map nekaj_map = {std::pair("a"s, "2"s), std::pair("b"s, "3"s), std::pair("c"s, "4"s)};
  const std::map nekaj_drugega_map = {std::pair("d"s, "5"s), std::pair("c"s, "6"s), std::pair("e"s, "10"s)};
  const std::map nekaj_tretjega_map = {std::pair("abs"s, "79"s)};

	compare(nekaj, nekaj_map);
  compare(nekaj_drugega, nekaj_drugega_map);
  compare(nekaj_tretjega, nekaj_tretjega_map);
	EXPECT_THROW(parser("txt/block.txt", "zadeva"), std::runtime_error); 
}

TEST(misc, skip_comments){
  auto file = safe_open_for_reading("txt/nextline.txt");
  const std::vector result = {"zdravo"s, "nekaj"s, "adijo"s};
	std::string line;
  for(int i = 0; i < 3; i++){
		skip_comments(file);
		std::getline(file,line);
    EXPECT_EQ(line,result[i]);
  }
  skip_comments(file);
  std::getline(file,line);
  EXPECT_EQ(nextline(file),std::nullopt);
}

TEST(misc, sortfirst){
	std::vector a = {std::pair(12,5), std::pair(12,4), std::pair(99, 122), std::pair(-10, 41), std::pair(42, 21)};
  std::sort(a.begin(), a.end(), sortfirst());
	const std::vector mpkey = {-10, 12, 12, 42, 99};
	const std::vector mpval = {41,  5,  4,  21, 122};
	auto i = 0;
	for (const auto& [key, val]: a) {
		EXPECT_EQ(mpkey[i], key);
		EXPECT_EQ(mpval[i], val);
		i++;
	}
}

TEST(misc, vector_of_keys){
	const std::map a = {std::pair(12,5), std::pair(99, 122), std::pair(-10, 41), std::pair(42, 21)};
	const auto b = vector_of_keys(a);
	const std::vector expected = {-10,12,42,99};
  ASSERT_EQ(b, expected);
}






TEST(misc, ublas_to_eigen){
  auto ublas_matrix = read_matrix("txt/matrix.txt");
  auto eigen_matrix = ublas_to_eigen(ublas_matrix);
  compare(eigen_matrix, ublas_matrix);

  ublas::vector<int> ublas_vector(5);
  for(int i = 0; i < 5; i++) ublas_vector(i) = i;
  auto eigen_vector = ublas_to_eigen(ublas_vector);
  compare(eigen_vector, ublas_vector);
}

TEST(misc, eigen_to_ublas){
  Eigen::Vector4i eigen_vector(3,5,8,14);
  auto ublas_vector = eigen_to_ublas_vector(eigen_vector);
  compare(eigen_vector, ublas_vector);
  {
    Eigen::Matrix3i eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }
  
  {
    Eigen::Matrix<int, 3, 3, Eigen::RowMajor> eigen_matrix;
    for (int i = 0; i < eigen_matrix.size(); i++) eigen_matrix(i) = i;
    auto ublas_matrix = eigen_to_ublas_matrix(eigen_matrix);
    compare(eigen_matrix, ublas_matrix);
  }
}



int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
