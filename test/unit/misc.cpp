#include <gtest/gtest.h>
#include <complex>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <fstream>

#include <misc.hpp>

using namespace NRG;

TEST(containers, misc) {
  {
    std::list l = {1, 2, 3, 4};
    auto b = get_back(l);
    EXPECT_EQ(b, 4);
    auto f = get_front(l);
    EXPECT_EQ(f, 1);
  }
}

TEST(strings, misc) {
  {
    auto str = "123  ";
    auto res = strip_trailing_whitespace(str);
    EXPECT_EQ(res, "123");
  }
}

TEST(tokenizer, misc) {
  {
    auto str = "1 2 3 4";
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
}

TEST(misc, nextline) {
	auto file = safe_open_for_reading("nextline.txt");
  const std::vector<std::string> result = {"zdravo", "nekaj", "adijo"};
	for(int i = 0; i < 3; i++){
		EXPECT_EQ(nextline(file),result[i]);
	}
	EXPECT_EQ(nextline(file), std::nullopt);	
}

TEST(misc, strip_trailing_whitespace) {
	const auto a = " test  \t \n  ";
	EXPECT_EQ(strip_trailing_whitespace(a), " test");
}

template<class T1, class T2>
void compare_maps(std::map<T1,T2> a, std::map<T1,T2> b) {
  ASSERT_EQ(a.size(), b.size());
  for(auto const& [key, value] : a)
    EXPECT_EQ(b[key], value);
}

TEST(misc, block) {
	const auto nekaj = parser("block.txt", "nekaj");
  const auto nekaj_drugega = parser("block.txt", "nekaj drugega");
  const auto nekaj_tretjega = parser("block.txt", "nekaj tretjega");
	
  std::map nekaj_map = {std::pair<std::string,std::string>("a", "2"), std::pair<std::string,std::string>("b", "3"), std::pair<std::string,std::string>("c", "4")};
  std::map nekaj_drugega_map = {std::pair<std::string,std::string>("d", "5"), std::pair<std::string,std::string>("c", "6"), std::pair<std::string,std::string>("e", "10")};
  std::map nekaj_tretjega_map = {std::pair<std::string,std::string>("abs", "79")};

	compare_maps(nekaj, nekaj_map);
  compare_maps(nekaj_drugega, nekaj_drugega_map);
  compare_maps(nekaj_tretjega, nekaj_tretjega_map);
}

TEST(misc, skip_comments){
  auto file = safe_open_for_reading("nextline.txt");
  const std::vector<std::string> result = {"zdravo", "nekaj", "adijo"};
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

/*TEST(misc, sortfirst){

	std::map<int,int> sortf = {std::pair<int,int>(12,5), std::pair<int,int>(99, 122), std::pair<int,int>(-10, 41), std::pair<int,int>(42, 21)};
	sortf.value_compare = sortfirst;
	int mpkey[4] = {-10, 12, 42, 99};
	int mpval[4] = {41, 5, 21, 122};
	int i = 0;
	for (const auto& [key, val]: sortf) {
		EXPECT_EQ(mpkey[i], key);
		EXPECT_EQ(mpval[i], val);
		i++;
	}
}*/

TEST(misc, vector_of_keys){
	std::map<int,int> a = {std::pair<int,int>(12,5), std::pair<int,int>(99, 122), std::pair<int,int>(-10, 41), std::pair<int,int>(42, 21)};
	std::vector<int> b = vector_of_keys(a);
	std::vector<int> test = {-10,12,42,99};
	ASSERT_EQ(b.size(), test.size());
	for(int i = 0; i < b.size(); i++)
		EXPECT_EQ(b[i], test[i]);
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
