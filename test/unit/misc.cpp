#include <gtest/gtest.h>
#include <misc.hpp>
#include <complex>
#include <list>
#include <sstream>

using namespace NRG;

TEST(containers, misc) {
  {
    std::list<> l = {1, 2, 3, 4};
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
    auto st = string_tokenizer(str);
    EXPECT_EQ(st.find("1"), true);
    EXPECT_EQ(st.find("2"), true);
    EXPECT_EQ(st.find("3"), true);
    EXPECT_EQ(st.find("4"), true);
    EXPECT_EQ(st.find("5"), false);
  }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
