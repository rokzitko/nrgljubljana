#include <gtest/gtest.h>
#include <portabil.hpp>

#include <sstream>

TEST(parse_string, stream) {
   std::stringstream ss;
   ss << "keyword: value\n" << "keyword2: value2\n";
   ss << "Xkeyword3:Xvalue3Y\n";
     {
	auto res = parse_string(ss, "keyword"); // find the first occurrence!
	EXPECT_TRUE(res);
	EXPECT_EQ(res.value(), ": value");
     }
     {
	auto res = parse_string(ss, "keyword2: ");
	EXPECT_TRUE(res);
	EXPECT_EQ(res.value(), "value2");
     }
     {
	auto res = parse_string(ss, "keyword3");
	EXPECT_TRUE(res);
	EXPECT_EQ(res.value(), "X:Xvalue3Y");
     }
     {
	auto res = parse_string(ss, "not_present");
	EXPECT_FALSE(res);
     }
}

int main(int argc, char **argv) {
   ::testing::InitGoogleTest(&argc, argv);
   return RUN_ALL_TESTS();
}
