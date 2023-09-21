#include <string>
#include <vector>

#include <gtest/gtest.h>

#include "MeshKernel/Formatting.hpp"

TEST(Formatting, FormattedVectorOfChars)
{
    std::string const formatted_str =
        fmt_ns::format("{}", std::vector<char>{'G', 'i', 'u', 's', 'e', 'p', 'p', 'e'});
    std::string const expected_str = "{G, i, u, s, e, p, p, e}";
    EXPECT_EQ(formatted_str, expected_str);
}

TEST(Formatting, FormattedVectorOfStrings)
{
    std::string const formatted_str =
        fmt_ns::format("{}", std::vector<std::string>{"Pomodori", "della", "nonna"});
    std::string const expected_str = "{Pomodori, della, nonna}";
    EXPECT_EQ(formatted_str, expected_str);
}

TEST(Formatting, FormattedVectorOfIntegers)
{
    std::string const formatted_str =
        fmt_ns::format("{}", std::vector<int>{-1, -2, 0, 1, 2});
    std::string const expected_str = "{-1, -2, 0, 1, 2}";
    EXPECT_EQ(formatted_str, expected_str);
}

TEST(Formatting, FormattedVectorOfFloats)
{
    std::string const formatted_str =
        fmt_ns::format("{:.5f}", std::vector<float>{-1.5f, -2.2f, 0.0f, 1.0f, 2.75f});
    std::string const expected_str = "{-1.50000, -2.20000, 0.00000, 1.00000, 2.75000}";
    EXPECT_EQ(formatted_str, expected_str);
}
