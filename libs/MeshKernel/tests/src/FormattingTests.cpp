//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2025.
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation version 3.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
// contact: delft3d.support@deltares.nl
// Stichting Deltares
// P.O. Box 177
// 2600 MH Delft, The Netherlands
//
// All indications and logos of, and references to, "Delft3D" and "Deltares"
// are registered trademarks of Stichting Deltares, and remain the property of
// Stichting Deltares. All rights reserved.
//
//------------------------------------------------------------------------------

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
