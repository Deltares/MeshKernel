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

#include "MeshKernel/Constants.hpp"
#include "MeshKernelApi/Utils.hpp"
#include <gtest/gtest.h>

TEST(UtilsTests, MinValidElementTest_ReturnsMinimumWhenSomeValuesAreValid)
{
    std::vector<double> edgeLengths = {
        5.0,
        meshkernel::constants::missing::doubleValue,
        3.0,
        7.5};

    auto result = meshkernelapi::MinValidElement(edgeLengths);

    ASSERT_TRUE(result.has_value());
    EXPECT_DOUBLE_EQ(result.value(), 3.0);
}

TEST(UtilsTests, MinValidElementTest_ReturnsNulloptWhenAllValuesAreInvalid)
{
    std::vector<double> edgeLengths = {
        meshkernel::constants::missing::doubleValue,
        meshkernel::constants::missing::doubleValue};

    auto result = meshkernelapi::MinValidElement(edgeLengths);

    EXPECT_FALSE(result.has_value());
}
