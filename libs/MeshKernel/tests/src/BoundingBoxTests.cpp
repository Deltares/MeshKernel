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

#include <gtest/gtest.h>

#include "MeshKernel/BoundingBox.hpp"
#include "MeshKernel/Entities.hpp"

TEST(BoundingBox, DefaultIntialization_MustIntializeCornersToNumericLimits)
{
    // Setup
    const auto boundingBox = meshkernel::BoundingBox();

    // Assert
    ASSERT_EQ(boundingBox.lowerLeft().x, std::numeric_limits<double>::lowest());
    ASSERT_EQ(boundingBox.lowerLeft().y, std::numeric_limits<double>::lowest());
    ASSERT_EQ(boundingBox.upperRight().x, std::numeric_limits<double>::max());
    ASSERT_EQ(boundingBox.upperRight().y, std::numeric_limits<double>::max());
}

TEST(BoundingBox, Contains_WhenPointInside_MustReturnTrue)
{
    // Setup
    const auto boundingBox = meshkernel::BoundingBox({0.0, 0.0}, {10.0, 10.0});

    std::vector<meshkernel::Point> points;
    points.push_back({0.0, 0.0});
    points.push_back({5.5, 5.0});
    points.push_back({10.0, 10.0});
    points.push_back({10.1, 10.1});

    // Assert
    ASSERT_TRUE(boundingBox.Contains(points[0]));
    ASSERT_TRUE(boundingBox.Contains(points[1]));
    ASSERT_TRUE(boundingBox.Contains(points[2]));
    ASSERT_FALSE(boundingBox.Contains(points[3]));
}

TEST(BoundingBox, Contains_WhenPointOutside_MustReturnFalse)
{
    // 1 Setup
    const auto boundingBox = meshkernel::BoundingBox({0.0, 0.0}, {10.0, 10.0});

    std::vector<meshkernel::Point> points;
    points.push_back({-10.0, -10.0});
    points.push_back({20.0, 20.0});

    // Assert
    ASSERT_FALSE(boundingBox.Contains(points[0]));
    ASSERT_FALSE(boundingBox.Contains(points[1]));
}
