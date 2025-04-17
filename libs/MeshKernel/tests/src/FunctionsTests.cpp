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

#include "MeshKernel/Cartesian3DPoint.hpp"

#include <chrono>
#include <gtest/gtest.h>
#include <random>

#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Operations.hpp>

TEST(FunctionsTest, NormalVectorInsideTestCartesian)
{
    // 1 Setup
    meshkernel::Point firstPoint{0.0, 0.0};
    meshkernel::Point secondPoint{1.0, 0.0};
    meshkernel::Point insidePoint{0.0, -1.0};
    meshkernel::Point normal;
    bool flippedNormal;

    // 2 Execute
    NormalVectorInside(firstPoint, secondPoint, insidePoint, normal, flippedNormal, meshkernel::Projection::cartesian);

    // 3 Validation
    ASSERT_EQ(normal.x, 0.0);
    ASSERT_EQ(normal.y, 1.0);
    ASSERT_EQ(flippedNormal, true);
}

TEST(SphericalCoordinatesTest, NorthPole)
{
    meshkernel::Point north_pole{0.0, 90.0}; // lon, lat
    auto [x, y, z] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(north_pole);

    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(x, 0.0, tolerance);
    EXPECT_NEAR(y, 0.0, tolerance);
    EXPECT_NEAR(z, meshkernel::constants::geometric::earth_radius, tolerance);
}

TEST(SphericalCoordinatesTest, SouthPole)
{
    meshkernel::Point south_pole{0.0, -90.0}; // lon, lat
    auto [x, y, z] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(south_pole);

    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(x, 0.0, tolerance);
    EXPECT_NEAR(y, 0.0, tolerance);
    EXPECT_NEAR(z, -meshkernel::constants::geometric::earth_radius, tolerance);
}

TEST(SphericalCoordinatesTest, LongitudeNegative180)
{
    meshkernel::Point point{-180.0, 0.0}; // Equator, -180 longitude
    auto [x, y, z] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(point);

    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(x, -meshkernel::constants::geometric::earth_radius, tolerance);
    EXPECT_NEAR(y, 0.0, tolerance);
    EXPECT_NEAR(z, 0.0, tolerance);
}

TEST(SphericalCoordinatesTest, LongitudePositive180)
{
    meshkernel::Point point{180.0, 0.0}; // Equator, 180 longitude
    auto [x, y, z] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(point);

    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(x, -meshkernel::constants::geometric::earth_radius, tolerance);
    EXPECT_NEAR(y, 0.0, tolerance);
    EXPECT_NEAR(z, 0.0, tolerance);
}

TEST(SphericalCoordinatesTest, LongitudeNegative180EqualsPositive180)
{
    meshkernel::Point point_neg{-180.0, 0.0};
    meshkernel::Point point_pos{180.0, 0.0};

    auto [xPos, yPos, zPos] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(point_neg);
    auto [xNeg, yNeg, zNeg] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(point_pos);

    constexpr double tolerance = 1e-6;
    EXPECT_NEAR(xPos, xNeg, tolerance);
    EXPECT_NEAR(yPos, yNeg, tolerance);
    EXPECT_NEAR(zPos, zNeg, tolerance);
}

TEST(SphericalToCartesianRoundTripTest, ConvertsAndBackCorrectly)
{
    std::vector<meshkernel::Point> testPoints = {
        {0.0, 0.0},     // Equator, Prime Meridian
        {180.0, 0.0},   // Equator, International Date Line
        {-180.0, 0.0},  // Equator, -180 deg
        {0.0, 90.0},    // North Pole
        {0.0, -90.0},   // South Pole
        {45.0, 45.0},   // NE Hemisphere
        {-90.0, 30.0},  // Western Hemisphere
        {135.0, -45.0}, // SE Hemisphere
    };

    constexpr double tolerance = 1e-6;
    for (const auto& original : testPoints)
    {
        auto [x, y, z] = ComputeSphericalCoordinatesFromLatitudeAndLongitude(original);
        meshkernel::Cartesian3DPoint cartesian{x, y, z};
        meshkernel::Point roundTripped = Cartesian3DToSpherical(cartesian, /*referenceLongitude=*/0.0);

        EXPECT_NEAR(original.x, roundTripped.x, tolerance);
        EXPECT_NEAR(original.y, roundTripped.y, tolerance);
    }
}
