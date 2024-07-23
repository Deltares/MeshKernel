//---- GPL ---------------------------------------------------------------------
//
// Copyright (C)  Stichting Deltares, 2011-2024.
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

#include <cmath>
#include <numbers>
#include <random>
#include <utility>

#include <gtest/gtest.h>

#include "MeshKernel/CurvilinearGrid/CurvilinearGrid.hpp"
#include "MeshKernel/CurvilinearGrid/CurvilinearGridGenerateCircularGrid.hpp"
#include "MeshKernel/Entities.hpp"
#include "MeshKernel/Operations.hpp"
#include "MeshKernel/Parameters.hpp"
#include "MeshKernel/UndoActions/UndoActionStack.hpp"
#include "MeshKernel/Utilities/LinearAlgebra.hpp"
#include "TestUtils/MakeCurvilinearGrids.hpp"

namespace mk = meshkernel;

TEST(CurvilinearGridGenerationTests, FullyUniformGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 12,
        .angle = 0.0,
        .block_size_x = 10.0,
        .block_size_y = 20.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 1.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;
    double y = 0.0;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        double x = 0.0;

        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            EXPECT_NEAR(x, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(y, grid.GetNode(i, j).y, tolerance);
            x += parameters.block_size_x;
        }

        y += parameters.block_size_y;
    }
}

TEST(CurvilinearGridGenerationTests, FullyUniformGridRotated)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 19,
        .num_rows = 7,
        .angle = 23.0, // degrees
        .block_size_x = 7.0,
        .block_size_y = 15.0,
        .radius_curvature = 0.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 1.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    constexpr double tolerance = 1.0e-10;

    double cs = std::cos(parameters.angle * mk::constants::conversion::degToRad);
    double sn = std::sin(parameters.angle * mk::constants::conversion::degToRad);

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    double y = 0.0;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        double x = 0.0;

        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            double expectedX = x * cs - y * sn;
            double expectedY = x * sn + y * cs;
            EXPECT_NEAR(expectedX, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(expectedY, grid.GetNode(i, j).y, tolerance);
            x += parameters.block_size_x;
        }

        y += parameters.block_size_y;
    }
}

TEST(CurvilinearGridGenerationTests, UniformXGradedYGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 12,
        .angle = 0.0,
        .block_size_x = 10.0,
        .block_size_y = 20.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 0.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    std::vector<double> gradedYValues{0.0,
                                      31.6763921753316,
                                      71.5410984881053,
                                      121.710789550376,
                                      184.849293106268,
                                      264.309033576453,
                                      364.309033576453,
                                      490.158928640636,
                                      648.540889517293,
                                      847.864421081162,
                                      1098.71287639251,
                                      1414.40539417197,
                                      1811.7040965229};

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        double x = 0.0;

        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            EXPECT_NEAR(x, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(gradedYValues[i], grid.GetNode(i, j).y, tolerance);
            x += parameters.block_size_x;
        }
    }
}

TEST(CurvilinearGridGenerationTests, GradedXUniformYGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 12,
        .angle = 0.0,
        .block_size_x = 10.0,
        .block_size_y = 20.0,
        .uniform_columns_fraction = 0.0,
        .uniform_rows_fraction = 1.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    std::vector<double> gradedXValues{0.0,
                                      38.2362245665865,
                                      67.4764019487152,
                                      89.8370817237131,
                                      106.93684119048,
                                      120.013446050598,
                                      133.090050910717,
                                      150.189810377484,
                                      172.550490152481,
                                      201.79066753461,
                                      240.026892101197};

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;
    double y = 0.0;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {

        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            EXPECT_NEAR(gradedXValues[j], grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(y, grid.GetNode(i, j).y, tolerance);
        }

        y += parameters.block_size_y;
    }
}

TEST(CurvilinearGridGenerationTests, UniformCircularGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 10,
        .angle = 0.0,
        .radius_curvature = 10.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 1.0,
        .maximum_uniform_size_columns = 10.0,
        .maximum_uniform_size_rows = 10.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    std::vector<double> radiusValues{10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0, 110.0};
    std::vector<double> thetaValues{1.5707963267949,
                                    0.942477796076938,
                                    0.314159265358979,
                                    -0.314159265358979,
                                    -0.942477796076938,
                                    -1.5707963267949,
                                    -2.19911485751286,
                                    -2.82743338823081,
                                    -3.45575191894877,
                                    -4.08407044966673,
                                    -4.71238898038469};

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            double expectedX = parameters.origin_x + radiusValues[i] * std::cos(thetaValues[j]);
            double expectedY = parameters.origin_y + radiusValues[i] * std::sin(thetaValues[j]);
            EXPECT_NEAR(expectedX, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(expectedY, grid.GetNode(i, j).y, tolerance);
        }
    }
}

TEST(CurvilinearGridGenerationTests, HalfUniformCircularGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 10,
        .num_rows = 10,
        .angle = 0.0,
        .radius_curvature = 10.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 0.5,
        .maximum_uniform_size_columns = 10.0,
        .maximum_uniform_size_rows = 10.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    std::vector<double> radiusValues{20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 101.544346900319, 147.960235236447, 247.960235236447, 463.403704239635};
    std::vector<double> thetaValues{1.5707963267949,
                                    0.942477796076938,
                                    0.314159265358979,
                                    -0.314159265358979,
                                    -0.942477796076938,
                                    -1.5707963267949,
                                    -2.19911485751286,
                                    -2.82743338823081,
                                    -3.45575191894877,
                                    -4.08407044966673,
                                    -4.71238898038469};

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            double expectedX = parameters.origin_x + radiusValues[i] * std::cos(thetaValues[j]);
            double expectedY = parameters.origin_y + radiusValues[i] * std::sin(thetaValues[j]);
            EXPECT_NEAR(expectedX, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(expectedY, grid.GetNode(i, j).y, tolerance);
        }
    }
}

TEST(CurvilinearGridGenerationTests, ShiftedRotatedCircularGrid)
{

    const mk::MakeGridParameters parameters = {
        .num_columns = 14,
        .num_rows = 10,
        .angle = 32.0,
        .origin_x = 17.0,
        .origin_y = -23.0,
        .radius_curvature = 10.0,
        .uniform_columns_fraction = 1.0,
        .uniform_rows_fraction = 0.0};

    auto grid = mk::CurvilinearGridGenerateCircularGrid::Compute(parameters, mk::Projection::cartesian);

    std::vector<double> radiusValues{20.0, 33.0766048601183, 50.1763643268853,
                                     72.5370441018832, 101.777221484012, 140.013446050598,
                                     190.013446050598, 255.39647035119, 340.895267685025,
                                     452.698666560014, 598.899553470658};

    std::vector<double> thetaValues{2.12930168743308, 1.68050273692025, 1.23170378640743, 0.782904835894599,
                                    0.334105885381772, -0.114693065131056, -0.563492015643884, -1.01229096615671,
                                    -1.46108991666954, -1.90988886718237, -2.35868781769519, -2.80748676820802,
                                    -3.25628571872085, -3.70508466923368, -4.1538836197465};

    ASSERT_EQ(grid.NumM(), static_cast<mk::UInt>(parameters.num_columns + 1));
    ASSERT_EQ(grid.NumN(), static_cast<mk::UInt>(parameters.num_rows + 1));

    constexpr double tolerance = 1.0e-10;

    for (mk::UInt i = 0; i < grid.NumN(); ++i)
    {
        for (mk::UInt j = 0; j < grid.NumM(); ++j)
        {
            double expectedX = parameters.origin_x + radiusValues[i] * std::cos(thetaValues[j]);
            double expectedY = parameters.origin_y + radiusValues[i] * std::sin(thetaValues[j]);
            EXPECT_NEAR(expectedX, grid.GetNode(i, j).x, tolerance);
            EXPECT_NEAR(expectedY, grid.GetNode(i, j).y, tolerance);
        }
    }
}
