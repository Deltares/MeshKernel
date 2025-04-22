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

#include <gmock/gmock.h>
#include <gtest/gtest.h>

#include <vector>

#include "MeshKernel/SampleAveragingInterpolator.hpp"
#include "MeshKernel/SampleTriangulationInterpolator.hpp"
#include "TestUtils/Definitions.hpp"
#include "TestUtils/MakeMeshes.hpp"
#include "TestUtils/SampleFileReader.hpp"

namespace mk = meshkernel;

TEST(SampleInterpolationTests, AveragingInterpolationWithPoints)
{
    const mk::UInt numberOfPointsX = 11;
    const mk::UInt numberOfPointsY = 11;
    const mk::UInt numberOfPoints = numberOfPointsX * numberOfPointsY;

    std::vector<double> xPoints(numberOfPoints);
    std::vector<double> yPoints(numberOfPoints);
    std::vector<double> data(numberOfPoints);

    // Generate sample data points
    double delta = 1000.0;
    double x = 0.0;
    double y = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < numberOfPointsY; ++i)
    {
        x = 0.0;

        for (size_t j = 0; j < numberOfPointsX; ++j)
        {
            xPoints[count] = x;
            yPoints[count] = y;
            data[count] = x;
            x += delta;
            ++count;
        }

        y += delta;
    }

    mk::InterpolationParameters params{.absolute_search_radius = 3.0 * delta};
    mk::SampleAveragingInterpolator interpolator(xPoints, yPoints, mk::Projection::cartesian, params);
    ASSERT_EQ(static_cast<size_t>(interpolator.Size()), xPoints.size());

    int propertyId = 1;
    interpolator.SetData(propertyId, data);

    // Execute
    ASSERT_EQ(interpolator.Size(), numberOfPoints);

    std::vector<mk::Point> interpolationPoints{{0.5 * delta, 0.5 * delta}, {9.5 * delta, 9.5 * delta}};
    const double initialValue = -1.0e20;
    std::vector<double> interpolationResult(interpolationPoints.size(), initialValue);
    std::vector<double> expectedResult{1400.0, 8600.0};

    interpolator.Interpolate(propertyId, interpolationPoints, interpolationResult);

    const double tolerance = 1.0e-8;

    for (size_t i = 0; i < expectedResult.size(); ++i)
    {
        EXPECT_NEAR(expectedResult[i], interpolationResult[i], tolerance);
    }
}

TEST(SampleInterpolationTests, AveragingInterpolationWithMesh)
{
    const mk::UInt numberOfPointsX = 11;
    const mk::UInt numberOfPointsY = 11;
    const mk::UInt numberOfPoints = numberOfPointsX * numberOfPointsY;

    std::vector<mk::Point> samplePoints(numberOfPoints);
    std::vector<double> data(numberOfPoints);

    // Generate sample data points
    double delta = 1000.0;
    double x = 0.0;
    double y = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < numberOfPointsY; ++i)
    {
        x = 0.0;

        for (size_t j = 0; j < numberOfPointsX; ++j)
        {
            samplePoints[i].x = x;
            samplePoints[i].y = y;
            data[count] = x;
            x += delta;
            ++count;
        }

        y += delta;
    }

    mk::InterpolationParameters params{.absolute_search_radius = 3.0 * delta, .minimum_number_of_samples = 1};
    mk::SampleAveragingInterpolator interpolator(samplePoints, mk::Projection::cartesian, params);

    ASSERT_EQ(static_cast<size_t>(interpolator.Size()), samplePoints.size());

    //--------------------------------

    const mk::UInt meshPointsX = 5;
    const mk::UInt meshPointsY = 5;

    const auto mesh = MakeRectangularMeshForTesting(meshPointsX,
                                                    meshPointsY,
                                                    8.0 * delta,
                                                    8.0 * delta,
                                                    mk::Projection::cartesian,
                                                    {0.5 * delta, 0.5 * delta});

    //--------------------------------

    int propertyId = 1;
    interpolator.SetData(propertyId, data);

    // Execute
    ASSERT_EQ(interpolator.Size(), numberOfPoints);

    // Test node data

    const double initialValue = -1.0e20;
    std::vector<double> interpolationNodeResult(mesh->GetNumNodes(), initialValue);
    std::vector<double> expectedNodeResult{0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0};

    interpolator.Interpolate(propertyId, *mesh, mk::Location::Nodes, interpolationNodeResult);

    const double tolerance = 1.0e-8;

    for (size_t i = 0; i < interpolationNodeResult.size(); ++i)
    {
        EXPECT_NEAR(expectedNodeResult[i], interpolationNodeResult[i], tolerance);
    }

    // Test edge data

    std::vector<double> interpolationEdgeResult(mesh->GetNumEdges(), initialValue);
    std::vector<double> expectedEdgeResult{0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           0.0, 2000.0, 4000.0, 6000.0, 8000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0};

    interpolator.Interpolate(propertyId, *mesh, mk::Location::Edges, interpolationEdgeResult);

    for (size_t i = 0; i < interpolationEdgeResult.size(); ++i)
    {
        EXPECT_NEAR(expectedEdgeResult[i], interpolationEdgeResult[i], tolerance);
    }

    // Test face data

    std::vector<double> interpolationFaceResult(mesh->GetNumFaces(), initialValue);
    std::vector<double> expectedFaceResult{1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0,
                                           1000.0, 3000.0, 5000.0, 7000.0};

    interpolator.Interpolate(propertyId, *mesh, mk::Location::Faces, interpolationFaceResult);

    for (size_t i = 0; i < interpolationFaceResult.size(); ++i)
    {
        EXPECT_NEAR(expectedFaceResult[i], interpolationFaceResult[i], tolerance);
    }
}

TEST(SampleInterpolationTests, TriangulationInterpolationWithPoints)
{
    const mk::UInt numberOfPointsX = 11;
    const mk::UInt numberOfPointsY = 11;
    const mk::UInt numberOfPoints = numberOfPointsX * numberOfPointsY;

    std::vector<double> xPoints(numberOfPoints);
    std::vector<double> yPoints(numberOfPoints);
    std::vector<double> data(numberOfPoints);

    // Generate sample data points
    double delta = 1000.0;
    double x = 0.0;
    double y = 0.0;
    size_t count = 0;

    for (size_t i = 0; i < numberOfPointsY; ++i)
    {
        x = 0.0;

        for (size_t j = 0; j < numberOfPointsX; ++j)
        {
            xPoints[count] = x;
            yPoints[count] = y;
            data[count] = x;
            x += delta;
            ++count;
        }

        y += delta;
    }

    mk::InterpolationParameters params{.absolute_search_radius = 3.0 * delta, .minimum_number_of_samples = 1};
    mk::SampleAveragingInterpolator interpolator(xPoints, yPoints, mk::Projection::cartesian, params);

    const mk::UInt meshPointsX = 3;
    const mk::UInt meshPointsY = 3;

    const auto mesh = MakeRectangularMeshForTesting(meshPointsX,
                                                    meshPointsY,
                                                    3.0 * delta,
                                                    3.0 * delta,
                                                    mk::Projection::cartesian,
                                                    {0.5 * delta, 0.5 * delta});

    ASSERT_EQ(static_cast<size_t>(interpolator.Size()), xPoints.size());

    int propertyId = 1;
    interpolator.SetData(propertyId, data);

    // Execute
    ASSERT_EQ(interpolator.Size(), numberOfPoints);

    const double initialValue = -1.0e20;
    std::vector<double> interpolationResult(mesh->GetNumNodes(), initialValue);
    std::vector<double> expectedResult{1000.0, 1000.0, 1000.0, 2000.0, 2000.0, 2000.0, 3000.0, 3000.0, 3000.0};

    interpolator.Interpolate(propertyId, *mesh, mk::Location::Nodes, interpolationResult);

    const double tolerance = 1.0e-8;

    for (size_t i = 0; i < expectedResult.size(); ++i)
    {
        EXPECT_NEAR(expectedResult[i], interpolationResult[i], tolerance);
    }
}
