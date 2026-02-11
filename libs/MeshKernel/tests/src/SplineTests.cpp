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

#include <MeshKernel/Constants.hpp>
#include <MeshKernel/Entities.hpp>
#include <MeshKernel/Exceptions.hpp>
#include <MeshKernel/LandBoundary.hpp>
#include <MeshKernel/Operations.hpp>
#include <MeshKernel/SplineAlgorithms.hpp>
#include <MeshKernel/Splines.hpp>

#include <TestUtils/Definitions.hpp>
#include <TestUtils/SplineReader.hpp>

TEST(Splines, SetSpline)
{
    // a valid spline
    std::vector<meshkernel::Point> splineNodes({{212.001953125000, 155.627197265625},
                                                {529.253906250000, 432.379974365234},
                                                {930.506469726562, 453.380187988281},
                                                {1030.506469726562, 653.380187988281}});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(splineNodes, 0, static_cast<meshkernel::UInt>(splineNodes.size()));

    ASSERT_EQ(1, splines.GetNumSplines());
    ASSERT_EQ(4, splines.m_splineNodes[0].size());
}

TEST(Splines, CubicSplineInterpolation)
{
    // One gets the edges
    std::vector<meshkernel::Point> splineNodes;
    splineNodes.push_back(meshkernel::Point{212.001953125000, 155.627197265625});
    splineNodes.push_back(meshkernel::Point{529.253906250000, 432.379974365234});
    splineNodes.push_back(meshkernel::Point{930.506469726562, 453.380187988281});

    int pointsBetweenNodes = 20;
    auto coordinatesDerivatives = meshkernel::SplineAlgorithms::SecondOrderDerivative(splineNodes, 0, static_cast<meshkernel::UInt>(splineNodes.size()) - 1);
    std::vector<meshkernel::Point> splineCoordinates;

    for (size_t n = 0; n < splineNodes.size() - 1; n++)
    {
        for (auto p = 0; p <= pointsBetweenNodes; p++)
        {
            const double pointAdimensionalCoordinate = n + double(p) / double(pointsBetweenNodes);
            auto pointCoordinate = ComputePointOnSplineAtAdimensionalDistance(splineNodes, coordinatesDerivatives, pointAdimensionalCoordinate);
            ASSERT_TRUE(pointCoordinate.IsValid());

            splineCoordinates.push_back({pointCoordinate.x, pointCoordinate.y});
        }
    }

    const double tolerance = 1e-6;
    ASSERT_NEAR(226.817168170929, splineCoordinates[1].x, tolerance);
    ASSERT_NEAR(241.648133331299, splineCoordinates[2].x, tolerance);
    ASSERT_NEAR(256.510598720551, splineCoordinates[3].x, tolerance);
    ASSERT_NEAR(271.420314453125, splineCoordinates[4].x, tolerance);
    ASSERT_NEAR(286.393030643463, splineCoordinates[5].x, tolerance);
    ASSERT_NEAR(930.506469726562, splineCoordinates.back().x, tolerance);

    ASSERT_NEAR(172.653750896454, splineCoordinates[1].y, tolerance);
    ASSERT_NEAR(189.632350921631, splineCoordinates[2].y, tolerance);
    ASSERT_NEAR(206.515043735504, splineCoordinates[3].y, tolerance);
    ASSERT_NEAR(223.253875732422, splineCoordinates[4].y, tolerance);
    ASSERT_NEAR(453.380187988281, splineCoordinates.back().y, tolerance);
}

TEST(Splines, SplineIntersection)
{
    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{152.001571655273, 86.6264953613281});
    firstSpline.push_back(meshkernel::Point{374.752960205078, 336.378997802734});
    firstSpline.push_back(meshkernel::Point{850.255920410156, 499.130676269531});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{72.5010681152344, 391.129577636719});
    secondSpline.push_back(meshkernel::Point{462.503479003906, 90.3765411376953});
    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    double crossProductIntersection;
    meshkernel::Point dimensionalIntersection;
    double intersectionAngle;
    double firstSplineRatio;
    double secondSplineRatio;

    splines.GetSplinesIntersection(0, 1, crossProductIntersection, intersectionAngle, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    const double tolerance = 1e-6;
    ASSERT_NEAR(261.736770097059, dimensionalIntersection.x, tolerance);
    ASSERT_NEAR(245.199166962145, dimensionalIntersection.y, tolerance);
    ASSERT_NEAR(0.601498208554790, firstSplineRatio, tolerance);
    ASSERT_NEAR(0.485216749175026, secondSplineRatio, tolerance);
    ASSERT_NEAR(-0.996215079635043, crossProductIntersection, tolerance);
}

TEST(Splines, SnapToLandBoundaryTest)
{
    // Test the algorithm for snapping splines to land boundaries.

    constexpr double tolerance = 1.0e-8;

    constexpr int iterationCountForTest = 5;

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{257.002197, 442.130066},
                                                      {518.753845, 301.128662},
                                                      {938.006470, 416.629822}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> splinePoints{{281.0023, 447.3801},
                                                {367.2529, 401.6296},
                                                {461.7534, 354.3792},
                                                {517.2538, 318.3788},
                                                {614.0045, 338.629},
                                                {720.5051, 377.6294},
                                                {827.7558, 417.3798},
                                                {923.7563, 424.1299}};

    // The expected spline values after snapping to land boundary.
    std::vector<meshkernel::Point> expectedSplinePoints{{273.5868719643935, 434.2730022174478},
                                                        {359.5998304717778, 386.1712239047134},
                                                        {451.5303458337523, 338.3551703843473},
                                                        {517.7962262926076, 306.3259738916997},
                                                        {616.7325138813335, 327.9627689164845},
                                                        {725.7358644094627, 358.0902879743862},
                                                        {836.262785315633, 388.6415116416172},
                                                        {923.500177844106, 412.5818685325169}};

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(splinePoints, 0, static_cast<meshkernel::UInt>(splinePoints.size()));

    // Snap the spline to the land boundary
    splines.SnapSpline(0, landBoundary, iterationCountForTest);

    // The number of splines should be unchanged.
    ASSERT_EQ(splines.GetNumSplines(), 1) << ", expected number of splines to be unchenged, there should be only 1 spline";

    // The number of points in the spline should be unchanged
    ASSERT_EQ(splines.m_splineNodes[0].size(), expectedSplinePoints.size()) << ", expected the number of points to be unchanged.";

    // Now test the snapped spline points are close to the expected values
    for (size_t i = 0; i < splinePoints.size(); ++i)
    {
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].x, splines.m_splineNodes[0][i].x, tolerance))
            << "Expected x-value: " << expectedSplinePoints[i].x << ", actual: " << splines.m_splineNodes[0][i].x << ", relative tolerance: " << tolerance;
        EXPECT_TRUE(meshkernel::IsEqual(expectedSplinePoints[i].y, splines.m_splineNodes[0][i].y, tolerance))
            << "Expected y-value: " << expectedSplinePoints[i].y << ", actual: " << splines.m_splineNodes[0][i].y << ", relative tolerance: " << tolerance;
    }
}

TEST(Splines, SplineSnappingIndexOutOfRangeTest)
{
    // Test the algorithm for snapping splines to land boundaries.

    // The land boundary to which the spline is to be snapped.
    std::vector<meshkernel::Point> landBoundaryPoints{{257.002197, 442.130066},
                                                      {518.753845, 301.128662},
                                                      {938.006470, 416.629822}};

    meshkernel::LandBoundary landBoundary(landBoundaryPoints);

    // The original spline points.
    std::vector<meshkernel::Point> firstSpline{{281.0023, 447.3801},
                                               {367.2529, 401.6296},
                                               {461.7534, 354.3792},
                                               {517.2538, 318.3788},
                                               {614.0045, 338.629},
                                               {720.5051, 377.6294},
                                               {827.7558, 417.3798},
                                               {923.7563, 424.1299}};

    std::vector<meshkernel::Point> secondSpline({{212.001953125000, 155.627197265625},
                                                 {529.253906250000, 432.379974365234},
                                                 {930.506469726562, 453.380187988281},
                                                 {1030.506469726562, 653.380187988281}});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);
    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));
    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    // Spline 2 is out of range, so should throw an exception.
    EXPECT_THROW(splines.SnapSpline(2, landBoundary), meshkernel::ConstraintError);
}

TEST(Splines, SplineMultiIntersection)
{
    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{0.0, 0.0});
    firstSpline.push_back(meshkernel::Point{1.0, 2.0});
    firstSpline.push_back(meshkernel::Point{3.0, 3.0});

    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{0.25, 1.0});
    secondSpline.push_back(meshkernel::Point{1.5, 1.5});
    secondSpline.push_back(meshkernel::Point{2.5, 3.0});
    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    double crossProductIntersection;
    meshkernel::Point dimensionalIntersection;
    double intersectionAngle;
    double firstSplineRatio;
    double secondSplineRatio;

    bool crossing = splines.GetSplinesIntersection(0, 1, crossProductIntersection, intersectionAngle, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    ASSERT_TRUE(crossing);
}

TEST(Splines, MultiSplineConstructionTest)
{

    namespace mk = meshkernel;
    auto splines = std::make_shared<mk::Splines>(LoadSplines(TEST_FOLDER + "/data/CurvilinearGrids/seventy_splines.spl"));

    // Concatenate the points of all splines into single array separated by a separator (invalid) point
    std::vector<mk::Point> allPoints;

    for (size_t i = 0; i < splines->GetNumSplines(); ++i)
    {
        allPoints.insert (allPoints.end (), splines->m_splineNodes [i].begin (), splines->m_splineNodes [i].end ());

        if (i != splines->GetNumSplines() - 1)
        {
            allPoints.push_back (mk::Point (mk::constants::missing::doubleValue, mk::constants::missing::doubleValue));
        }

    }

    // Construct the other splines
    mk::Splines otherSplines (allPoints, splines->m_projection);

    ASSERT_EQ (splines->GetNumSplines (), otherSplines.GetNumSplines ());

    for (size_t i = 0; i < splines->GetNumSplines (); ++i)
    {
        const std::vector<mk::Point>& spline1 (splines->m_splineNodes [i]);
        const std::vector<mk::Point>& spline2 (otherSplines.m_splineNodes [i]);

        ASSERT_EQ (spline1.size (), spline2.size ());

        for (size_t i = 0; i < spline1.size (); ++i)
        {
            EXPECT_EQ (spline1[i].x, spline2[i].x);
            EXPECT_EQ (spline1[i].y, spline2[i].y);
        }

    }

}

TEST(Splines, SimpleSplineIntersectionAngle)
{
    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{0.0, 0.0});
    firstSpline.push_back(meshkernel::Point{10.0, 10.0});

    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{0.0, 3.0});
    secondSpline.push_back(meshkernel::Point{10.0, 3.0});

    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    std::vector<meshkernel::Point> thirdSpline;
    thirdSpline.push_back(meshkernel::Point{0.0, 8.0});
    thirdSpline.push_back(meshkernel::Point{10.0, 8.0});

    splines.AddSpline(thirdSpline, 0, static_cast<meshkernel::UInt>(thirdSpline.size()));

    double crossProductIntersection;
    meshkernel::Point dimensionalIntersection;
    double intersectionAngle;
    double firstSplineRatio;
    double secondSplineRatio;

    constexpr double tolerance = 1.0e-10;
    constexpr double expectedSpline0Spline1IntersectionAngle = 45.0;
    constexpr double expectedSpline0Spline2IntersectionAngle = 45.0;
    constexpr double expectedSpline1Spline2IntersectionAngle = meshkernel::constants::missing::doubleValue;

    bool crossing = splines.GetSplinesIntersection(0, 1, crossProductIntersection, intersectionAngle, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    EXPECT_TRUE (crossing);
    EXPECT_NEAR (intersectionAngle, expectedSpline0Spline1IntersectionAngle, tolerance);

    crossing = splines.GetSplinesIntersection(0, 2, crossProductIntersection, intersectionAngle, dimensionalIntersection, firstSplineRatio, secondSplineRatio);

    EXPECT_TRUE (crossing);
    EXPECT_NEAR (intersectionAngle, expectedSpline0Spline2IntersectionAngle, tolerance);

    crossing = splines.GetSplinesIntersection(1, 2, crossProductIntersection, intersectionAngle, dimensionalIntersection, firstSplineRatio, secondSplineRatio);
    EXPECT_FALSE (crossing);
    EXPECT_NEAR (intersectionAngle, expectedSpline1Spline2IntersectionAngle, tolerance);
}


TEST(Splines, MultiSegmentSplineIntersectionAngle)
{
    meshkernel::Splines splines(meshkernel::Projection::cartesian);

    std::vector<meshkernel::Point> firstSpline;
    firstSpline.push_back(meshkernel::Point{-100.0, 7.5});
    firstSpline.push_back(meshkernel::Point{20.0, 7.5});

    splines.AddSpline(firstSpline, 0, static_cast<meshkernel::UInt>(firstSpline.size()));

    std::vector<meshkernel::Point> secondSpline;
    secondSpline.push_back(meshkernel::Point{0.0, 0.0});
    secondSpline.push_back(meshkernel::Point{5.0, 5.0});
    secondSpline.push_back(meshkernel::Point{5.0, 10.0});
    secondSpline.push_back(meshkernel::Point{0.0, 15.0});

    splines.AddSpline(secondSpline, 0, static_cast<meshkernel::UInt>(secondSpline.size()));

    double crossProductIntersection;
    meshkernel::Point intersectionPoint;
    double intersectionAngle;
    double firstSplineRatio;
    double secondSplineRatio;

    constexpr double tolerance = 1.0e-10;
    constexpr double expectedSplineIntersectionAngle = 90.0;

    bool crossing = splines.GetSplinesIntersection(0, 1, crossProductIntersection, intersectionAngle, intersectionPoint, firstSplineRatio, secondSplineRatio);

    constexpr double expectedIntersectionPointX = 5.7470703125;
    // Should be halfway between the second and third spline support points
    constexpr double expectedIntersectionPointY = 7.5;
    constexpr double expectedFirstSplineRatio = 0.8812255859375;
    // Should be half way along the second spline segment
    constexpr double expectedSecondSplineRatio = 1.5;

    EXPECT_TRUE (crossing);
    EXPECT_NEAR (intersectionAngle, expectedSplineIntersectionAngle, tolerance);
    EXPECT_NEAR (intersectionPoint.x, expectedIntersectionPointX, tolerance);
    EXPECT_NEAR (intersectionPoint.y, expectedIntersectionPointY, tolerance);
    EXPECT_NEAR (firstSplineRatio, expectedFirstSplineRatio, tolerance);
    EXPECT_NEAR (secondSplineRatio, expectedSecondSplineRatio, tolerance);
}


TEST(Splines, MultiSplineIntersectionAngle)
{
    namespace mk = meshkernel;
    auto splines = std::make_shared<mk::Splines>(LoadSplines(TEST_FOLDER + "/data/CurvilinearGrids/seventy_splines.spl"));

    std::vector<mk::Point> mySpline{{30825.9994380052, 377851.802873428},
                                    {31503.2022933559, 379206.208584129},
                                    {33217.3720209622, 381915.020005531},
                                    {33153.884253273, 383481.05160853},
                                    {32053.4296133283, 384454.530713096},
                                    {32730.6324686789, 385576.147942271},
                                    {34085.0381793802, 385999.399726865},
                                    {34338.9892501366, 387438.455794485}};

    std::vector<int> splineIndices;
    std::vector<double> angles;
    std::vector<double> xCrossOver;
    std::vector<double> yCrossOver;

    std::vector<int> expectedIntersectedSplines{0, 1, 2, 3, 4, 10, 23};

    std::vector<double> expectedIntersectedAngles{9.369091291645886e+01, 1.128522130510987e+02, 1.133073718558457e+02, 1.146340134566103e+02,
                                                  4.215466692286142e+01, 6.936103716736207e+01, 8.094729661593756e+01};

    std::vector<double> expectedIntersectedCoordX{3.345969673585541e+04,3.291179368046583e+04,3.199660809941513e+04,3.204744958459182e+04,
                                                  3.270633719945683e+04,3.431874691644472e+04,3.108691001430429e+04};

    std::vector<double> expectedIntersectedCoordY{3.826926563132096e+05,3.813522268790855e+05,3.799435990928682e+05,3.848083038591905e+05,
                                                  3.838236381325046e+05,3.865335746218469e+05,3.784979991789844e+05};

    splines->GetAllIntersections (mySpline, splineIndices, angles, xCrossOver, yCrossOver);

    ASSERT_EQ (splineIndices.size (), expectedIntersectedSplines.size ());

    constexpr double tolerance = 1.0e-5;

    for (size_t i = 0; i < splineIndices.size(); ++i)
    {
        EXPECT_EQ(splineIndices[i], expectedIntersectedSplines[i]);
        EXPECT_NEAR(angles[i], expectedIntersectedAngles[i], tolerance);
        EXPECT_NEAR(xCrossOver[i], expectedIntersectedCoordX[i], tolerance);
        EXPECT_NEAR(yCrossOver[i], expectedIntersectedCoordY[i], tolerance);
    }

}
